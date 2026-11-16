"""Reconstruct saved FSL-MRS molarity with NumPy equations, without quantifying."""

import argparse
import csv
import json
from pathlib import Path

import fsl_mrs
import numpy as np
from fsl_mrs.utils.mrs_io.main import read_basis
from fsl_mrs.utils.constants import H2O_MOLARITY, H2O_PROTONS, TISSUE_WATER_DENSITY

from simulation.frozen_basis_normalization import basis_identity
from simulation.mrscloud_quantitative_pair import REFERENCE_BASIS


def spectral_integral(fid, ppm, limits):
    """Independent reproduction of FSL's first-point/FFT/integration convention."""
    data = np.asarray(fid).copy()
    data[0] *= .5
    spectrum = np.fft.fftshift(np.fft.fft(data, norm='ortho'))
    lo, hi = sorted(np.argmin(abs(ppm-limit)) for limit in limits)
    return float(np.trapezoid(abs(spectrum[lo:hi+1].real)))


def relaxation(te, tr, t1, t2, flip_angle_deg=90):
    saturation = np.exp(-tr/t1)
    return float((1-saturation)/(1-np.cos(np.deg2rad(flip_angle_deg))*saturation)*np.exp(-te/t2))


def derive_case(case_dir):
    case_dir = Path(case_dir)
    fit = json.loads((case_dir/'fitting_results.json').read_text())
    config = json.loads((case_dir/'simulation_config.json').read_text())
    truth = json.loads((case_dir/'ground_truth.json').read_text())
    frozen_meta = config['synthetic_calibration']['frozen_basis']
    if basis_identity(REFERENCE_BASIS) != frozen_meta['reference_identity']:
        raise ValueError('saved reference basis identity mismatch')
    basis = read_basis(str(REFERENCE_BASIS))
    raw = basis.get_formatted_basis(config['bandwidth_hz'],config['points'])
    k = float(np.mean(abs(raw)))
    if k != frozen_meta['normalization_scalar']:
        raise ValueError('saved normalization mismatch')
    b = raw[:,basis.names.index('NAA')]/k
    q = float(abs(raw[0,basis.names.index('Ref0ppm')])/k)
    t = np.arange(config['points'])*config['dwell_time_s']
    ppm = np.fft.fftshift(np.fft.fftfreq(config['points'],config['dwell_time_s']))/config['spectrometer_frequency_mhz']+config['frequency_origin_ppm']
    params = fit['fit_parameters']
    scale_d, scale_b = fit['scaling']['FID'], fit['scaling']['basis']
    coefficient_mol_l = params['NAA']*scale_b/scale_d*k
    # RefIntegral(no_phase=True) removes Phi0/Phi1, retains fitted gamma/eps.
    metabolite_shape = b*np.exp(-(params['gamma_0']+1j*params['eps_0'])*t)
    _, gamma_w, sigma_w, _, _ = fit['water_optimizer']['params']
    water_shape = np.exp(-gamma_w*t-sigma_w*t*t).astype(complex)
    fitted_water_amp = fit['water_optimizer']['params'][0]/fit['water_internal_fit_scaling']/scale_d
    j_m = spectral_integral(metabolite_shape,ppm,(1.8,2.2))
    j_w = spectral_integral(water_shape,ppm,(1.65,7.65))
    i_m = scale_d*coefficient_mol_l*j_m
    i_w = scale_d*fitted_water_amp*j_w
    quant = fit['quantification']
    fractions = {name:quant[name] for name in ['WM','GM','CSF']}
    rw = {name:relaxation(quant['TE_s'],quant['TR_s'],quant['T1_s']['H2O_'+name],quant['T2_s']['H2O_'+name]) for name in fractions}
    rm = relaxation(quant['TE_s'],quant['TR_s'],quant['T1_s']['METAB'],quant['T2_s']['METAB'])
    density = sum(fractions[name]*TISSUE_WATER_DENSITY[name] for name in fractions)
    water_relaxed_mM = float(H2O_MOLARITY*sum(fractions[name]*TISSUE_WATER_DENSITY[name]*rw[name] for name in fractions))
    csf = 1/(1-fractions['CSF'])
    additional = quant['additional_factor']
    n_m = quant['reference_protons']
    q_eff = water_relaxed_mM/rm*csf*additional
    prediction = i_m/i_w*(H2O_PROTONS/n_m)*q_eff
    expected_mM = truth['metabolite_concentrations_mM']['NAA']
    c_w = truth['water_concentration_mol_l']
    g = config['receiver_factor']
    factors = {
        'coefficient_recovery': coefficient_mol_l/(expected_mM*.001),
        'water_amplitude_recovery': 2*q*g*c_w/fitted_water_amp,
        'reference_integral_vs_Ref0ppm': j_m/(n_m*q*j_w),
        'tissue_water_vs_40M': H2O_MOLARITY*density/(1000*c_w),
        'water_to_metabolite_relaxation': water_relaxed_mM/(H2O_MOLARITY*density)/rm,
        'CSF_correction': csf, 'additional_correction': additional,
    }
    # Proposed compatibility constants use only the frozen basis, acquisition
    # parameters and declared software model, not fitted concentrations/ratio.
    j_m_ideal = spectral_integral(b,ppm,(1.8,2.2))
    j_w_ideal = spectral_integral(np.ones(config['points'],dtype=complex),ppm,(1.65,7.65))
    response_ideal = j_m_ideal/(n_m*q*j_w_ideal)
    required_water_default = g*(j_m_ideal/j_w_ideal)*(H2O_PROTONS/n_m)*(q_eff/1000)
    required_water_no_relaxation_40M = g*(j_m_ideal/j_w_ideal)*(H2O_PROTONS/n_m)*c_w
    config_only_additional_scale = (1000*c_w/q_eff)/response_ideal
    ideal_prediction = expected_mM*response_ideal*q_eff/(1000*c_w)
    recorded = fit['comparison']['fsl_molarity_mM']
    return {'ground_truth_mM':expected_mM,'recorded_fsl_mM':recorded,'reconstructed_fsl_mM':prediction,
            'absolute_error_mM':prediction-recorded,'relative_error':prediction/recorded-1,
            'ideal_prediction_mM':ideal_prediction,'ideal_absolute_error_mM':ideal_prediction-recorded,
            'factors':factors,'factor_product':float(np.prod(list(factors.values()))),
            'J_metabolite':j_m,'J_water':j_w,'I_metabolite':i_m,'I_water':i_w,
            'saved_integral_differences':{'metabolite':i_m-fit['metabolite_integral'],'water':i_w-fit['water_integral']},
            'basis_identity':frozen_meta['reference_identity'],'frozen_scalar':k,'q':q,
            'fitted_water_amplitude':fitted_water_amp,'coefficient_mol_l':coefficient_mol_l,
            'relaxation_water':rw,'relaxation_metabolite':rm,
            'water_relaxed_mM':water_relaxed_mM,'effective_water_mM':q_eff,
            'tissue_water_density':density,'constants':{'H2O_MOLARITY_mM':H2O_MOLARITY,'water_protons':H2O_PROTONS,'reference_protons':n_m},
            'ideal_response_factor':response_ideal,
            'J_metabolite_ideal':j_m_ideal,'J_water_ideal':j_w_ideal,
            'proposed_water_amplitude_default_FSL':required_water_default,
            'proposed_water_amplitude_no_relaxation_40M':required_water_no_relaxation_40M,
            'proposed_configuration_additional_scale':config_only_additional_scale,
            'configuration_only_predicted_mM':prediction*config_only_additional_scale,
            'proposal_applied':False}


def run_derivation(source_dir, output_dir):
    source, out = Path(source_dir), Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    rows=[derive_case(source/f'NAA_{c:02d}mM') for c in [5,10,15,20]]
    result={'FSL_MRS_version':fsl_mrs.__version__,'cases':rows}
    (out/'derivation.json').write_text(json.dumps(result,indent=2)+'\n')
    columns=['ground_truth_mM','recorded_fsl_mM','reconstructed_fsl_mM','absolute_error_mM','relative_error','ideal_prediction_mM','ideal_absolute_error_mM','configuration_only_predicted_mM']
    with (out/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        writer.writerows({key:row[key] for key in columns} for row in rows)
    with (out/'factors_10mM.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['factor','value','cumulative_predicted_mM'])
        cumulative=10.
        for name,value in rows[1]['factors'].items():
            cumulative*=value;writer.writerow([name,value,cumulative])
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir')
    parser.add_argument('--source-dir',default='outputs/fsl_mrs_water_scaling_phase5b')
    args=parser.parse_args()
    run_derivation(args.source_dir,args.output_dir)
