"""Opt-in deterministic voxel mixing on the frozen Phase 6a response scale."""

from numbers import Integral

import numpy as np

from simulation.mrscloud_full_basis import synthesize
from simulation.tissue_relaxation import relaxation_metadata
from simulation.mrscloud_quantitative_noise import add_pair_noise

TISSUES=('WM','GM','CSF')


def simulate_voxel(tissue_fractions, *, calibration_bundle, metabolite_concentrations_mM,
                   water_concentrations_mol_l, receiver_factor=1., normalize_fractions=False,
                   bandwidth_hz=2000, points=2048, expected_calibration=None, relaxation=None, noise=None):
    """Mix explicit tissue concentrations; fractions are volume fractions.

    Tissue profiles must specify all three tissues and identical metabolite
    names (including explicit zeros). Fractions may omit zero-volume tissues.
    Tissue-only concentrations mean the WM/GM average, excluding CSF signal.
    Optional relaxation uses explicit second-valued T1/T2 profiles and the
    actual sequence TE; None preserves the original no-attenuation pathway.
    """
    if set(tissue_fractions)-set(TISSUES):raise ValueError('only WM, GM and CSF fractions are supported')
    original={t:float(tissue_fractions.get(t,0.)) for t in TISSUES}
    total=sum(original.values())
    if not all(np.isfinite(v) and v>=0 for v in original.values()) or not np.isfinite(total) or total<=0:
        raise ValueError('fractions must be finite, nonnegative and have a positive sum')
    if not normalize_fractions and not np.isclose(total,1.,rtol=0,atol=1e-8):
        raise ValueError('fractions must sum to one; use normalize_fractions=True explicitly for weights')
    fractions={t:v/total for t,v in original.items()}
    if set(metabolite_concentrations_mM)!=set(TISSUES) or set(water_concentrations_mol_l)!=set(TISSUES):
        raise ValueError('explicit concentration profiles for WM, GM and CSF are required')
    names=list(metabolite_concentrations_mM['WM'])
    if not names or any(set(metabolite_concentrations_mM[t])!=set(names) for t in TISSUES):
        raise ValueError('all tissue profiles must contain the same nonempty metabolite set; specify zeros explicitly')
    concentrations={t:{m:float(metabolite_concentrations_mM[t][m]) for m in names} for t in TISSUES}
    water={t:float(water_concentrations_mol_l[t]) for t in TISSUES}
    if not all(np.isfinite(v) and v>=0 for t in TISSUES for v in [water[t],*concentrations[t].values()]):
        raise ValueError('tissue concentrations must be finite and nonnegative')
    average={m:sum(fractions[t]*concentrations[t][m] for t in TISSUES) for m in names}
    effective_water=sum(fractions[t]*water[t] for t in TISSUES)
    relax=relaxation_metadata(relaxation,names) if relaxation is not None else None
    signal_average=average if relax is None else {
        m:sum(fractions[t]*concentrations[t][m]*relax['metabolite_factors'][t][m] for t in TISSUES)
        for m in names}
    signal_water=effective_water if relax is None else sum(
        fractions[t]*water[t]*relax['water_factors'][t] for t in TISSUES)
    pair=synthesize(calibration_bundle,signal_average,water_concentration_mol_l=signal_water,
        receiver_factor=receiver_factor,bandwidth_hz=bandwidth_hz,points=points,
        expected_calibration=expected_calibration)
    parenchyma=fractions['WM']+fractions['GM']
    tissue_only={m:sum(fractions[t]*concentrations[t][m] for t in ['WM','GM'])/parenchyma
                 if parenchyma>0 else None for m in names}
    calibration=pair['simulation_config']['synthetic_calibration']
    if relax is not None and not np.isclose(relax['TE_s'],calibration['sequence']['TE_ms']/1000,rtol=0,atol=1e-12):
        raise ValueError('relaxation TE_s must match the frozen sequence TE; changing TE requires a new basis')
    pair['ground_truth'].update(
        metabolite_concentrations_mM=average,water_concentration_mol_l=effective_water,
        concentration_convention='tissue profiles: nominal concentration per total volume of that tissue compartment; voxel averages: per full selected voxel volume; metabolites mM, water mol/L',
        tissue_fractions=fractions,input_tissue_fractions=original,
        fraction_normalization={'input_sum':total,'explicitly_requested':bool(normalize_fractions),
                                'applied_divisor':total,'default_sum_tolerance':1e-8},
        tissue_specific_metabolite_concentrations_mM=concentrations,
        tissue_specific_water_concentrations_mol_l=water,
        voxel_average_metabolite_concentrations_mM=average,
        tissue_only_metabolite_concentrations_mM=tissue_only,
        tissue_only_definition='sum(WM/GM fraction * concentration) / (f_WM+f_GM); excludes CSF contribution; null for pure CSF',
        effective_voxel_water_concentration_mol_l=effective_water,
        calibration_identity=calibration['calibration_identity'],
        shared_normalization_scalar=calibration['normalization']['scalar'])
    pair['simulation_config'].update(
        pathway='synthetic_mrscloud_quantitative_single_voxel',tissue_fractions=fractions,
        tissue_mixing='volume_fraction_weighted_nominal_concentrations',
        anatomical_volume_multiplier_applied=False,tissue_specific_relaxation=False,
        B0_inhomogeneity=False,noise=False,spatial_simulation=False)
    if relax is not None:
        pair['ground_truth'].update(relaxation=relax,
            relaxation_weighted_metabolite_concentrations_mM=signal_average,
            relaxation_weighted_water_concentration_mol_l=signal_water)
        pair['simulation_config'].update(relaxation=relax,TE_s=relax['TE_s'],TR_s=relax['TR_s'],
            tissue_specific_relaxation=relax['mode']=='tissue_specific')
    pair['time_s']=np.arange(points)/bandwidth_hz
    return add_pair_noise(pair, **noise) if noise is not None else pair


def simulate_phantom_voxel(phantom, voi_bounds, **kwargs):
    """Read a CPU phantom label-map ROI; never call its stochastic simulator.

    Bounds use DigitalPhantom.simulate_data's exclusive-stop six-index tuple.
    Equal-volume anatomical cells produce the same one-hot mean fractions.
    """
    labels=np.asarray(phantom.get_label_data())
    if labels.ndim!=3 or not np.all(np.isfinite(labels)) or not np.all(labels==np.floor(labels)):
        raise ValueError('phantom labels must be a finite integer-coded 3D map')
    bounds=tuple(voi_bounds)
    if len(bounds)!=6 or any(not isinstance(x,Integral) or isinstance(x,bool) for x in bounds):
        raise ValueError('VOI bounds must contain six integer start/stop indices')
    starts=bounds[::2];stops=bounds[1::2]
    if any(a<0 or b<=a or b>n for a,b,n in zip(starts,stops,labels.shape)):
        raise ValueError('VOI must be nonempty and lie inside the phantom label map')
    names=list(phantom.label_names)
    if len(set(names))!=len(names) or not set(TISSUES)<=set(names):
        raise ValueError('phantom label names must identify WM, GM and CSF unambiguously')
    region=labels[tuple(slice(a,b) for a,b in zip(starts,stops))]
    ids={t:names.index(t) for t in TISSUES}
    if not np.all(np.isin(region,list(ids.values()))):
        raise ValueError('quantitative VOI contains background or an unsupported tissue; it is not renormalized away')
    fractions={t:float(np.mean(region==i)) for t,i in ids.items()}
    spacing=np.asarray(phantom.spacing,dtype=float);affine=np.asarray(phantom.affine,dtype=float)
    if (spacing.shape!=(3,) or not np.all(np.isfinite(spacing)) or np.any(spacing<=0)
            or affine.shape!=(4,4) or not np.all(np.isfinite(affine))):
        raise ValueError('phantom geometry requires positive 3D spacing and finite 4x4 affine')
    pair=simulate_voxel(fractions,**kwargs)
    pair['ground_truth']['voxel_geometry']={
        'exclusive_stop_bounds_index':[int(v) for v in bounds],
        'anatomical_spacing_mm':spacing.tolist(),'anatomical_cell_count':int(region.size),
        'voxel_volume_mm3':float(region.size*np.prod(spacing)),
        'anatomical_affine':affine.tolist()}
    pair['simulation_config']['fraction_source']='DigitalPhantom.get_label_data; equal-volume label-indicator means'
    return pair
