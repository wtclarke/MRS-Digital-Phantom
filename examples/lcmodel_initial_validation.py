"""Phase 7b: exact Phase 7a conversion and one local, noise-free LCModel fit."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np
from fsl_mrs.utils import mrs_io

from utils.lcmodel_quantitative_io import (
    NAMES, compatible_basis, control_text, independent_water_scaling,
    parse_basis, parse_coord, parse_raw, parse_table, write_raw,
)

SOURCE = Path('outputs/mrscloud_noise_fitting_validation/full/SNR_infinite/realization_00')
BUNDLE = Path('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative')
OUTPUT = Path('outputs/lcmodel_validation')
BINARY = Path(os.environ.get('LCMODEL_BINARY', 'lcmodel'))


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def array_errors(actual, expected):
    return dict(max_absolute_complex_error=float(np.max(abs(actual-expected))),
        relative_complex_error=float(np.linalg.norm(actual-expected)/np.linalg.norm(expected)),
        initial_amplitude_error=float(abs(actual[0])-abs(expected[0])),
        max_phase_error_rad=float(np.max(abs(np.angle(actual*np.conj(expected))))),
        sample_order_preserved=bool(np.array_equal(actual, expected)))


def run(output=OUTPUT, *, reuse_fit=False):
    if not BINARY.is_file() or not os.access(BINARY, os.X_OK):
        raise RuntimeError(f'LCModel executable unavailable: {BINARY}; set LCMODEL_BINARY to your licensed executable')
    output=Path(output).resolve(); directory=output/'noise_free'
    directory.mkdir(parents=True, exist_ok=True)
    (directory/'lcmodel_outputs').mkdir(exist_ok=True)
    protected=list(SOURCE.iterdir()) + list(BUNDLE.iterdir()) + [BINARY]
    protected=[p for p in protected if p.is_file()]
    before={str(p.resolve()):digest(p) for p in protected}
    config=json.loads((SOURCE/'simulation_config.json').read_text())
    truth=json.loads((SOURCE/'ground_truth.json').read_text())
    fsl=json.loads((SOURCE/'fit_result.json').read_text())
    normalization=json.loads((BUNDLE/'normalization.json').read_text())
    if normalization!=config['synthetic_calibration']['normalization']:
        raise ValueError('source acquisition and basis calibration metadata differ')
    for filename,expected in fsl['input_hashes'].items():
        if digest(SOURCE/filename)!=expected:
            raise ValueError('source acquisition changed since Phase 7a fitting')
    basis=BUNDLE/'metabolite_basis.BASIS'
    entries, basis_meta=parse_basis(basis)
    if config['noise'] or config['points']!=2048 or config['dwell_time_s']!=.0005:
        raise ValueError('expected exact Phase 7a noise-free 2048-point acquisition')
    if basis_meta!=dict(points=8192, dwell_time_s=.00025, frequency_mhz=127.731):
        raise ValueError('incompatible Phase 6a basis sampling')
    water_concentration=truth['effective_voxel_water_concentration_mol_l']
    if water_concentration!=39.0165:
        raise ValueError('expected the existing Phase 7a WM water concentration')
    validation={}
    for name, target in [('metabolites','metabolites.RAW'),('water_reference','water_reference.H2O')]:
        reader=mrs_io.read_FID(str(SOURCE/f'{name}.nii.gz')); mrs=reader.mrs()
        fid=mrs.FID
        if len(fid)!=2048 or not np.isclose(mrs.dwellTime,.0005) or not np.isclose(mrs.centralFrequency/1e6,127.731):
            raise ValueError('source NIfTI-MRS grid mismatch')
        path=directory/target
        if not path.exists():write_raw(path,fid,frequency_mhz=config['spectrometer_frequency_mhz'])
        stored,header=parse_raw(path)
        checks=array_errors(stored.conj(),fid)
        if checks['max_absolute_complex_error']!=0:raise ValueError('RAW does not preserve source FID')
        if header['frequency_mhz']!=127.731:raise ValueError('RAW frequency metadata mismatch')
        checks['initial_complex_FSL']=[float(fid[0].real),float(fid[0].imag)]
        checks['initial_complex_RAW']=[float(stored[0].real),float(stored[0].imag)]
        validation[name]=checks
    water_raw,_=parse_raw(directory/'water_reference.H2O')
    met_raw,_=parse_raw(directory/'metabolites.RAW')
    source_met=mrs_io.read_FID(str(SOURCE/'metabolites.nii.gz')).mrs().FID
    source_water=mrs_io.read_FID(str(SOURCE/'water_reference.nii.gz')).mrs().FID
    validation['initial_complex_ratio_error']=float(abs((met_raw[0]/water_raw[0]).conjugate()-source_met[0]/source_water[0]))
    validation['water_FSL_frequency_hz']=float(np.angle(source_water[1]/source_water[0])/(2*np.pi*.0005))
    validation['water_RAW_frequency_hz']=float(np.angle(water_raw[1]/water_raw[0])/(2*np.pi*.0005))
    validation['water_ppm']=4.68-validation['water_RAW_frequency_hz']/127.731
    adapted=directory/'lcmodel_compatible.BASIS'
    if not adapted.exists():
        save(directory/'basis_compatibility.json',compatible_basis(basis,adapted,data_origin_ppm=4.68))
    adapted_entries,_=parse_basis(adapted)
    if [e[0] for e in entries]!=[e[0] for e in adapted_entries]:raise ValueError('basis names changed')
    for original, copy in zip(entries,adapted_entries):
        if not np.array_equal(original[2],copy[2]):raise ValueError('basis amplitudes changed')
        if int(re.search(r'ISHIFT\s*=\s*([-+]?\d+)',copy[1])[1])!=8:raise ValueError('basis referencing mismatch')
    if "FMTBAS='(2E25.16)'" not in adapted.read_text():raise ValueError('missing LCModel fixed-width format')
    text=control_text(directory,adapted,names=[e[0] for e in entries],points=2048,
        dwell_time_s=.0005,frequency_mhz=127.731,origin_ppm=4.68,water_mol_l=water_concentration)
    control=directory/'lcmodel.CONTROL'
    if reuse_fit:
        if control.read_text()!=text:raise ValueError('cached control differs from declared configuration')
        execution=json.loads((directory/'execution.json').read_text())
    else:
        if (directory/'execution.json').exists():raise FileExistsError('fit already exists; use --reuse-fit')
        control.write_text(text)
        process=subprocess.run([str(BINARY)],input=text,text=True,cwd=directory,capture_output=True,timeout=300)
        execution=dict(returncode=process.returncode,stdout=process.stdout,stderr=process.stderr)
        save(directory/'execution.json',execution)
    fits=directory/'lcmodel_outputs'
    result=parse_table(fits/'fit.TABLE')
    if execution['returncode']!=0 or result['fatal'] or set(result['concentrations'])!=set(NAMES):
        raise RuntimeError(f'LCModel fit failed despite possible zero exit: {result}')
    for suffix in ['TABLE','COORD','PRINT','PS']:
        if (fits/f'fit.{suffix}').stat().st_size==0:raise RuntimeError('empty LCModel output')
    coord=parse_coord(fits/'fit.COORD')
    result['relative_phased_real_spectrum_residual']=float(np.linalg.norm(coord['data']-coord['fit'])/np.linalg.norm(coord['data']))
    result['complex_residual_available']=False
    printed=(fits/'fit.PRINT').read_text()
    result['printed_scaling']={key:float(re.search(pattern,printed)[1]) for key,pattern in {
        'normalized_basis_area':r'Normalized area of reference Basis singlet =\s*([\d.E+-]+)',
        'water_area':r'Area of unsuppressed water peak =\s*([\d.E+-]+)',
        'FCALIB':r'FCALIB =\s*([\d.E+-]+)'}.items()}
    result['independent_scaling']=independent_water_scaling(basis,water_raw)
    result['source_calibration']=config['synthetic_calibration']
    result['water_mM']=water_concentration*1000
    result['execution']=execution
    result['protected_input_sha256']=before
    result['validation']=validation
    rows=[]
    for name in NAMES:
        ground=fsl['coefficients'][name]['truth_mM']; concentration=result['concentrations'][name]['concentration']
        rows.append(dict(metabolite=name,ground_truth_mM=ground,
            FSL_recovered_source_coefficient_mM=fsl['coefficients'][name]['estimate_mM'],
            FSL_standard_mM=fsl['coefficients'][name]['standard_FSL_mM'],LCModel_mM=concentration,
            LCModel_minus_truth_mM=concentration-ground,LCModel_relative_error=(concentration-ground)/ground,
            LCModel_percent_SD=result['concentrations'][name]['percent_SD'],
            normalization_only_prediction_mM=ground*result['independent_scaling']['ideal_coefficient_mM_to_report_factor']))
    with (output/'comparison.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    command=f'cd {shlex.quote(str(directory))} && {shlex.quote(str(BINARY))} < lcmodel.CONTROL'
    log=command+'\n'+json.dumps(execution,indent=2)+'\n'
    if reuse_fit:log+='Cached fit reused; original first-run stdout/stderr and exit code retained.\n'
    (output/'execution_log.txt').write_text(log)
    save(directory/'conversion_validation.json',validation)
    for device,suffix,extra in [('pdfwrite','pdf',[]),('png16m','png',['-r120','-dFirstPage=1','-dLastPage=1'])]:
        path=fits/f'fit.{suffix}'
        if not path.exists():
            subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE',f'-sDEVICE={device}',*extra,
                            f'-sOutputFile={path}',str(fits/'fit.PS')],check=True,capture_output=True)
    after={str(p.resolve()):digest(p) for p in protected}
    if before!=after:raise RuntimeError('protected source changed')
    result['source_files_unchanged']=True
    save(directory/'fitting_results.json',result)
    print(json.dumps(dict(version=result['version'],conversion=validation,comparison=rows,
                         residual=result['relative_phased_real_spectrum_residual'],
                         independent_scaling=result['independent_scaling'],diagnostics=result['diagnostics']),indent=2))
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUTPUT)
    parser.add_argument('--reuse-fit',action='store_true',help='Validate existing first-fit artifacts without executing LCModel again')
    args=parser.parse_args();run(args.output,reuse_fit=args.reuse_fit)
