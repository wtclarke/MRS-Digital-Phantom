"""LCModel water scaling for saved headless acquisitions; no fitted calibration."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

import numpy as np
from fsl_mrs.utils import mrs_io

from utils.lcmodel_quantitative_io import compatible_basis, parse_basis, parse_coord, write_raw


BINARY = Path('/Users/wclarke/Documents/LCModel/test_lcm/lcmodel')


def fit_case(source: Path, basis: Path, output: Path) -> dict:
    """Return voxel-average mM estimates for every non-water BASIS entry.

    Unit water/metabolite attenuation and voxel-weighted water concentration are
    explicit protocol assumptions. Source acquisitions and BASIS are never changed.
    """
    source, basis, output = source.resolve(), basis.resolve(), output.resolve()
    cfg = json.loads((source / 'simulation_config.json').read_text())
    truth = json.loads((source / 'ground_truth.json').read_text())
    provenance = json.loads((source / 'calibration_provenance.json').read_text())
    entries, metadata = parse_basis(basis)
    names = [entry[0] for entry in entries if entry[0] not in {'H2O', 'Ref0ppm'}]
    if len(names) != len(set(names)) or 'NAA' not in names:
        raise ValueError('LCModel protocol requires unique names and an NAA scaling reference')
    if any("'" in str(path) or '\n' in str(path) for path in (source, basis, output)):
        raise ValueError('unsupported quote/newline in LCModel paths')
    if any("'" in name or '\n' in name for name in names):
        raise ValueError('unsupported BASIS metabolite name')
    if not BINARY.is_file() or not os.access(BINARY, os.X_OK):
        raise RuntimeError('licensed LCModel executable unavailable')
    protected = [source / name for name in (
        'total.nii.gz', 'water_reference.nii.gz', 'ground_truth.json',
        'simulation_config.json', 'calibration_provenance.json')] + [basis, BINARY]
    def hashes():
        return {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in protected}
    before = hashes()
    if before[str(basis)] != provenance['basis_sha256']:
        raise ValueError('BASIS differs from saved acquisition provenance')
    readers = [mrs_io.read_FID(str(source / name)) for name in ('total.nii.gz', 'water_reference.nii.gz')]
    fids = [np.asarray(reader.mrs().FID) for reader in readers]
    dt = float(readers[0].dwelltime)
    cf = float(readers[0].spectrometer_frequency[0])
    origin = float(provenance['frequency_origin_ppm'])
    te = float(cfg['basis_set_settings']['TE'])
    text = basis.read_text()
    def field(name):
        match = re.search(rf'\b{name}\s*=\s*([-+\d.eEdD]+)', text, re.I)
        if match is None:
            raise ValueError(f'BASIS requires explicit {name}')
        value = float(match[1].replace('D', 'E').replace('d', 'e'))
        if not np.isfinite(value):
            raise ValueError(f'nonfinite BASIS {name}')
        return value
    fwhmba = field('FWHMBA')
    if fwhmba <= 0 or not np.isclose(te, field('ECHOT')) or not np.isclose(cf, metadata['frequency_mhz']):
        raise ValueError('invalid linewidth or mismatched acquisition/BASIS TE/frequency')
    for reader, fid in zip(readers, fids):
        if (fid.ndim != 1 or fid.shape != fids[0].shape or not np.iscomplexobj(fid)
                or not np.all(np.isfinite(fid)) or not np.isclose(reader.dwelltime, dt, rtol=1e-7)
                or not np.isclose(reader.spectrometer_frequency[0], cf)
                or not np.isclose(reader.SpecFreqChemShift, origin)):
            raise ValueError('paired acquisition samples or metadata mismatch')
    water_mM = 1000 * float(truth['water_reference']['nominal_water_concentration_mol_l'])
    if not np.isfinite(water_mM) or water_mM <= 0:
        raise ValueError('positive voxel-weighted water concentration required')
    output.mkdir(parents=True, exist_ok=False)
    native = output / 'lcmodel_outputs'
    native.mkdir()
    adapted = output / 'lcmodel_compatible.BASIS'
    conversion = compatible_basis(basis, adapted, data_origin_ppm=origin)
    for filename, fid in zip(('metabolites.RAW', 'water_reference.H2O'), fids):
        write_raw(output / filename, fid, frequency_mhz=cf, echo_time_ms=te)
    omitted = [entry[0] for entry in entries if entry[0] not in names]
    control = f""" $LCMODL
 KEY=210387309
 TITLE='Headless paired water reference; unit attenuation'
 NUNFIL={len(fids[0])}, DELTAT={dt:.16g}, HZPPPM={cf:.16g}, ECHOT={te:.16g}
 PPMCEN={origin:.16g}, PPMH2O=4.70, PPMST=4.2, PPMEND=.2
 FILRAW='{output / 'metabolites.RAW'}'
 FILH2O='{output / 'water_reference.H2O'}'
 FILBAS='{adapted}'
 FILPS='{native / 'fit.PS'}'
 FILPRI='{native / 'fit.PRINT'}', LPRINT=6
 FILTAB='{native / 'fit.TABLE'}', LTABLE=7
 FILCOO='{native / 'fit.COORD'}', LCOORD=9
 DOWS=T, DOECC=F, WCONC={water_mM:.16g}, ATTH2O=1., ATTMET=1.
 WSMET='NAA', WSPPM=2.01, N1HMET=3, FWHMBA={fwhmba:.16g}
 NSIMUL=0, NRATIO=0, NOMIT={len(omitted)}
""" + ''.join(f" CHOMIT({index})='{name}'\n" for index, name in enumerate(omitted, 1)) + ' $END\n'
    (output / 'lcmodel.CONTROL').write_text(control)
    (output / 'input_hashes.json').write_text(json.dumps(before, indent=2) + '\n')
    started = time.monotonic()
    try:
        process = subprocess.run([str(BINARY)], input=control, text=True, cwd=output,
                                 capture_output=True, timeout=300)
        execution = dict(returncode=process.returncode, stdout=process.stdout, stderr=process.stderr)
    except subprocess.TimeoutExpired:
        execution = dict(returncode=None, stdout='', stderr='LCModel exceeded 300 seconds')
    execution['elapsed_s'] = time.monotonic() - started
    (output / 'execution.json').write_text(json.dumps(execution, indent=2) + '\n')
    table_path = native / 'fit.TABLE'
    table = table_path.read_text() if table_path.exists() else ''
    values = {}
    for line in table.splitlines():
        match = re.match(r'^\s*([\d.Ee+-]+)\s+(\d+(?:\.\d+)?)%\s+([\d.Ee+-]+)\s+(\S+)\s*$', line)
        if match:
            values[match[4]] = {'mM': float(match[1]), 'percent_SD': float(match[2])}
    valid = (execution['returncode'] == 0 and not re.search(r'\bFATAL\b', table, re.I)
             and all(name in values and np.isfinite(values[name]['mM']) for name in names))
    diagnostics = {'execution': execution, 'table_messages': [line.strip() for line in table.splitlines()
        if re.search(r'\b(fatal|warning|info|error)\b', line, re.I)],
        'missing_metabolites': [name for name in names if name not in values],
        'percent_SD': {name: values[name]['percent_SD'] for name in names if name in values},
        'water_concentration_mM': water_mM, 'FWHMBA_ppm': fwhmba,
        'basis_conversion': conversion, 'attenuation': 'ATTH2O=ATTMET=1; voxel-average mM'}
    if valid:
        coord = parse_coord(native / 'fit.COORD')
        norm = np.linalg.norm(coord['data'])
        diagnostics['relative_phased_real_residual'] = float(np.linalg.norm(coord['data'] - coord['fit']) / norm) if norm else None
    if before != hashes():
        raise RuntimeError('source inputs changed during LCModel fit')
    result = {'estimates_mM': {name: values[name]['mM'] for name in names if name in values},
              'status': 'ok' if valid else 'failed', 'diagnostics': diagnostics}
    (output / 'fitting_results.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result
