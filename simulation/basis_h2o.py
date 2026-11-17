"""Explicitly approximate water referencing from a supplied, unchanged LCModel BASIS."""
import hashlib
import json
from numbers import Integral
from pathlib import Path
import re

import fsl_mrs
import numpy as np
from fsl_mrs.utils import mrs_io

ASSUMPTION='assume_common_scale_protons_encoded'


def _field(header,name, *, string=False):
    pattern=rf"\b{name}\s*=\s*'([^']*)'" if string else rf'\b{name}\s*=\s*([-+\d.eEdD]+)'
    match=re.search(pattern,header,re.I)
    if match is None:raise ValueError(f'BASIS requires explicit {name}')
    return match[1].strip() if string else float(match[1].replace('D','E').replace('d','e'))


def load_responses(path, *, frequency_origin_ppm, bandwidth_hz=2000,points=2048,
                   water_component_name='H2O',expected_basis_sha256=None):
    """Format the entire fixed reference without any additional amplitude scaling."""
    path=Path(path);content=path.read_bytes();digest=hashlib.sha256(content).hexdigest()
    if expected_basis_sha256 is not None and digest!=expected_basis_sha256:
        raise ValueError('BASIS SHA256 mismatch')
    if (isinstance(points,bool) or not isinstance(points,Integral) or points<2
            or isinstance(bandwidth_hz,bool) or not np.isfinite(bandwidth_hz) or bandwidth_hz<=0
            or not np.isfinite(frequency_origin_ppm) or not 0<=frequency_origin_ppm<=10):
        raise ValueError('finite Hz bandwidth, integer points >=2 and declared ppm origin are required')
    text=content.decode('ascii')
    blocks=re.findall(r'\$(SEQPAR|BASIS1|BASIS)\b(.*?)\$END',text,re.I|re.S)
    headers={kind.upper():header for kind,header in blocks if kind.upper()!='BASIS'}
    entries=[header for kind,header in blocks if kind.upper()=='BASIS']
    if any(sum(kind.upper()==key for kind,_ in blocks)!=1 for key in ['SEQPAR','BASIS1']) or not entries:
        raise ValueError('missing LCModel BASIS headers/components')
    native_points=_field(headers['BASIS1'],'NDATAB');dt=_field(headers['BASIS1'],'BADELT')
    cf=_field(headers['SEQPAR'],'HZPPPM');te=_field(headers['SEQPAR'],'ECHOT')
    if not all(np.isfinite(v) and v>0 for v in [native_points,dt,cf,te]) or native_points!=int(native_points):
        raise ValueError('inconsistent BASIS sampling/sequence metadata')
    component_headers={}
    for header in entries:
        name=_field(header,'METABO',string=True)
        if name in component_headers:raise ValueError('duplicate BASIS component names')
        values={key:_field(header,key) for key in ['CONC','TRAMP','VOLUME','ISHIFT']}
        if any(values[key]!=1. for key in ['CONC','TRAMP','VOLUME']):
            raise ValueError('unsupported CONC/TRAMP/VOLUME: common scale cannot be established')
        if not np.isfinite(values['ISHIFT']) or values['ISHIFT']!=int(values['ISHIFT']):
            raise ValueError('invalid component frequency shift')
        component_headers[name]=values
    basis=mrs_io.read_basis(str(path))
    names=list(map(str,basis.names));native=basis.original_basis_array
    if (len(set(names))!=len(names) or set(names)!=set(component_headers)
            or basis.original_points!=native_points or native.shape!=(int(native_points),len(names))
            or not np.isclose(basis.original_dwell,dt,rtol=1e-12,atol=0)
            or not np.isclose(basis.cf,cf,rtol=1e-12,atol=0)
            or not np.iscomplexobj(native) or not np.all(np.isfinite(native))):
        raise ValueError('inconsistent or nonfinite BASIS complex samples/sampling')
    if water_component_name not in names:raise ValueError(f'BASIS has no {water_component_name} component')
    if bandwidth_hz>1/dt or (points-1)/bandwidth_hz>(native_points-1)*dt:
        raise ValueError('requested sampling exceeds native BASIS bandwidth/time coverage')
    # Full reference always includes water; no legacy mean-magnitude or per-column scaling.
    formatted=basis.get_formatted_basis(bandwidth_hz,points,scale_factor=None)
    water=formatted[:,names.index(water_component_name)].copy()
    if not np.all(np.isfinite(formatted)) or not np.any(abs(water)>0):
        raise ValueError('H2O response must be finite and nonzero')
    frequencies=np.fft.fftshift(np.fft.fftfreq(points,1/bandwidth_hz))
    peak=float(frequencies[np.argmax(abs(np.fft.fftshift(np.fft.fft(water))))])
    descriptor={'basis_sha256':digest,'bandwidth_hz':float(bandwidth_hz),'points':int(points),
        'frequency_origin_ppm':float(frequency_origin_ppm),'water_component_name':water_component_name,
        'reference_names':names,'normalization':'identity_preserve_supplied_common_scale',
        'amplitude_assumption':ASSUMPTION,'fsl_mrs_version':fsl_mrs.__version__,'numpy_version':np.__version__}
    identity='basis-h2o-approximate-v1:'+hashlib.sha256(json.dumps(descriptor,sort_keys=True).encode()).hexdigest()
    metadata={**descriptor,'calibration_identity':identity,'basis_filename':path.name,
        'physical_calibration':False,'approximate_calibration':True,
        'normalization':{'scalar':1.,'operation':'divide all formatted responses by 1; no extra normalization',
            'statistic':'identity_preserve_supplied_common_scale','reference_names':names,
            'selection_independent':True},'component_headers':component_headers,
        'sequence':{'sequence':_field(headers['SEQPAR'],'SEQ',string=True),'TE_ms':te,
            'spectrometer_frequency_mhz':cf,'frequency_origin_ppm':float(frequency_origin_ppm),
            'native_points':int(native_points),'native_bandwidth_hz':1/dt,
            'water_frequency_offset_hz':peak,'water_phase_rad':float(np.angle(water[0])),
            'water_frequency_estimate':'magnitude FFT maximum bin, not a fitted carrier',
            'water_proton_multiplier_applied':False},
        'limitations':['Common molar response and encoded proton counts are explicitly assumed, not proven by CONC=1.',
            'RF carrier/refocusing sensitivity, intrinsic phase, DC correction, first point and readout are retained.',
            'No off-resonance compensation or per-metabolite correction; not a water-centred acquisition.']}
    return {n:formatted[:,i].copy() for i,n in enumerate(names) if n not in {water_component_name,'Ref0ppm'}},water,metadata
