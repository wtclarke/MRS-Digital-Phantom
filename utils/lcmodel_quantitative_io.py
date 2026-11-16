"""Explicit LCModel text adapters; no simulation or calibration changes."""
from pathlib import Path
import re

import numpy as np

NAMES=('NAA','Cr','GPC','PCh','Glu')


def write_raw(path,fid,*,frequency_mhz,echo_time_ms=30.):
    """LCModel stores conjugates of the validated positive-frequency FSL FID."""
    fid=np.asarray(fid)
    if fid.ndim!=1 or not np.iscomplexobj(fid) or not np.all(np.isfinite(fid)):
        raise ValueError('RAW requires a finite 1D complex FID')
    with Path(path).open('x') as stream:
        stream.write(f" $SEQPAR\n HZPPPM={frequency_mhz:.16g}, ECHOT={echo_time_ms:.16g}, SEQ='PRESS'\n $END\n")
        stream.write(" $NMID\n ID='Phase7a noise-free', FMTDAT='(2E25.16)', TRAMP=1., VOLUME=1.\n $END\n")
        for z in fid.conj():stream.write(f'{z.real:25.16E}{z.imag:25.16E}\n')


def parse_raw(path):
    """Independent ASCII parser returns STORED values, without hidden conjugation."""
    text=Path(path).read_text()
    header,body=text.rsplit('$END',1)
    values=np.fromstring(body.replace('D','E'),sep=' ')
    if not len(values) or len(values)%2 or not np.all(np.isfinite(values)):
        raise ValueError('invalid RAW complex pairs')
    data=values[::2]+1j*values[1::2]
    frequency=re.search(r'HZPPPM\s*=\s*([\d.eE+-]+)',header,re.I)
    return data,dict(header=header,frequency_mhz=float(frequency[1]) if frequency else None)


def parse_basis(path):
    text=Path(path).read_text()
    blocks=list(re.finditer(r'\$(\w+)\b(.*?)\$END',text,re.I|re.S))
    entries=[]
    for index,block in enumerate(blocks):
        if block[1].upper()!='BASIS':continue
        end=blocks[index+1].start() if index+1<len(blocks) else len(text)
        tokens=np.fromstring(text[block.end():end].replace('D','E'),sep=' ')
        if not len(tokens) or len(tokens)%2 or not np.all(np.isfinite(tokens)):raise ValueError('invalid basis spectrum')
        name=re.search(r"METABO\s*=\s*'([^']+)'",block[2],re.I)[1]
        entries.append((name,block[2],tokens[::2]+1j*tokens[1::2]))
    n=int(re.search(r'NDATAB\s*=\s*(\d+)',text,re.I)[1])
    dt=float(re.search(r'BADELT\s*=\s*([\d.eE+-]+)',text,re.I)[1])
    cf=float(re.search(r'HZPPPM\s*=\s*([\d.eE+-]+)',text,re.I)[1])
    if not entries or any(len(e[2])!=n for e in entries):raise ValueError('basis sample count mismatch')
    return entries,dict(points=n,dwell_time_s=dt,frequency_mhz=cf)


def compatible_basis(source,target,*,data_origin_ppm):
    """Preserve all spectra; supply FMTBAS and cancel LCModel's implicit 4.65 shift."""
    entries,metadata=parse_basis(source);text=Path(source).read_text()
    step=1/(metadata['points']*metadata['dwell_time_s']*metadata['frequency_mhz'])
    auto_shift=int(np.copysign(np.floor(abs((4.65-data_origin_ppm)/step)+.5),4.65-data_origin_ppm))
    prefix=text[:re.search(r'\$BASIS\b',text,re.I).start()]
    prefix=re.sub(r"\bFMTBAS\s*=\s*'[^']*'\s*,?",'',prefix,flags=re.I)
    prefix=re.sub(r'(\$BASIS1\b)',r"\1\n FMTBAS='(2E25.16)',",prefix,flags=re.I)
    with Path(target).open('x') as stream:
        stream.write(prefix)
        for name,header,values in entries:
            old=int(re.search(r'ISHIFT\s*=\s*([-+]?\d+)',header,re.I)[1])
            header=re.sub(r'(ISHIFT\s*=\s*)[-+]?\d+',lambda m:m[1]+str(old-auto_shift),header,flags=re.I)
            stream.write(' $BASIS'+header+'$END\n')
            for z in values:stream.write(f'{z.real:25.16E}{z.imag:25.16E}\n')
    after,_=parse_basis(target)
    error=max(float(abs(a[2]-b[2]).max()) for a,b in zip(entries,after))
    if error!=0:raise ValueError('basis serialization changed spectral samples')
    return dict(**metadata,spectrum_max_error=error,implicit_LCModel_origin_ppm=4.65,
        declared_source_and_data_origin_ppm=data_origin_ppm,basis_grid_ppm=step,
        added_ISHIFT=-auto_shift,net_LCModel_referencing_shift=0,
        reason='missing required FMTBAS and LCModel implicit 4.65-ppm basis reference; no amplitude change')


def control_text(directory,basis,*,names,points,dwell_time_s,frequency_mhz,origin_ppm,water_mol_l):
    directory=Path(directory).resolve();output=directory/'lcmodel_outputs'
    omit=[n for n in names if n not in NAMES]
    if not set(NAMES)<=set(names):raise ValueError('basis missing required metabolites')
    for path in [directory,basis]:
        if "'" in str(path):raise ValueError('single quotes are unsupported in LCModel paths')
    return f""" $LCMODL
 KEY=210387309
 TITLE='Phase7b exact Phase7a noise-free acquisition'
 NUNFIL={points}, DELTAT={dwell_time_s:.16g}, HZPPPM={frequency_mhz:.16g}, ECHOT=30.
 PPMCEN={origin_ppm:.16g}, PPMH2O=4.70, PPMST=4.2, PPMEND=.2
 FILRAW='{directory/'metabolites.RAW'}'
 FILH2O='{directory/'water_reference.H2O'}'
 FILBAS='{Path(basis).resolve()}'
 FILPS='{output/'fit.PS'}'
 FILPRI='{output/'fit.PRINT'}', LPRINT=6
 FILTAB='{output/'fit.TABLE'}', LTABLE=7
 FILCOO='{output/'fit.COORD'}', LCOORD=9
 DOWS=T, DOECC=F, WCONC={water_mol_l*1000:.16g}, ATTH2O=1., ATTMET=1.
 WSMET='NAA', WSPPM=2.01, N1HMET=3
 FWHMBA={1/frequency_mhz:.16g}
 NSIMUL=0, NRATIO=0, NOMIT={len(omit)}
"""+''.join(f" CHOMIT({i})='{name}'\n" for i,name in enumerate(omit,1))+" $END\n"


def parse_table(path):
    text=Path(path).read_text();concentrations={}
    for line in text.splitlines():
        match=re.match(r'^\s*([\d.Ee+-]+)\s+(\d+(?:\.\d+)?)%\s+([\d.Ee+-]+)\s+(\S+)\s*$',line)
        if match and match[4] in NAMES:
            concentrations[match[4]]=dict(concentration=float(match[1]),percent_SD=float(match[2]),relative_to_reference=float(match[3]))
    diagnostics=[]
    for line in text.splitlines():
        match=re.match(r'^\s*(\d+)\s+(fatal|warning|info|error)\s+(\w+)\s+(\d+)\s*$',line,re.I)
        if match:diagnostics.append(dict(count=int(match[1]),severity=match[2].lower(),routine=match[3],code=int(match[4])))
    version=re.search(r'Version\s+([^\)]+)',text)
    return dict(version=version[1] if version else None,concentrations=concentrations,diagnostics=diagnostics,
                fatal=bool(re.search(r'\bFATAL\b',text)))


def parse_coord(path):
    text=Path(path).read_text();n=int(re.search(r'(\d+) points on ppm-axis = NY',text)[1])
    markers={'ppm':'points on ppm-axis = NY','data':'NY phased data points follow',
             'fit':'NY points of the fit to the data follow','baseline':'NY background values follow'}
    result={}
    for key,marker in markers.items():
        body=text.split(marker,1)[1].splitlines()[1:]
        tokens=[]
        for line in body:
            tokens.extend(line.split())
            if len(tokens)>=n:break
        result[key]=np.array(tokens[:n],dtype=float)
        if len(result[key])!=n or not np.all(np.isfinite(result[key])):raise ValueError('invalid COORD array')
    return result


def local_peak_area(real,center,half_width,baseline_points,step):
    """Independent port of symmetric-minima/bordering-mean AREABA/INTEGRATE."""
    low=center-half_width;high=center+half_width
    left=low+int(np.argmin(real[low:center]))
    right=center+1+int(np.argmin(real[center+1:high+1]))
    distance=max(center-left,right-center);low=center-distance;high=center+distance
    background=(np.mean(real[low-baseline_points:low])+np.mean(real[high+1:high+1+baseline_points]))/2
    return float((np.sum(real[low:high+1])-(high-low+1)*background)*step)


def independent_water_scaling(basis_path,water_raw,*,origin_ppm=4.68,water_mM=39016.5):
    """Source-derived default integration for this pure, zero-phase water control.

    This does not reimplement LCModel's full GETPHA search or nonlinear fit.
    Known zero initial water phase is used, not a fitted calibration multiplier.
    """
    entries,meta=parse_basis(basis_path);spectrum=dict((e[0],e[2]) for e in entries)['NAA']
    step=1/(meta['points']*meta['dwell_time_s']*meta['frequency_mhz'])
    center=int(np.floor((origin_ppm-2.01)/step+.5))
    search=int(np.floor(20*.004/step+.5))
    center=center-search+int(np.argmax(spectrum.real[center-search:center+search+1]))
    half=int(np.floor(.5*10*(1/meta['frequency_mhz'])/step+.5))
    area_met=local_peak_area(spectrum.real,center,half,int(np.floor(.1/step+.5)),step)/3
    n=len(water_raw);dt=.0005;data_step=1/(2*n*dt*meta['frequency_mhz']);water_step=2*data_step
    water_spectrum=np.fft.fftshift(np.fft.fft(water_raw,norm='ortho'))/np.sqrt(2)
    ppm=origin_ppm-np.fft.fftshift(np.fft.fftfreq(n,dt))/meta['frequency_mhz']
    center=int(np.argmax(np.where(abs(ppm-4.70)<=1,abs(water_spectrum),0)))
    area_water=local_peak_area(water_spectrum.real,center,int(np.floor(2/water_step+.5)),
                              int(np.floor(.2/data_step+.5)),water_step)
    factor=area_met*2*water_mM/area_water
    nfreq=round(meta['points']*meta['dwell_time_s']/dt)
    # LCModel uses unitary transforms, while the stored source FFT is unnormalized.
    raw_basis_factor=meta['points']/np.sqrt(nfreq)
    return dict(normalized_basis_area=area_met,water_area=area_water,water_norm=area_water/(2*water_mM),
        FCALIB=factor,derived_basis_time_normalization=raw_basis_factor,
        ideal_coefficient_mM_to_report_factor=factor/raw_basis_factor/1000,
        scope='independent area/FFT-normalization prediction before nonlinear fit; zero-phase water control only')
