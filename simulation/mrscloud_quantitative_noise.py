"""Optional acquisition-level circular complex Gaussian noise; no calibration changes."""
from copy import deepcopy
from numbers import Integral

import numpy as np


def complex_gaussian_noise(points, *, sigma, seed):
    """sigma is the SD of EACH real/imaginary FID component, in signal units."""
    if isinstance(points,bool) or not isinstance(points,Integral) or points<=0:
        raise ValueError('points must be a positive integer')
    if not np.isfinite(sigma) or sigma<0:raise ValueError('sigma must be finite and nonnegative')
    if isinstance(seed,bool) or not isinstance(seed,(Integral,np.random.SeedSequence)) or (isinstance(seed,Integral) and seed<0):
        raise ValueError('seed must be a nonnegative integer or explicit SeedSequence')
    if sigma==0:return np.zeros(points,complex)
    values=np.random.Generator(np.random.PCG64(seed)).normal(0,sigma,size=(2,points))
    return values[0]+1j*values[1]


def reference_peak(fid, config, window_ppm):
    """Clean magnitude peak of an unweighted orthonormal FFT; no first-point half weighting."""
    fid=np.asarray(fid)
    if fid.ndim!=1 or not np.iscomplexobj(fid) or not np.all(np.isfinite(fid)):
        raise ValueError('reference must be a finite 1D complex FID')
    ppm=np.fft.fftshift(np.fft.fftfreq(len(fid),config['dwell_time_s']))/config['spectrometer_frequency_mhz']+config['frequency_origin_ppm']
    mask=(ppm>=window_ppm[0]) & (ppm<=window_ppm[1])
    if not mask.any():raise ValueError('reference window lies outside the spectral grid')
    return float(np.max(abs(np.fft.fftshift(np.fft.fft(fid,norm='ortho'))[mask])))


def add_pair_noise(pair, *, seed, metabolite_sigma=None, metabolite_snr=None,
                   water_sigma=None, water_snr=None):
    """Return a separate pair; never inject per cell or modify truth/components.

    SNR uses the clean isolated NAA magnitude peak at 1.8--2.2 ppm divided
    by the per-component spectral noise SD (equal to FID sigma for this FFT).
    Water uses its own clean magnitude peak at 4.1--5.3 ppm. Omitted noise
    specs disable that acquisition. Infinite target SNR means zero noise.
    """
    if isinstance(seed,bool) or not isinstance(seed,Integral) or seed<0:
        raise ValueError('noise seed must be an explicit nonnegative integer')
    if 'noise_parameters' in pair['simulation_config']:raise ValueError('noise has already been applied')
    result=deepcopy(pair);config=result['simulation_config']
    streams=np.random.SeedSequence(int(seed)).spawn(2)
    metadata={'root_seed':int(seed),'generator':'PCG64','numpy_version':np.__version__,
        'convention':'clean isolated reference magnitude peak / per-component spectral SD',
        'FFT':'fftshift(fft(fid,norm=ortho)); no first-point weighting',
        'domain':'frequency','sigma_domain':'FID real or imaginary component; orthonormal spectral SD is identical',
        'complex_noise_power':'E[abs(eta)^2]=2*sigma^2',
        'injection_stage':'after deterministic summation, relaxation and spatial B0; once per acquisition',
        'components_are_clean':True}
    for index,(label,key,sigma,snr,window) in enumerate([
        ('metabolite','metabolite_fid',metabolite_sigma,metabolite_snr,(1.8,2.2)),
        ('water','water_reference_fid',water_sigma,water_snr,(4.1,5.3))]):
        if sigma is not None and snr is not None:raise ValueError('specify sigma or SNR, not both')
        clean=result[key].copy()
        if clean.ndim!=1 or not np.iscomplexobj(clean) or not np.all(np.isfinite(clean)):
            raise ValueError('noise requires finite one-dimensional complex FIDs')
        reference=result['components'].get('NAA') if label=='metabolite' else clean
        peak=reference_peak(reference,config,window) if reference is not None else None
        if snr is not None:
            if np.isnan(snr) or snr<=0:raise ValueError('SNR must be positive or infinite')
            if peak is None or peak<=0:raise ValueError('target SNR requires a nonzero clean reference (NAA for metabolites)')
            sigma=peak/float(snr)
        sigma=0. if sigma is None else float(sigma)
        noise=complex_gaussian_noise(len(clean),sigma=sigma,seed=streams[index])
        result['clean_'+key]=clean;result[key]=clean+noise
        spectrum_key='metabolite_spectrum' if label=='metabolite' else 'water_reference_spectrum'
        if spectrum_key in result:result[spectrum_key]=np.fft.fftshift(np.fft.fft(result[key]))
        # This is a realized statistical estimate, not a renormalization of the noise draw.
        spec_noise=np.fft.fftshift(np.fft.fft(noise,norm='ortho'))
        estimated_sd=float(np.sqrt((np.var(spec_noise.real,ddof=1)+np.var(spec_noise.imag,ddof=1))/2))
        metadata[label]={'sigma_per_component':sigma,'requested_snr':('infinite' if np.isinf(snr) else float(snr)) if snr is not None else None,
            'reference_peak':peak,'reference_window_ppm':list(window),'stream_spawn_key':list(streams[index].spawn_key),
            'estimated_spectral_component_sd':estimated_sd,
            'realized_snr':peak/estimated_sd if peak is not None and estimated_sd>0 else None}
    config.update(noise=bool(metadata['metabolite']['sigma_per_component'] or metadata['water']['sigma_per_component']),
                  noise_parameters=metadata)
    return result
