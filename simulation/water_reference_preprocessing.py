"""Optional deterministic alignment of a separately acquired water reference."""

import numpy as np


def align_water_reference(fid, *, dwell_time_s, frequency_offset_hz, initial_phase_rad):
    """Return a new water FID demodulated to DC and zero initial phase.

    Corrections refer to interpreted FID samples, before NIfTI-MRS storage
    conjugation. They must be explicitly supplied; no noisy-data estimator is
    implied. Input samples and their magnitudes are preserved to roundoff.
    """
    data = np.asarray(fid)
    if data.ndim != 1 or not data.size or not np.iscomplexobj(data) or not np.all(np.isfinite(data)):
        raise ValueError('water must be a finite, nonempty, one-dimensional complex FID')
    if not np.isfinite(dwell_time_s) or dwell_time_s <= 0:
        raise ValueError('dwell time must be finite and positive')
    if not np.isfinite(frequency_offset_hz) or abs(frequency_offset_hz) >= .5 / dwell_time_s:
        raise ValueError('frequency correction must be finite and inside the sampled bandwidth')
    if not np.isfinite(initial_phase_rad):
        raise ValueError('phase correction must be finite')
    t = np.arange(data.size) * dwell_time_s
    aligned = data * np.exp(-1j * (2 * np.pi * frequency_offset_hz * t + initial_phase_rad))
    metadata = {
        'method': 'known-frequency-demodulation-and-zero-order-phase-v1',
        'applied_frequency_correction_hz': -float(frequency_offset_hz),
        'applied_phase_correction_rad': -float(initial_phase_rad),
        'input_frequency_offset_hz': float(frequency_offset_hz),
        'input_initial_phase_rad': float(initial_phase_rad),
        'dwell_time_s': float(dwell_time_s), 'points': int(data.size),
        'time_origin_s': 0.0, 'magnitude_scaling': 1.0,
        'aligned_frequency_offset_hz': 0.0, 'aligned_initial_phase_rad': 0.0,
        'interpretation': 'digitally aligned analysis reference; physical water chemical shift unchanged',
    }
    return aligned, metadata
