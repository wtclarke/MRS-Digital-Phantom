"""Explicit pre-acquisition attenuation; all times are in seconds."""
import numpy as np

TISSUES = ('WM', 'GM', 'CSF')


def attenuation(*, TE_s, TR_s, T1_s, T2_s):
    """90-degree excitation, independent saturation and echo attenuation."""
    values = np.asarray([TE_s, TR_s, T1_s, T2_s], dtype=float)
    TE_s, TR_s, T1_s, T2_s = values
    if not np.all(np.isfinite(values)) or TE_s < 0 or TR_s <= TE_s or T1_s <= 0 or T2_s <= 0:
        raise ValueError('times must be finite seconds: TE >= 0, TR > TE, T1 > 0, T2 > 0')
    return float(-np.expm1(-TR_s/T1_s)*np.exp(-TE_s/T2_s))


def relaxation_metadata(config, names):
    """Validate complete explicit profiles, including zero-volume tissues."""
    if not isinstance(config, dict) or config.get('mode') not in ('none', 'tissue_specific'):
        raise ValueError('relaxation mode must be none or tissue_specific')
    required = {'mode', 'TE_s', 'TR_s'}
    active = config['mode'] == 'tissue_specific'
    if active:
        required |= {'water', 'metabolites'}
    if set(config) != required:
        raise ValueError(f'relaxation requires exactly {sorted(required)}; time units are seconds')
    te, tr = float(config['TE_s']), float(config['TR_s'])
    attenuation(TE_s=te, TR_s=tr, T1_s=1., T2_s=1.)
    meta = {'mode': config['mode'], 'TE_s': te, 'TR_s': tr,
            'equation': '(1-exp(-TR/T1))*exp(-TE/T2)', 'times_unit': 's',
            'application': 'pre_acquisition_amplitude_only',
            'basis_assumption': 'coherent PRESS TE evolution and readout decay retained; no T1/TR or TE/T2 attenuation encoded',
            'water_factors': dict.fromkeys(TISSUES, 1.),
            'metabolite_factors': {t: dict.fromkeys(names, 1.) for t in TISSUES}}
    if not active:
        return meta
    for species in ('water', 'metabolites'):
        if not isinstance(config[species], dict) or set(config[species]) != set(TISSUES):
            raise ValueError('relaxation requires explicit WM, GM and CSF profiles')
    meta.update(water={}, metabolites={})
    def times(profile):
        if not isinstance(profile, dict) or set(profile) != {'T1_s', 'T2_s'}:
            raise ValueError('each relaxation profile requires T1_s and T2_s in seconds')
        result = {k: float(v) for k, v in profile.items()}
        factor = attenuation(TE_s=te, TR_s=tr, **result)
        return result, factor
    for tissue in TISSUES:
        meta['water'][tissue], meta['water_factors'][tissue] = times(config['water'][tissue])
        profiles = config['metabolites'][tissue]
        if not isinstance(profiles, dict) or set(profiles) != set(names):
            raise ValueError('relaxation profiles must exactly match the selected metabolite names')
        meta['metabolites'][tissue] = {}
        for name in names:
            meta['metabolites'][tissue][name], meta['metabolite_factors'][tissue][name] = times(profiles[name])
    return meta
