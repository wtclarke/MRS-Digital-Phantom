"""Independent checks of mixture, perturbation and experiment coverage controls."""
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from research.validation.headless_water.run import (
    GROUPS, MIXTURES, SCALES, TARGETS, concentration_table, design, subject,
)


def test_group_perturbation_preserves_nuisance_metabolites_and_pcr_reference():
    names = ['NAA', 'PCr', 'Cr', 'GPC', 'PCh', 'bHG']
    base = concentration_table(names, 'baseline', 1).set_index(['Tissue', 'Metabolite'])
    changed = concentration_table(names, 'tCho', 1.5).set_index(['Tissue', 'Metabolite'])
    for tissue in ['WM', 'GM']:
        for name in names:
            multiplier = 1.5 if name in GROUPS['tCho'] else 1
            assert changed.loc[(tissue, name), 'Conc_mean'] == pytest.approx(multiplier * base.loc[(tissue, name), 'Conc_mean'])
    assert (changed.loc[['Background', 'CSF', 'Fat', 'Skull'], 'Conc_mean'] == 0).all()
    zero = concentration_table(names, 'NAA', 0)
    assert (zero.loc[zero.Metabolite == 'NAA', 'Conc_mean'] == 0).all()
    assert (zero.loc[(zero.Metabolite == 'PCr') & zero.Tissue.isin(['WM', 'GM']), 'Conc_mean'] > 0).all()


@pytest.mark.parametrize('fractions', MIXTURES)
def test_all_mixtures_and_fixed_spatial_field_shape(fractions):
    case = dict(zip(['wm', 'gm', 'csf'], fractions), b0_sd=5.)
    a = subject(case)
    roi = a.labels.data[0, 1:11, 1:11, 1:11].numpy()
    for label, expected in enumerate(fractions, 1):
        assert (roi == label).mean() == pytest.approx(expected)
    field = a.B0.data[0, 1:11, 1:11, 1:11].numpy()
    assert field.std() == pytest.approx(5., rel=1e-6)
    assert abs(field.mean()) < 1e-6
    b = subject(dict(case, b0_sd=10.))
    np.testing.assert_array_equal(a.labels.data, b.labels.data)
    np.testing.assert_array_equal(2 * a.B0.data, b.B0.data)
    assert a.spacing == (2., 2., 2.)


def test_matrix_coverage_without_external_anatomy(tmp_path, monkeypatch):
    for vendor in ['GE', 'PH', 'SH', 'UI']:
        for te in [35, 97, 144]:
            (tmp_path / f'{vendor}_PRESS_TE{te}.BASIS').touch()
    labels = torch.ones(1, 64, 64, 64, dtype=torch.int16)
    labels[:, 32:] = 2
    labels[:, :, 48:] = 3
    anatomy = SimpleNamespace(spacing=(2., 2., 2.), labels=SimpleNamespace(data=labels))
    monkeypatch.setattr('utils.loading.load_skeleton_as_subject', lambda *a: (anatomy, []))
    cases = design(tmp_path)
    for basis_name in ['GE_PRESS_TE35', 'SH_PRESS_TE97']:
        clean = [c for c in cases if c['basis_name'] == basis_name and c['stage'] == 'clean']
        assert len(clean) == 7 + 6 * len(TARGETS) * (len(SCALES) - 1)
        assert len([c for c in clean if c['csf'] == 1]) == 1
        assert all(c['b0_sd'] == 0 and np.isinf(c['snr']) for c in clean)
        for target in TARGETS:
            assert {c['scale'] for c in clean if c['target'] == target} == set(SCALES) - {1}
    transfer = [c for c in cases if c['stage'] == 'transfer']
    assert len({c['basis_name'] for c in transfer}) == 10
    assert len(transfer) == 50
    assert all(c['wm'] + c['gm'] + c['csf'] == pytest.approx(1.) for c in cases)
