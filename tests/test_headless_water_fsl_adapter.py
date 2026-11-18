"""Native water scaling must be converted back from tissue-only to voxel mM."""
import hashlib
import json

import numpy as np
import pytest
from nifti_mrs.create_nmrs import gen_nifti_mrs

from research.validation.headless_water.fsl_adapter import fit_case


@pytest.mark.parametrize('concentration', [0., 6.])
def test_native_quantification_matches_water_profile_and_voxel_target(tmp_path, basis, concentration):
    path, responses = basis
    path.write_text(path.read_text().replace("METABO='Cr'", "METABO='PCr'"))
    source = tmp_path / 'source'
    source.mkdir()
    water_molar = .5 * 39.0165 + .25 * 47.2305 + .25 * 55.5
    for name, fid in [('total', concentration * responses['NAA'] + 4 * responses['Cr']),
                      ('water_reference', 1000 * water_molar * responses['H2O'])]:
        obj = gen_nifti_mrs(fid[None, None, None, :], .0005, 100.)
        obj.add_hdr_field('SpecFreqChemShift', 4.68)
        obj.save(str(source / f'{name}.nii.gz'))
    truth = {'tissue_fractions': {'WM': .5, 'GM': .25, 'CSF': .25}}
    (source / 'ground_truth.json').write_text(json.dumps(truth))
    provenance = {'basis_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                  'water_options': {'water_concentrations_mol_l': {'WM': 39.0165, 'GM': 47.2305, 'CSF': 55.5}}}
    (source / 'calibration_provenance.json').write_text(json.dumps(provenance))
    output = tmp_path / 'fit'
    record = fit_case(source, path, output)
    assert record['status'] in ('ok', 'optimizer_warning'), record
    assert record['diagnostics']['water_molarity_mM'] == pytest.approx(1000 * water_molar)
    for name, value in record['native_estimates_mM'].items():
        assert record['estimates_mM'][name] == pytest.approx(.75 * value)
    assert np.isfinite(record['estimates_mM']['NAA'])
    if concentration > 0:
        assert record['estimates_mM']['NAA'] > 0
    assert record['diagnostics']['water_ref_metab'] == 'PCr'
    repeated = fit_case(source, path, tmp_path / 'repeated')
    assert repeated['estimates_mM'] == record['estimates_mM']
    assert (output / 'fit.pdf').is_file()
    assert (output / 'fit_arrays.npz').is_file()
    with pytest.raises(FileExistsError):
        fit_case(source, path, output)


def test_failure_is_recorded_without_fabricated_estimates(tmp_path, basis):
    source = tmp_path / 'source'
    source.mkdir()
    (source / 'ground_truth.json').write_text('{}')
    (source / 'calibration_provenance.json').write_text('{"basis_sha256": "wrong"}')
    record = fit_case(source, basis[0], tmp_path / 'failed')
    assert record['status'] == 'failed'
    assert record['estimates_mM'] == {}
    assert 'BASIS does not match' in record['error']
    assert 'Traceback' in (tmp_path / 'failed/execution.log').read_text()
