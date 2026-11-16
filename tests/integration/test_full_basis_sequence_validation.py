import json
from pathlib import Path
import subprocess
import shutil

import numpy as np

ROOT=Path('outputs/mrscloud_full_basis/mrscloud_press_te30_quantitative')


def test_smoke_bundle_and_raw_cache_reuse_without_simulation():
    from examples.mrscloud_full_basis_water_reference import prepare_capture
    from simulation.mrscloud_full_basis import load_library, sha256
    source=Path('outputs/mrscloud_full_basis_smoke_final')
    cfg,_,entries,_,command=prepare_capture(source,smoke=True)
    files=[source/'raw'/label/'response.mat' for label,_,_ in entries]
    hashes={str(p):sha256(p) for p in files}
    completed=subprocess.run(command,capture_output=True,text=True,check=True)
    assert completed.stdout.count('CACHE ')==36 and 'COMPLETE ' not in completed.stdout
    assert hashes=={str(p):sha256(p) for p in files}
    r,w,m=load_library(source/'mrscloud_press_te30_quantitative',points=128)
    assert len(r)==27 and w.shape==(128,)
    assert m['sequence']['spatial_points_per_axis']==3


def test_actual_full_resolution_determinism_proton_controls_and_cross_calibration():
    result=json.loads((ROOT/'validation_results.json').read_text())
    for case in result['determinism'].values():
        assert all(error==0 for error in case.values())
    assert abs(result['proton_controls']['on_resonance_methyl_per_proton_vs_water']-1)<1e-12
    assert abs(result['proton_controls']['RF3_methyl_per_proton_vs_water']-1)<3e-6
    cross=result['cross_calibration']
    assert cross['NAA_relative_error']<3e-7
    assert cross['water_relative_error']<1e-14
    assert max(result['selection_max_errors'].values())==0
    assert result['mixture_superposition_max_error']==0
    assert result['mixture_export_max_error']==0


def test_checked_basis_comparisons_use_single_alignment_and_report_every_column():
    result=json.loads((ROOT/'validation_results.json').read_text())
    rows=result['checked_basis_comparison']
    assert len(rows)==29
    for row in rows:
        assert np.isfinite(row['relative_complex_error'])
        assert np.isfinite(row['spectral_integral_ratio'])
        assert row['distance_to_checked_near_max_ppm']<=2000/2048/127.731
        assert abs(complex(row['correlation_real'],row['correlation_imag']))>.99


def test_cached_raw_rebuild_reproduces_normalization_and_arrays_without_octave(tmp_path,monkeypatch):
    from examples.mrscloud_full_basis_water_reference import build
    from simulation.mrscloud_full_basis import load_library
    source=ROOT.parent/'raw'
    for path in source.glob('*/response.mat'):
        destination=tmp_path/'raw'/path.parent.name/'response.mat'
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,destination)
    def forbidden(*args,**kwargs):
        raise AssertionError('cached assembly must not invoke Octave')
    monkeypatch.setattr('examples.mrscloud_full_basis_water_reference.subprocess.run',forbidden)
    rebuilt,_=build(tmp_path,reuse_raw=True)
    a,wa,ma=load_library(ROOT)
    b,wb,mb=load_library(rebuilt)
    assert ma==mb
    for name in a:np.testing.assert_array_equal(a[name],b[name])
    np.testing.assert_array_equal(wa,wb)
