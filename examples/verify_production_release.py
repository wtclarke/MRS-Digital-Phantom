"""Verify the proposed production commit without staging or changing Git history.

Export the existing tracked baseline (omitting generated outputs/external anatomy),
overlay only the explicit production manifest, then test with the active Pixi Python.
Optional research files are absent. No fitting or regeneration is performed.
"""
import argparse
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile


def verify(output):
    repo=Path(__file__).resolve().parents[1];output=Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False);checkout=output/'production_checkout';checkout.mkdir()
    archive=subprocess.check_output(['git','archive','--format=tar','HEAD','--','.',
        ':(exclude)outputs',':(exclude)data/skeletons'],cwd=repo)
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:tar.extractall(checkout,filter='data')
    production=(repo/'docs/commit_manifests/production.txt').read_text().splitlines()
    research=set((repo/'docs/commit_manifests/research.txt').read_text().splitlines())
    assert set(production).isdisjoint(research)
    for name in production:
        target=checkout/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(repo/name,target)
    assert not (checkout/'outputs').exists()
    assert not any((checkout/name).exists() for name in research)
    environment=dict(os.environ,PYTHONPATH='',MPLCONFIGDIR=str(output/'matplotlib'))
    records=[]
    def execute(args):
        p=subprocess.run([sys.executable,*args],cwd=checkout,env=environment,capture_output=True,text=True,timeout=300)
        records.append({'command':[sys.executable,*args],'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr})
        (output/'verification.json').write_text(json.dumps(records,indent=2)+'\n')
        print(p.stdout,flush=True)
        if p.returncode:raise RuntimeError(p.stderr)
    execute(['-m','pytest','-q','tests/test_quantitative_production.py',
        'tests/test_quantitative_production_legacy.py','tests/test_calibration_resources.py'])
    execute(['-m','simulation.headless','--config','examples/configs/quantitative/pure_WM.json',
        '--output-dir','outputs/independent_cli'])
    # Build from this production-only source tree; no optional modules can leak in.
    execute(['-c',"import setuptools.build_meta; setuptools.build_meta.build_wheel('dist')"])
    import zipfile
    wheel=next((checkout/'dist').glob('*.whl'))
    with zipfile.ZipFile(wheel) as built:
        assert 'calibration/mrscloud_press_te30_quantitative/provenance_sources.zip' in built.namelist()
        assert 'simulation/quantitative.py' in built.namelist()
        assert 'simulation/mrscloud_water_calibration.py' not in built.namelist()
    print('Production-only checkout, CLI and wheel verified; Git index untouched.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    verify(parser.parse_args().output)
