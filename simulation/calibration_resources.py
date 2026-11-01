"""Path-independent lookup and exact frozen generator-provenance verification."""
import hashlib
from importlib.resources import files
import json
from pathlib import Path
import zipfile

BUNDLE_NAME='mrscloud_press_te30_quantitative'


def default_calibration_bundle():
    return Path(str(files('calibration').joinpath(BUNDLE_NAME)))


def calibration_bundle(directory=None):
    root=default_calibration_bundle() if directory is None else Path(directory).expanduser().resolve()
    if not (root/'reference_manifest.json').is_file():
        raise FileNotFoundError(f'Frozen PRESS TE30 calibration missing at {root}. Restore calibration/{BUNDLE_NAME} from this repository or set quantitative.calibration_bundle explicitly. No download or regeneration is performed.')
    return root


def verify_source_provenance(directory,normalization):
    """Every historical source byte must match; sources are never imported/run."""
    root=Path(directory);archive=root/'provenance_sources.zip'
    if not archive.exists():
        # An explicit legacy bundle can share the same immutable source snapshot.
        stable=default_calibration_bundle()
        snapshot=json.loads((stable/'normalization.json').read_text())
        for key in ['calibration_identity','source_hashes','processing_source_hashes']:
            if snapshot[key]!=normalization[key]:
                raise ValueError('no matching frozen source provenance; provide provenance_sources.zip with this calibration')
        archive=stable/'provenance_sources.zip'
    expected={**normalization['source_hashes'],**normalization['processing_source_hashes']}
    try:
        with zipfile.ZipFile(archive) as source:
            if len(source.namelist())!=len(expected) or set(source.namelist())!=set(expected):
                raise ValueError('frozen source provenance inventory mismatch')
            for path,digest in expected.items():
                if hashlib.sha256(source.read(path)).hexdigest()!=digest:
                    raise ValueError(f'frozen source provenance identity mismatch: {path}')
    except (FileNotFoundError,zipfile.BadZipFile,KeyError) as error:
        raise ValueError(f'frozen source provenance unavailable or corrupt: {archive}') from error
