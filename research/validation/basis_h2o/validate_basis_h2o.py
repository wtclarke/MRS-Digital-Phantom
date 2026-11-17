"""Compare a supplied BASIS with its companion MAT; no fitting or amplitude scaling."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils import mrs_io


def compare_source_mat(basis_path,source_mat):
    """One declared acquisition-frame shift, no amplitude alignment/fitting."""
    from scipy.io import loadmat
    mat=loadmat(source_mat,simplify_cells=True)['BASIS']
    basis=mrs_io.read_basis(str(basis_path))
    columns=np.column_stack([mat['fids'][:,list(mat['name']).index(n)] for n in basis.names])
    time=np.arange(basis.original_points)*basis.original_dwell
    shift=(float(mat['centerFreq'])-4.68)*basis.cf
    transformed=columns*np.exp(2j*np.pi*shift*time)[:,None]
    errors=np.linalg.norm(transformed-basis.original_basis_array,axis=0)/np.linalg.norm(transformed,axis=0)
    return {'source_mat_sha256':hashlib.sha256(source_mat.read_bytes()).hexdigest(),
        'source_global_scale':float(mat['scale']),'rf_carrier_ppm':float(mat['centerFreq']),
        'acquisition_origin_ppm':4.68,'frame_shift_hz':shift,'amplitude_factor_applied':1.,
        'relative_errors':dict(zip(map(str,basis.names),map(float,errors))),
        'max_relative_error':float(errors.max()),
        'normalized_source_max_real_spectrum':float(np.max(mat['specs'].real)),
        'individual_max_real_spectrum_range':[float(np.max(mat['specs'].real,axis=0).min()),float(np.max(mat['specs'].real,axis=0).max())]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--basis', type=Path, required=True)
    parser.add_argument('--source-mat', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New JSON comparison file')
    args = parser.parse_args()
    result = compare_source_mat(args.basis, args.source_mat)
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
