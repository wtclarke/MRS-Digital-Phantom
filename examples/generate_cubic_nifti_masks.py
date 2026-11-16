"""Generate numbered, axis-aligned cubic masks on a reference NIfTI grid."""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


DEFAULT_REFERENCE = Path("data/skeletons/MRiLab/T1_map.nii.gz")


def generate_masks(reference, output_dir, size_mm=(20.0, 20.0, 20.0),
                   start=(50, 50, 50), stop=(120, 120, 120), stride=10):
    reference = nib.load(str(reference))
    if len(reference.shape) != 3:
        raise ValueError("Reference image must be 3D")
    size_mm = np.asarray(size_mm, dtype=float)
    spacing = np.linalg.norm(reference.affine[:3, :3], axis=0)
    if np.any(size_mm <= 0) or np.any(spacing <= 0):
        raise ValueError("Cube dimensions and voxel spacing must be positive")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    starts, stops = np.asarray(start), np.asarray(stop)
    if stride <= 0 or np.any(stops < starts):
        raise ValueError("Stride must be positive and stop must be >= start")
    centers = [tuple((starts + i * stride).tolist())
               for i in range(int(np.min((stops - starts) // stride)) + 1)]
    if centers[-1] != tuple(stops.tolist()):
        raise ValueError("Stop center must be reachable from start using stride")

    shape = np.asarray(reference.shape)
    grid = np.ogrid[tuple(slice(0, n) for n in shape)]
    for label, center in enumerate(centers, 1):
        center = np.asarray(center)
        # Convert the requested physical extent into half-open voxel-center bounds.
        low = center - size_mm / (2 * spacing)
        high = center + size_mm / (2 * spacing)
        if np.any(low < -0.5) or np.any(high > shape - 0.5):
            raise ValueError(f"Cube at center {tuple(center)} extends outside reference image")
        mask = np.ones(shape, dtype=bool)
        for axis in range(3):
            mask &= (grid[axis] >= low[axis]) & (grid[axis] < high[axis])
        if not mask.any():
            raise ValueError(f"Cube at center {tuple(center)} selects no voxels")
        data = np.zeros(shape, dtype=np.uint16)
        data[mask] = label
        path = output_dir / f"mask_{label:03d}_{center[0]}_{center[1]}_{center[2]}.nii.gz"
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite {path}")
        nib.save(nib.Nifti1Image(data, reference.affine, header=reference.header.copy()), str(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--size-mm", type=float, nargs=3, default=(20, 20, 20), metavar=("X", "Y", "Z"))
    parser.add_argument("--start", type=int, nargs=3, default=(50, 50, 50), metavar=("I", "J", "K"))
    parser.add_argument("--stop", type=int, nargs=3, default=(120, 120, 120), metavar=("I", "J", "K"))
    parser.add_argument("--stride", type=int, default=10)
    args = parser.parse_args()
    generate_masks(args.reference, args.output_dir, args.size_mm, args.start, args.stop, args.stride)


if __name__ == "__main__":
    main()
