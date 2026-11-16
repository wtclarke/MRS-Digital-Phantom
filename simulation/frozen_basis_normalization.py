"""Selection-independent synthetic normalization for formatted MRSCloud bases."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from fsl_mrs.utils.mrs_io.main import read_basis


@dataclass(frozen=True)
class FrozenBasisCalibration:
    reference_identity: str
    bandwidth_hz: float
    points: int
    normalization_scalar: float
    statistic: str = "mean_absolute_sample_over_all_reference_fids"

    def validate(self, *, reference_identity, bandwidth_hz, points):
        if (reference_identity, float(bandwidth_hz), int(points)) != (
            self.reference_identity, self.bandwidth_hz, self.points
        ):
            raise ValueError("basis calibration metadata mismatch")

    def to_dict(self):
        return {"reference_identity": self.reference_identity,
                "bandwidth_hz": self.bandwidth_hz, "points": self.points,
                "normalization_scalar": self.normalization_scalar,
                "statistic": self.statistic}


def basis_identity(path):
    path = Path(path)
    return f"{path.name}:{hashlib.sha256(path.read_bytes()).hexdigest()}"


def format_frozen_basis(path, *, bandwidth_hz, points, selected_names=None, calibration=None):
    """Format a fixed reference basis, normalize once, then select a subset."""
    identity = basis_identity(path)
    basis = read_basis(str(path))
    names = list(basis.names)
    formatted_all = basis.get_formatted_basis(bandwidth_hz, points)
    if calibration is None:
        scalar = float(np.mean(np.abs(formatted_all)))
        calibration = FrozenBasisCalibration(identity, float(bandwidth_hz), int(points), scalar)
    else:
        calibration.validate(reference_identity=identity, bandwidth_hz=bandwidth_hz, points=points)
        scalar = calibration.normalization_scalar
    if not np.isfinite(scalar) or scalar <= 0:
        raise ValueError("normalization scalar must be finite and positive")
    wanted = names if selected_names is None else list(selected_names)
    unknown = set(wanted) - set(names)
    if unknown:
        raise ValueError(f"metabolites absent from reference basis: {sorted(unknown)}")
    idx = [names.index(name) for name in wanted]
    return formatted_all[:, idx] / scalar, wanted, calibration
