"""Opt-in tissue-weighted analytical synthetic MRS pair generation/export."""

import json
from pathlib import Path

import numpy as np
from nifti_mrs.create_nmrs import gen_nifti_mrs

from simulation.analytical_water_reference import concentration_mol_per_litre

CALIBRATION_ID = "synthetic-protons-per-mol-v1"
TISSUE_WATER_MOL_L = {"WM": 39.0165, "GM": 47.2305, "CSF": 55.5}


def simulate_analytical_voxel_pair(
    tissue_fractions,
    *,
    metabolite_concentrations_mM,
    water_concentrations_mol_l=None,
    n_points=256,
    dwell_time_s=5e-4,
    receiver_factor=1.0,
):
    """Generate same-grid metabolite and unsuppressed water FIDs by tissue.

    Concentrations are per total tissue volume. Metabolite response is 3
    equivalent protons; water is 2. No noise, relaxation, phase or offset.
    """
    fractions = {str(k): float(v) for k, v in tissue_fractions.items()}
    if not np.isclose(sum(fractions.values()), 1.0):
        raise ValueError("tissue fractions must sum to one")
    water_concentrations_mol_l = water_concentrations_mol_l or TISSUE_WATER_MOL_L
    missing = set(fractions) - set(metabolite_concentrations_mM)
    if missing:
        raise ValueError(f"Missing metabolite concentrations for tissues: {sorted(missing)}")
    t = np.arange(n_points) * dwell_time_s
    met = np.zeros(n_points, dtype=np.complex128)
    water = np.zeros_like(met)
    gt_met = {}
    gt_water = {}
    for tissue, fraction in fractions.items():
        cm = concentration_mol_per_litre(metabolite_concentrations_mM[tissue], "mM")
        cw = float(water_concentrations_mol_l[tissue])
        gt_met[tissue], gt_water[tissue] = cm, cw
        met += fraction * cm * 3 * receiver_factor
        water += fraction * cw * 2 * receiver_factor
    metadata = {
        "CalibrationIdentifier": CALIBRATION_ID,
        "SignalNormalization": "SRU = receiver_factor * concentration_mol_per_litre * equivalent_protons",
        "CalibrationStatus": "synthetic_internal_consistency_only_not_scanner_calibrated",
        "ConcentrationConvention": "nominal mol/L per total tissue volume; metabolite inputs mM",
        "TissueFractions": fractions,
        "GroundTruthMetaboliteConcentration_mol_L": gt_met,
        "GroundTruthWaterConcentration_mol_L": gt_water,
        "MetaboliteEquivalentProtons": 3,
        "WaterEquivalentProtons": 2,
        "ReceiverFactor": float(receiver_factor),
        "DwellTime_s": float(dwell_time_s),
        "SpectralPoints": int(n_points),
        "ResidualWaterIncluded": False,
        "RelaxationIncluded": False,
    }
    return {"metabolite_fid": met, "water_reference_fid": water,
            "time_s": t, "metadata": metadata}


def export_analytical_pair(pair, output_dir, *, spectrometer_frequency_mhz, affine=None):
    """Write paired 1D NIfTI-MRS files using nifti_mrs library convention."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    affine = np.eye(4) if affine is None else np.asarray(affine)
    paths = {}
    for key, filename in (("metabolite_fid", "metabolites.nii.gz"),
                          ("water_reference_fid", "water_reference.nii.gz")):
        obj = gen_nifti_mrs(pair[key][None, None, None, :, None],
                            pair["metadata"]["DwellTime_s"], spectrometer_frequency_mhz,
                            nucleus="1H", affine=affine, dim_tags=["DIM_DYN", None, None])
        path = output_dir / filename
        obj.save(str(path))
        paths[key] = str(path)
    with open(output_dir / "ground_truth.json", "w") as stream:
        json.dump(pair["metadata"], stream, indent=2)
    return paths
