"""Headless helpers for generating spectra from the MRiLab phantom."""

import os
import math
import json
from numbers import Integral

import numpy as np
import torch
from nifti_mrs.create_nmrs import gen_nifti_mrs

from gui.defaults import DEFAULT_METAB_DF_PATH
from simulation.mrs_phantom import DigitalPhantom
from utils.auxillary import spec2fid
from utils.nii_processing import save_ground_truth_concentrations, save_nifti_mrs


def _position_and_bounds(position, voxel_counts, volume_shape):
    position = tuple(position)
    if len(position) != 3 or any(not isinstance(value, Integral) for value in position):
        raise ValueError("Each position must be a three-integer (x, y, z) tuple.")
    position = tuple(int(value) for value in position)
    starts = tuple(center - count // 2 for center, count in zip(position, voxel_counts))
    ends = tuple(start + count for start, count in zip(starts, voxel_counts))
    if any(start < 0 or end > length for start, end, length in zip(starts, ends, volume_shape)):
        raise ValueError(f"Voxel centered at {position} extends beyond volume bounds {volume_shape}.")
    bounds = (starts[0], ends[0], starts[1], ends[1], starts[2], ends[2])
    return position, bounds


def _voxel_affine(affine, position, voxel_size_mm):
    affine = np.array(affine, dtype=float, copy=True)
    directions = affine[:3, :3]
    center_world = directions @ np.asarray(position) + affine[:3, 3]
    spacing = np.linalg.norm(directions, axis=0)
    if np.any(spacing == 0):
        raise ValueError("Phantom affine must have non-zero spatial axes.")
    affine[:3, :3] = directions / spacing * np.asarray(voxel_size_mm)
    affine[:3, 3] = center_world
    return affine


def simulate_batch(
    positions,
    settings,
    save_dir=None,
    total_components=None,
    individual_components=(),
    seed=None,
):
    """Simulate MRiLab spectra at voxel-center indices, reusing one phantom.

    ``settings`` contains ``path2basis``, ``simulation_params``,
    ``basis_set_settings`` (vendor, localization, and TE), and ``voxel_size_mm``.
    Saving is enabled by providing ``save_dir``. The result is one dictionary per
    requested position, including the in-memory spectra and concentration ground truth.
    """
    basis_path = os.fspath(settings["path2basis"])
    if not os.path.isfile(basis_path):
        raise FileNotFoundError(f"Basis set not found: {basis_path}")

    simulation_settings = dict(settings["simulation_params"])
    basis_settings = dict(settings["basis_set_settings"])
    voxel_size_mm = tuple(float(value) for value in settings["voxel_size_mm"])
    if len(voxel_size_mm) != 3 or any(not math.isfinite(value) or value <= 0 for value in voxel_size_mm):
        raise ValueError("voxel_size_mm must contain three positive values.")
    if seed is not None and (isinstance(seed, bool) or not isinstance(seed, Integral) or not 0 <= seed <= 2**32 - 1):
        raise ValueError("seed must be an integer between 0 and 2**32 - 1.")
    positions = list(positions)
    normalized_positions = []
    for position in positions:
        try:
            position = tuple(position)
        except TypeError as error:
            raise ValueError("Each position must be a three-integer (x, y, z) tuple.") from error
        if len(position) != 3 or any(not isinstance(value, Integral) for value in position):
            raise ValueError("Each position must be a three-integer (x, y, z) tuple.")
        normalized_positions.append(tuple(int(value) for value in position))
    positions = normalized_positions
    if save_dir is not None and not positions:
        raise ValueError("At least one position is required when saving NIfTI-MRS files.")
    if not positions:
        return []

    if save_dir is not None:
        total = ["Metabs", "MMs", "Lipids", "Water", "Noise"] if total_components is None else list(total_components)
        individual = list(individual_components)
        unknown = (set(total) | set(individual)) - {"Metabs", "MMs", "Lipids", "Water", "Noise"}
        if unknown:
            raise ValueError(f"Unknown component names: {sorted(unknown)}")
        if not total and not individual:
            raise ValueError("Select at least one component to save.")
        save_dir = os.fspath(save_dir)
        if os.path.exists(save_dir) and not os.path.isdir(save_dir):
            raise NotADirectoryError(f"NIfTI output path is not a directory: {save_dir}")

    simulation_params = {**simulation_settings, **basis_settings}
    phantom = DigitalPhantom(
        skeleton="MRiLab",
        path2metabs=DEFAULT_METAB_DF_PATH,
        path2basis=basis_path,
        simulation_params=simulation_params,
        metabs=[],
        verbose=False,
    )
    voxel_counts = tuple(round(size / spacing) for size, spacing in zip(voxel_size_mm, phantom.spacing))
    if any(count < 1 for count in voxel_counts):
        raise ValueError("Voxel size must cover at least one phantom voxel on each axis.")

    volume_shape = phantom.get_label_data().shape
    position_bounds = [
        _position_and_bounds(position, voxel_counts, volume_shape)
        for position in positions
    ]
    if save_dir is not None:
        output_names = [
            f"voxel_{index:04d}_{position[0]}_{position[1]}_{position[2]}"
            for index, (position, _) in enumerate(position_bounds)
        ]
        existing = [name for name in output_names if os.path.exists(os.path.join(save_dir, name))]
        if existing:
            raise FileExistsError(f"NIfTI output already exists: {os.path.join(save_dir, existing[0])}")
        os.makedirs(save_dir, exist_ok=True)

    numpy_state = np.random.get_state() if seed is not None else None
    torch_state = torch.random.get_rng_state() if seed is not None else None
    if seed is not None:
        np.random.seed(int(seed))
        torch.manual_seed(int(seed))

    try:
        results = []
        for index, (position, bounds) in enumerate(position_bounds):
            spectrum, components, time_axis, ppm_axis, _ = phantom.simulate_data(bounds)
            model = phantom.signal_model
            ground_truth = {
                "sampled_per_tissue": model.sampled_concs_per_tissue,
                "voxel_weighted": model.final_voxel_concentrations,
                "tissue_fractions": model.voxel_tissue_fractions,
            }
            sim_params = {
                "skeleton": "MRiLab",
                "path2metabs": DEFAULT_METAB_DF_PATH,
                "path2basis": basis_path,
                "metabs": phantom.metabs,
                "basis_set_settings": basis_settings,
                "simulation_params": simulation_settings,
                "seed": seed,
                "voxel_definitions": {
                    "coords": bounds,
                    "size": voxel_counts,
                    "size_mm": voxel_size_mm,
                },
            }
            nifti_metadata = {
                "affine": _voxel_affine(phantom.affine, position, voxel_size_mm).tolist(),
                "dwelltime": float(phantom.basis_set.dwelltime),
                "spec_freq": float(phantom.basis_set.cf),
                "nucleus": "1H",
                "sim_params": sim_params,
            }
            result = {
                "position": position,
                "spectrum": spectrum,
                "components": components,
                "time_axis": time_axis,
                "ppm_axis": ppm_axis,
                "ground_truth": ground_truth,
                "seed": seed,
                "nifti_metadata": nifti_metadata,
            }

            if save_dir is not None:
                base_name = output_names[index]
                output_path = os.path.join(save_dir, base_name)
                save_nifti_mrs(
                    save_dir,
                    total,
                    individual,
                    components,
                    _voxel_affine(phantom.affine, position, voxel_size_mm),
                    phantom.basis_set.dwelltime,
                    phantom.basis_set.cf,
                    sim_params=sim_params,
                    base_name=base_name,
                    sampled_per_tissue=ground_truth["sampled_per_tissue"],
                    voxel_weighted=ground_truth["voxel_weighted"],
                    tissue_fractions=ground_truth["tissue_fractions"],
                )
                result["saved_path"] = output_path

            results.append(result)
        return results
    finally:
        if seed is not None:
            np.random.set_state(numpy_state)
            torch.random.set_rng_state(torch_state)


def _as_nifti_mrs(spectrum, metadata):
    fid = spec2fid(spectrum, axis=-1)
    if isinstance(fid, torch.Tensor):
        fid = fid.detach().cpu().numpy()
    data = np.asarray(fid)[None, None, None, :, None].astype(np.complex64)
    nifti = gen_nifti_mrs(
        data,
        metadata["dwelltime"],
        metadata["spec_freq"],
        nucleus=metadata["nucleus"],
        affine=np.asarray(metadata["affine"], dtype=float),
        dim_tags=["DIM_DYN", None, None],
    )
    sim_params = metadata["sim_params"]
    nifti.add_hdr_field("EchoTime", sim_params["basis_set_settings"]["TE"] * 1e-3)
    nifti.add_hdr_field("RepetitionTime", sim_params["simulation_params"]["TR"] * 1e-3)
    nifti.add_hdr_field("Manufacturer", sim_params["basis_set_settings"]["vendor"])
    return nifti


def results_to_nifti_mrs(
    results,
    output_dir=None,
    total_components=None,
    individual_components=(),
):
    """Convert simulation results to plot-ready NIFTI_MRS objects.

    Returns one dictionary per voxel with a ``total`` object and selected
    ``components``. Supplying ``output_dir`` also saves these objects and the
    simulation/ground-truth sidecars; returned objects remain available to plot.
    """
    results = list(results)
    if not results:
        return []

    component_names = list(results[0]["components"])
    total = component_names if total_components is None else list(total_components)
    individual = list(individual_components)
    if not total and not individual:
        raise ValueError("Select at least one component to convert.")
    for result in results:
        unknown = (set(total) | set(individual)) - set(result["components"])
        if unknown:
            raise ValueError(f"Unknown component names: {sorted(unknown)}")

    output_dir = None if output_dir is None else os.fspath(output_dir)
    if output_dir is not None:
        if not output_dir:
            raise ValueError("output_dir must be a non-empty path.")
        if os.path.exists(output_dir) and not os.path.isdir(output_dir):
            raise NotADirectoryError(f"NIfTI output path is not a directory: {output_dir}")
        output_names = [
            f"voxel_{index:04d}_{result['position'][0]}_{result['position'][1]}_{result['position'][2]}"
            for index, result in enumerate(results)
        ]
        existing = [name for name in output_names if os.path.exists(os.path.join(output_dir, name))]
        if existing:
            raise FileExistsError(f"NIfTI output already exists: {os.path.join(output_dir, existing[0])}")
        os.makedirs(output_dir, exist_ok=True)

    converted = []
    for index, result in enumerate(results):
        components = result["components"]
        metadata = result["nifti_metadata"]
        total_object = _as_nifti_mrs(
            sum((components[name] for name in total), start=torch.zeros_like(result["spectrum"]))
            if total else result["spectrum"],
            metadata,
        ) if total else None
        component_objects = {
            name: _as_nifti_mrs(components[name], metadata)
            for name in individual
        }
        saved_paths = {}

        if output_dir is not None:
            voxel_dir = os.path.join(output_dir, output_names[index])
            os.makedirs(voxel_dir)
            if total_object is not None:
                path = os.path.join(voxel_dir, "total.nii.gz")
                total_object.save(path)
                saved_paths["total"] = path
            component_paths = {}
            for name, nifti in component_objects.items():
                path = os.path.join(voxel_dir, f"{name.replace(' ', '_').lower()}.nii.gz")
                nifti.save(path)
                component_paths[name] = path
            saved_paths["components"] = component_paths
            with open(os.path.join(voxel_dir, "sim_params.json"), "w") as stream:
                json.dump(metadata["sim_params"], stream, indent=2)
            ground_truth = result.get("ground_truth", {})
            save_ground_truth_concentrations(
                voxel_dir,
                sampled_per_tissue=ground_truth.get("sampled_per_tissue"),
                voxel_weighted=ground_truth.get("voxel_weighted"),
                tissue_fractions=ground_truth.get("tissue_fractions"),
            )

        converted.append({
            "position": result["position"],
            "total": total_object,
            "components": component_objects,
            "saved_paths": saved_paths,
        })
    return converted
