"""Opt-in normalized anatomical-cell integration on the Phase 6a signal scale."""
from numbers import Integral

import numpy as np

from simulation.mrscloud_full_basis import load_library
from simulation.mrscloud_quantitative_voxel import simulate_voxel, TISSUES
from simulation.mrscloud_quantitative_noise import add_pair_noise


def simulate_spatial(cell_labels, *, calibration_bundle, metabolite_concentrations_mM,
                     water_concentrations_mol_l, b0_offsets_hz=None, spatial_weights=None,
                     cell_volumes_mm3=None, relaxation=None, receiver_factor=1.,
                     bandwidth_hz=2000, points=2048, expected_calibration=None, noise=None):
    """Integrate a selected 1D cell list or 3D ROI, in C-order.

    Labels are explicit WM/GM/CSF names. Weights must already sum to one;
    volumes are instead normalized by their total, with no volume gain.
    Positive Hz offsets multiply both acquisitions by exp(+2j*pi*offset*t).
    """
    labels=np.asarray(cell_labels,dtype=str)
    if labels.ndim not in (1,3) or not labels.size or not np.all(np.isin(labels,TISSUES)):
        raise ValueError('cell_labels must be a nonempty 1D/3D array of WM, GM, CSF names')
    shape=labels.shape
    def grid(values, name):
        arr=np.asarray(values,dtype=float)
        if arr.shape!=shape or not np.all(np.isfinite(arr)):
            raise ValueError(f'{name} must be finite with the exact selected ROI shape {shape}')
        return arr.reshape(-1).copy()
    labels=labels.reshape(-1)
    offsets=np.zeros(labels.size) if b0_offsets_hz is None else grid(b0_offsets_hz,'b0_offsets_hz')
    if spatial_weights is not None and cell_volumes_mm3 is not None:
        raise ValueError('provide normalized spatial_weights or cell_volumes_mm3, not both')
    if cell_volumes_mm3 is not None:
        volumes=grid(cell_volumes_mm3,'cell_volumes_mm3')
        if np.any(volumes<=0):raise ValueError('cell volumes must be positive mm^3')
        weights=volumes.copy();mode='cell_volume';volume=float(np.sum(volumes))
        if not np.isfinite(volume):raise ValueError('total cell volume must be finite')
    elif spatial_weights is not None:
        weights=grid(spatial_weights,'spatial_weights');mode='explicit_normalized';volume=None
        if np.any(weights<0) or not np.isclose(np.sum(weights),1.,rtol=0,atol=1e-8):
            raise ValueError('spatial weights must be nonnegative and sum to one (absolute tolerance 1e-8)')
    else:
        weights=np.full(labels.size,1/labels.size);mode='equal_cell';volume=None
    divisor=float(np.sum(weights));weights/=divisor
    fractions={t:float(np.sum(weights[labels==t])) for t in TISSUES}
    pair=simulate_voxel(fractions,calibration_bundle=calibration_bundle,
        metabolite_concentrations_mM=metabolite_concentrations_mM,
        water_concentrations_mol_l=water_concentrations_mol_l,relaxation=relaxation,
        receiver_factor=receiver_factor,bandwidth_hz=bandwidth_hz,points=points,
        expected_calibration=expected_calibration)
    names=list(pair['components'])
    responses,water,_=load_library(calibration_bundle,selected_names=names,
        bandwidth_hz=bandwidth_hz,points=points,expected_calibration=pair['simulation_config']['synthetic_calibration'])
    gt=pair['ground_truth'];config=pair['simulation_config'];time=pair['time_s']
    envelopes={t:np.zeros(points,complex) for t in TISSUES}
    # Aggregate identical offsets; avoid an anatomical-cells by time-points allocation.
    for tissue in TISSUES:
        mask=labels==tissue
        frequencies,inverse=np.unique(offsets[mask],return_inverse=True)
        grouped=np.bincount(inverse,weights=weights[mask])
        for frequency,weight in zip(frequencies,grouped):
            envelopes[tissue]+=weight*np.exp(2j*np.pi*frequency*time)
    relax=gt.get('relaxation')
    components={}
    for name in names:
        coefficient=sum(envelopes[t]*gt['tissue_specific_metabolite_concentrations_mM'][t][name]*
            (relax['metabolite_factors'][t][name] if relax else 1.) for t in TISSUES)
        components[name]=float(receiver_factor)*.001*responses[name]*coefficient
    water_coefficient=sum(envelopes[t]*gt['tissue_specific_water_concentrations_mol_l'][t]*
        (relax['water_factors'][t] if relax else 1.) for t in TISSUES)
    pair.update(components=components,metabolite_fid=sum(components.values()),
                water_reference_fid=float(receiver_factor)*water*water_coefficient)
    mean=float(np.dot(weights,offsets))
    spatial={'integration':'normalized cell-weighted complex sum; no anatomical volume multiplier',
        'roi_shape':list(shape),'flattening_order':'C','cell_count':int(labels.size),
        'cell_labels':labels.tolist(),'cell_counts_by_tissue':{t:int(np.sum(labels==t)) for t in TISSUES},
        'normalized_weights':weights.tolist(),'weighting_mode':mode,'input_weight_sum':divisor,
        'weight_concentration_convention':'normalized weights define represented volume quadrature; no coil sensitivity model',
        'cell_count_fractions':{t:float(np.mean(labels==t)) for t in TISSUES},
        'applied_weight_divisor':divisor,'explicit_weight_sum_tolerance':1e-8,
        'total_anatomical_volume_mm3':volume,
        'cell_volumes_mm3':volumes.tolist() if cell_volumes_mm3 is not None else None,
        'b0_offsets_hz':offsets.tolist(),'b0_units':'Hz','b0_modulation':'exp(+2j*pi*offset_hz*time_s)',
        'b0_source':'explicit_Hz_offsets' if b0_offsets_hz is not None else 'zero',
        'b0_weighted_mean_hz':mean,'b0_weighted_std_hz':float(np.sqrt(np.dot(weights,(offsets-mean)**2))),
        'b0_min_hz':float(np.min(offsets)),'b0_max_hz':float(np.max(offsets)),
        'shared_metabolite_water_offsets':True}
    gt['spatial_integration']=spatial
    # Always expose effective signal concentrations, including the R=1 spatial case.
    if relax is None:
        gt['relaxation_weighted_metabolite_concentrations_mM']=dict(gt['voxel_average_metabolite_concentrations_mM'])
        gt['relaxation_weighted_water_concentration_mol_l']=gt['effective_voxel_water_concentration_mol_l']
    config.update(pathway='synthetic_mrscloud_quantitative_spatial_voxel',spatial_simulation=True,
        B0_inhomogeneity=bool(np.any(offsets!=0)),spatial_integration=spatial,
        water_frequency_offset_hz_scope='intrinsic unit-response offset before cell B0 modulation')
    return add_pair_noise(pair, **noise) if noise is not None else pair


def simulate_phantom_spatial(phantom, voi_bounds, *, use_phantom_b0=False,
                             b0_units=None, b0_convention=None, spatial_weights=None, **kwargs):
    """Crop an existing phantom. Loading B0 requires explicit units and sign.

    b0_convention is positive_frequency or legacy_shim; the latter negates
    values to reproduce the legacy exp(-i*2*pi*map*t) modulation.
    """
    labels=np.asarray(phantom.get_label_data())
    if labels.ndim!=3 or not np.all(np.isfinite(labels)) or not np.all(labels==np.floor(labels)):
        raise ValueError('phantom labels must be a finite integer-coded 3D map')
    bounds=tuple(voi_bounds)
    if len(bounds)!=6 or any(not isinstance(v,Integral) or isinstance(v,bool) for v in bounds):
        raise ValueError('VOI bounds require six integer start/stop indices')
    if any(a<0 or b<=a or b>n for a,b,n in zip(bounds[::2],bounds[1::2],labels.shape)):
        raise ValueError('VOI must be nonempty and inside the phantom label map')
    spacing=np.asarray(phantom.spacing,dtype=float);affine=np.asarray(phantom.affine,dtype=float)
    if (spacing.shape!=(3,) or np.any(spacing<=0) or not np.all(np.isfinite(spacing))
            or affine.shape!=(4,4) or not np.all(np.isfinite(affine))
            or not np.allclose(affine[3],[0,0,0,1],rtol=0,atol=1e-12)
            or not np.allclose(np.linalg.norm(affine[:3,:3],axis=0),spacing,rtol=1e-6,atol=1e-8)):
        raise ValueError('phantom spacing and finite homogeneous affine must agree')
    volume=float(abs(np.linalg.det(affine[:3,:3])))
    if not np.isfinite(volume) or volume<=0:raise ValueError('phantom affine must have nonzero cell volume')
    if not np.allclose(phantom.subject['labels'].affine,affine,rtol=0,atol=1e-6):
        raise ValueError('phantom and label-image affines disagree')
    region=tuple(slice(a,b) for a,b in zip(bounds[::2],bounds[1::2]))
    selected=labels[region];names=list(phantom.label_names)
    if len(set(names))!=len(names) or not set(TISSUES)<=set(names):
        raise ValueError('phantom label names must identify WM, GM and CSF unambiguously')
    ids={names.index(t):t for t in TISSUES}
    if not np.all(np.isin(selected,list(ids))):raise ValueError('ROI contains background or unknown tissue')
    cell_labels=np.empty(selected.shape,dtype='<U3')
    for index,tissue in ids.items():cell_labels[selected==index]=tissue
    offsets=None
    if use_phantom_b0:
        if b0_units!='Hz' or b0_convention not in ('positive_frequency','legacy_shim'):
            raise ValueError('phantom B0 requires b0_units=Hz and explicit positive_frequency or legacy_shim convention')
        if 'B0' not in phantom.subject:raise ValueError('phantom has no B0 map')
        if not np.allclose(phantom.subject['B0'].affine,affine,rtol=0,atol=1e-6):
            raise ValueError('B0 and label-map affines disagree; resampling is not implicit')
        b0=np.asarray(phantom.get_B0_data(),dtype=float)
        if b0.shape!=labels.shape:raise ValueError('B0 and label-map shapes disagree')
        offsets=b0[region].copy()
        if b0_convention=='legacy_shim':offsets*=-1
    elif b0_units is not None or b0_convention is not None:
        raise ValueError('B0 unit/sign declarations require use_phantom_b0=True')
    pair=simulate_spatial(cell_labels,b0_offsets_hz=offsets,spatial_weights=spatial_weights,
        cell_volumes_mm3=np.full(selected.shape,volume) if spatial_weights is None else None,**kwargs)
    pair['ground_truth']['voxel_geometry']={'exclusive_stop_bounds_index':[int(v) for v in bounds],
        'anatomical_spacing_mm':spacing.tolist(),'anatomical_affine':affine.tolist(),
        'cell_volume_mm3':volume,'cell_volume_definition':'abs(det(affine[:3,:3]))',
        'voxel_volume_mm3':volume*selected.size}
    spatial=pair['ground_truth']['spatial_integration']
    spatial.update(b0_source='phantom.get_B0_data' if use_phantom_b0 else 'zero; phantom B0 not loaded',
                   input_b0_units=b0_units,input_b0_convention=b0_convention,
                   b0_sign_conversion=-1 if b0_convention=='legacy_shim' else 1)
    return pair
