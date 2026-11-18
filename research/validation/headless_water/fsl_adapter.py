"""Fit saved notebook acquisitions using native FSL-MRS water quantification."""
from contextlib import redirect_stdout, redirect_stderr
import hashlib
import io
import json
from pathlib import Path
import time
import traceback
import warnings

import numpy as np
from fsl_mrs import __version__
from fsl_mrs.utils import fitting, mrs_io, quantify
from fsl_mrs.utils.constants import H2O_MOLARITY


def fit_case(source: Path, basis: Path, output: Path) -> dict:
    """Return voxel-average nominal mM; preserve native tissue-only estimates."""
    source, basis, output = Path(source), Path(basis), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    log = io.StringIO()
    record = {'status': 'failed', 'estimates_mM': {}, 'native_estimates_mM': {},
              'version': __version__, 'diagnostics': {}}
    with redirect_stdout(log), redirect_stderr(log), warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            truth = json.loads((source / 'ground_truth.json').read_text())
            provenance = json.loads((source / 'calibration_provenance.json').read_text())
            if hashlib.sha256(basis.read_bytes()).hexdigest() != provenance['basis_sha256']:
                raise ValueError('fitting BASIS does not match simulated BASIS')
            fid = mrs_io.read_FID(str(source / 'total.nii.gz'))
            water = mrs_io.read_FID(str(source / 'water_reference.nii.gz'))
            if (not np.isclose(fid.dwelltime, water.dwelltime, rtol=1e-7, atol=0)
                    or not np.allclose(fid.spectrometer_frequency, water.spectrometer_frequency)
                    or fid.SpecFreqChemShift != water.SpecFreqChemShift):
                raise ValueError('metabolite/water sampling or frequency metadata differ')
            mrs = fid.mrs(basis=mrs_io.read_basis(str(basis)), ref_data=water)
            # BASIS readers and NIfTI reader already define the common complex convention.
            mrs.ignore = [n for n in mrs.names if n in ('H2O', 'Ref0ppm')]
            mrs.rescaleForFitting()
            result = fitting.fit_FSLModel(
                mrs, method='Newton', model='voigt', baseline='off',
                ppmlim=(0.2, 4.2), metab_groups=[0] * len(mrs.names),
                capture_minimize_output=True)
            fractions = {t: float(truth['tissue_fractions'].get(t, 0.)) for t in ('WM', 'GM', 'CSF')}
            if not np.isclose(sum(fractions.values()), 1., atol=1e-6):
                raise ValueError('water quantification requires WM/GM/CSF fractions summing to one')
            if fractions['CSF'] >= 1.:
                raise ValueError('pure CSF has no tissue-only concentration or reference metabolite')
            # PCr is unswept, so zero-NAA/Cr cases retain a nonzero methyl reference.
            # Quantification is unrelaxed; intrinsic sequence evolution in BASIS is retained.
            qinfo = quantify.QuantificationInfo(
                0., np.inf, mrs.names, mrs.centralFrequency / 1e6,
                water_ref_metab='PCr', water_ref_metab_protons=3,
                water_ref_metab_limits=(2.9, 3.15),
                t1=dict.fromkeys(('H2O_WM', 'H2O_GM', 'H2O_CSF', 'METAB'), 1.),
                t2=dict.fromkeys(('H2O_WM', 'H2O_GM', 'H2O_CSF', 'METAB'), 1.))
            qinfo.set_fractions(fractions)
            profiles = provenance['water_options']['water_concentrations_mol_l']
            qinfo.set_densities({t: 1000 * profiles[t] / H2O_MOLARITY for t in fractions})
            result.calculateConcScaling(mrs, quant_info=qinfo, internal_reference='PCr')
            native = result.getConc(scaling='molarity')
            record['native_estimates_mM'] = dict(zip(result.metabs, map(float, native)))
            record['estimates_mM'] = {n: c * (1 - fractions['CSF']) for n, c in record['native_estimates_mM'].items()}
            if not all(np.isfinite(list(record['estimates_mM'].values()))):
                raise ValueError('nonfinite water-scaled concentrations')
            for what in ('concentrations', 'parameters', 'qc', 'summary'):
                result.to_file(output / f'{what}.csv', what=what)
            qinfo.summary_table.to_csv(output / 'quantification_info.csv')
            np.savez(output / 'fit_arrays.npz', ppm=mrs.getAxes('ppm'),
                     data=mrs.FID, predicted=result.pred, residual=result.residuals)
            opt = result.scipy_minimize_result
            (output / 'optimizer.txt').write_text(str(opt) + '\n')
            record['diagnostics'] = {
                'optimizer_success': bool(opt.success), 'optimizer_message': str(opt.message),
                'optimizer_status': int(opt.status),
                'relative_residual_norm': float(np.linalg.norm(result.residuals) / np.linalg.norm(mrs.FID)),
                'native_molarity_to_voxel_multiplier': 1 - fractions['CSF'],
                'water_molarity_mM': qinfo.relax_corr_water_molar,
                'raw_concentration_scaling': float(result.concScalings['raw']),
                'molarity_scaling': float(result.concScalings['molarity']),
                'water_ref_metab': 'PCr', 'water_ref_protons': 3,
                'water_ref_limits_ppm': [2.9, 3.15],
                'water_reference_reason': 'Fixed unswept PCr methyl reference remains present when NAA or Cr is zero',
                'quantification_relaxation': 'none; TE=0/TR=infinity limit',
                'model': 'voigt', 'baseline': 'off', 'fit_ppm': [0.2, 4.2]}
            record['status'] = 'ok' if opt.success else 'optimizer_warning'
            _plot(mrs, result, output)
        except Exception as exc:
            record['status'] = 'failed'
            record['estimates_mM'] = {}
            record['error'] = str(exc)
            traceback.print_exc()
        record['warnings'] = [str(w.message) for w in caught]
    record['elapsed_s'] = time.monotonic() - started
    (output / 'execution.log').write_text(log.getvalue())
    (output / 'result.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
    return record


def _plot(mrs, result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from fsl_mrs.utils.misc import FIDToSpec
    ppm = mrs.getAxes('ppm')
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(9, 5), height_ratios=(3, 1))
    axes[0].plot(ppm, FIDToSpec(mrs.FID).real, label='Data')
    axes[0].plot(ppm, result.pred_spec.real, label='Fit')
    axes[0].legend()
    axes[1].plot(ppm, FIDToSpec(result.residuals).real)
    axes[1].set(xlim=(4.2, .2), xlabel='Chemical shift (ppm)', ylabel='Residual')
    axes[0].set_ylabel('FSL fitting signal units')
    fig.tight_layout()
    for suffix in ('png', 'pdf'):
        fig.savefig(output / f'fit.{suffix}')
    plt.close(fig)
