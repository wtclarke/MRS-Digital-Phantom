import numpy as np
from research.validation.headless_water.report import self_check
from research.validation.headless_water.audit import linewidth


def test_biased_response_and_failure_reporting():
    self_check()


def test_measured_linewidth_matches_lorentzian():
    t=np.arange(8192)/2000
    measured=linewidth(np.exp(-np.pi*5*t+2j*np.pi*10*t),1/2000,127.)
    # Magnitude Lorentzian width is sqrt(3) times its absorption width.
    assert abs(measured-np.sqrt(3)*5)<.01


def test_reference_noise_uses_paired_controls():
    import pandas as pd
    from research.validation.headless_water.report import reference_noise_effect
    rows=[]
    for repeat,baseline in enumerate([2.,4.]):
        for stage,value in [('noise',baseline),('water_noise',baseline+1)]:
            rows.append(dict(basis_name='basis',fitter='fsl',metabolite='NAA',repeat=repeat,
                stage=stage,target='baseline',tissue='40/40/20',snr=30,water_snr=50,
                estimate_mM=value,valid=True))
    actual=reference_noise_effect(pd.DataFrame(rows)).iloc[0]
    assert actual.n==2 and actual.mean_delta_mM==1 and actual.sd_delta_mM==0
    assert actual.rmse_delta_mM==1


def test_residual_plot_preserves_pairing_and_error_sign(monkeypatch, tmp_path):
    import pandas as pd
    import matplotlib.pyplot as plt
    from research.validation.headless_water import report
    data=pd.DataFrame(dict(basis_name=['GE_PRESS_TE35']*3,target=['NAA']*3,
        metabolite=['NAA']*3,tissue=['100/0/0']*3,fitter=['fsl']*3,stage=['clean']*3,
        valid=[True,True,False],target_truth_mM=[1.,2.,3.],truth_mM=[1.,2.,3.],
        estimate_mM=[1.2,1.7,99.]))
    data=pd.concat([data,data.assign(fitter='lcmodel')],ignore_index=True)
    saved=[]
    def inspect(fig, root, name):
        ax=fig.axes[0]
        np.testing.assert_array_equal(ax.lines[0].get_xdata(),[1.,2.])
        np.testing.assert_allclose(ax.lines[0].get_ydata(),[.2,-.3])
        assert ax.lines[0].get_marker()=='o'
        assert ax.lines[1].get_marker()=='s'
        assert all(line.get_linestyle()=='None' for line in ax.lines[:2])
        np.testing.assert_array_equal(ax.lines[2].get_ydata(),[0.,0.])
        assert ax.lines[2].get_label()=='Zero error'
        assert 'Signed error' in ax.get_ylabel()
        saved.append(name)
        plt.close(fig)
    monkeypatch.setattr(report,'save_figure',inspect)
    report.response_figure(data,tmp_path,'residual')
    assert saved==['target_residual_GE_PRESS_TE35','target_residual']
