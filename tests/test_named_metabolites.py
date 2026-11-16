from types import SimpleNamespace
import numpy as np
import pandas as pd
import pytest
import torch
from simulation.signalModel import SignalModel


def model(rows,names=('NAA','Cr','bHG')):
    t=np.arange(256)/2000
    responses={n:np.exp(2j*np.pi*f*t+.2j) for n,f in [('NAA',-300),('Cr',-200),('bHG',-100)]}
    m=SignalModel.__new__(SignalModel);m.labels=['WM'];m.gui=None;m.verbose=False
    m.basis=SimpleNamespace(responses_by_name={n:responses[n] for n in names},n=256,t=t)
    m.metab_df=pd.DataFrame(rows,columns=['Metabolite','Conc_mean','Conc_std','T2']).assign(Tissue='WM')
    return m


@pytest.mark.parametrize('order',[[0,1,2],[2,0,1],[1,2,0]])
def test_named_concentrations_decay_and_appended_bhg(order):
    rows=[('Cr',2.,0.,200.),('NAA',8.,0.,300.),('bHG',4.,0.,np.nan)]
    m=model([rows[i] for i in order]);actual=m.simulate_metabs().numpy().ravel()
    expected=sum(c*m.basis.responses_by_name[n]*np.exp(-m.basis.t/(np.pi*t2*.001)) if np.isfinite(t2) else c*m.basis.responses_by_name[n] for n,c,sd,t2 in rows)
    np.testing.assert_allclose(actual,np.fft.fftshift(np.fft.fft(expected)),rtol=2e-6,atol=2e-5)
    reverse=model([rows[i] for i in order],names=('bHG','Cr','NAA'))
    np.testing.assert_allclose(reverse.simulate_metabs().numpy().ravel(),actual,rtol=2e-6,atol=2e-5)


def test_missing_duplicates_and_seed_draws():
    rows=[('Cr',2.,.1,200.),('NAA',8.,.2,300.)];m=model(rows)
    state=torch.random.get_rng_state()
    try:
        torch.manual_seed(7);expected=torch.normal(torch.tensor([2.,8.],dtype=torch.float64),torch.tensor([.1,.2],dtype=torch.float64))
        torch.manual_seed(7);m.simulate_metabs()
        assert m.sampled_concs_per_tissue['WM']=={'Cr':float(expected[0]),'NAA':float(expected[1]),'bHG':0.}
        assert m.metabolite_mapping_metadata['missing_entries_zero_filled']=={'WM':['bHG']}
        duplicate=model(rows+[rows[0]])
        before=torch.random.get_rng_state().clone()
        with pytest.raises(ValueError,match='duplicate'):duplicate.simulate_metabs()
        assert torch.equal(before,torch.random.get_rng_state())
    finally:torch.random.set_rng_state(state)


from test_basis_h2o import basis


def test_named_basis_reset_and_alias_collision(basis):
    from simulation.basis import Basis
    path,_=basis
    b=Basis(str(path),bw=2000,points=256,metabs=['NAA'])
    before=b.responses_by_name['NAA'].copy()
    assert b.names==['NAA']
    b.reset_basis(metabs=['NAA','Cr'])
    assert set(b.responses_by_name)=={'NAA','Cr'}
    assert b.fids.shape==(256,2)
    for i,n in enumerate(b.names):np.testing.assert_array_equal(b.fids[:,i],b.responses_by_name[n])
    path.write_text(path.read_text().replace("METABO='GPC'","METABO='Cre'"))
    with pytest.raises(ValueError,match='alias collision'):Basis(str(path),bw=2000,points=256)
