import numpy as np
from research.validation.headless_water.run import subject,concentration_table


def test_controlled_fractions_field_and_concentration_scaling():
    case=dict(wm=.4,gm=.4,csf=.2,b0_sd=5)
    s=subject(case);labels=s.labels.data.numpy()[0,1:11,1:11,1:11]
    np.testing.assert_allclose([(labels==i).mean() for i in [1,2,3]],[.4,.4,.2])
    assert len(np.unique(s.labels.data))==6
    assert np.isclose(s.B0.data.numpy()[0,1:11,1:11,1:11].std(),5)
    base=concentration_table(['NAA','Cr','bHG','Ins'],'baseline',1)
    doubled=concentration_table(['NAA','Cr','bHG','Ins'],'bHG',2)
    assert (base.Conc_std==0).all() and base.T2.isna().all()
    np.testing.assert_array_equal(doubled.loc[doubled.Metabolite=='bHG','Conc_mean'],2*base.loc[base.Metabolite=='bHG','Conc_mean'])
    assert 'mI' in set(base.Metabolite)
