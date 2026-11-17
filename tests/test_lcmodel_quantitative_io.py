"""Fast conversion tests: no proprietary executable required."""
import hashlib
from pathlib import Path

import numpy as np
import pytest

from utils.lcmodel_quantitative_io import (
    NAMES, compatible_basis, control_text, parse_basis, parse_coord, parse_raw, parse_table, write_raw,
)



def test_complex_raw_conjugation_precision_frequency_and_phase(tmp_path):
    t=np.arange(256)*.0005
    fid=1.234567891234567e-4*np.exp(-np.pi*3*t)*np.exp(1j*(2*np.pi*37*t+np.deg2rad(35)))
    first=tmp_path/'one.RAW';second=tmp_path/'two.RAW'
    for path in [first,second]:write_raw(path,fid,frequency_mhz=127.731)
    stored,meta=parse_raw(first)
    np.testing.assert_array_equal(stored,fid.conj())
    from fsl_mrs.utils.mrs_io.lcm_io import readLCModelRaw
    interpreted,_=readLCModelRaw(first,unpack_header=False)
    np.testing.assert_array_equal(interpreted,fid)
    assert first.read_bytes()==second.read_bytes()
    assert meta['frequency_mhz']==127.731
    assert np.angle(stored[0])==pytest.approx(-np.deg2rad(35))
    assert np.angle(stored[1]/stored[0])/(2*np.pi*.0005)==pytest.approx(-37)
    assert 'FMTDAT=\'(2E25.16)\'' in meta['header']
    with pytest.raises(FileExistsError):write_raw(first,fid,frequency_mhz=127.731)


def test_relative_water_amplitude_preserved(tmp_path):
    t=np.arange(256)*.0005
    metabolite=.003*np.exp(1j*(2*np.pi*37*t+.4))*np.exp(-np.pi*t)
    water=80*np.exp(1j*(2*np.pi*2.55462*t+.2))*np.exp(-np.pi*t)
    for name,fid in [('met',metabolite),('water',water)]:write_raw(tmp_path/name,fid,frequency_mhz=127.731)
    a,_=parse_raw(tmp_path/'met');b,_=parse_raw(tmp_path/'water')
    np.testing.assert_array_equal((a/b).conj(),metabolite/water)


@pytest.mark.parametrize('bad',[np.ones(8),np.ones((2,8),complex),np.array([complex(np.nan,0)])])
def test_invalid_raw_rejected(tmp_path,bad):
    with pytest.raises(ValueError):write_raw(tmp_path/'bad',bad,frequency_mhz=127.731)


def test_basis_fixed_format_reference_cancellation_and_source_immutability(tmp_path,basis):
    BASIS=basis[0]
    before=hashlib.sha256(BASIS.read_bytes()).hexdigest()
    target=tmp_path/'converted.BASIS'
    meta=compatible_basis(BASIS,target,data_origin_ppm=5.275)
    original,original_meta=parse_basis(BASIS);copy,copy_meta=parse_basis(target)
    assert copy_meta==original_meta
    assert len(copy)==6
    assert meta['added_ISHIFT']==8 and meta['net_LCModel_referencing_shift']==0
    assert meta['spectrum_max_error']==0
    assert "FMTBAS='(2E25.16)'" in target.read_text()
    for a,b in zip(original,copy):
        assert a[0]==b[0]
        np.testing.assert_array_equal(a[2],b[2])
        assert 'ISHIFT=8' in b[1].replace(' ','')
    assert hashlib.sha256(BASIS.read_bytes()).hexdigest()==before


def test_control_sampling_water_assumptions_and_explicit_selection(tmp_path):
    names=list(NAMES)+['H2O','Ref0ppm','Ins']
    text=control_text(tmp_path,tmp_path/'basis',names=names,points=2048,dwell_time_s=.0005,
        frequency_mhz=127.731,origin_ppm=4.68,water_mol_l=39.0165)
    for setting in ['NUNFIL=2048','DELTAT=0.0005','HZPPPM=127.731','WCONC=39016.5',
                    'ATTH2O=1.','ATTMET=1.',"WSMET='NAA'",'N1HMET=3','DOECC=F','NOMIT=3',"CHOMIT(1)='H2O'","CHOMIT(2)='Ref0ppm'",'NSIMUL=0','NRATIO=0']:
        assert setting in text
    with pytest.raises(ValueError,match='missing'):
        control_text(tmp_path,tmp_path/'basis',names=['Cr'],points=2048,dwell_time_s=.0005,
            frequency_mhz=127.731,origin_ppm=4.68,water_mol_l=39.0165)


def test_table_parser_detects_fatal_even_with_a_concentration(tmp_path):
    table=tmp_path/'test.TABLE'
    table.write_text(' (Version 6.3-1N)\n 1.0639E+01 0% 1.260 NAA\n 1 FATAL MYBASI 31\n')
    result=parse_table(table)
    assert result['version']=='6.3-1N' and result['fatal']
    assert result['concentrations']['NAA']['concentration']==10.639
    assert result['diagnostics'][0]['code']==31


def test_coord_parser_keeps_real_data_and_baseline_distinct(tmp_path):
    path=tmp_path/'test.COORD'
    path.write_text('3 points on ppm-axis = NY\n 3 2 1\nNY phased data points follow\n 4 5 6\n'
                    'NY points of the fit to the data follow\n 3 4 5\nNY background values follow\n .1 .2 .3\n')
    result=parse_coord(path)
    np.testing.assert_array_equal(result['data']-result['fit'],np.ones(3))
    np.testing.assert_array_equal(result['baseline'],[.1,.2,.3])


def test_basis_with_interleaved_nmused_and_existing_format(tmp_path):
    source=tmp_path/'GE.BASIS'
    source.write_text(" $SEQPAR\n HZPPPM=127.731, ECHOT=35\n $END\n $BASIS1\n NDATAB=2, BADELT=.00025, FMTBAS='(2E15.6)'\n $END\n"+''.join(f" $NMUSED\n XTRASH=0\n $END\n $BASIS\n METABO='{n}', ISHIFT=0\n $END\n 1 2\n 3 4\n" for n in ['NAA','H2O']))
    entries,_=parse_basis(source)
    assert [e[0] for e in entries]==['NAA','H2O']
    target=tmp_path/'converted.BASIS';compatible_basis(source,target,data_origin_ppm=4.68)
    assert target.read_text().count('FMTBAS')==1
    for a,b in zip(entries,parse_basis(target)[0]):np.testing.assert_array_equal(a[2],b[2])
