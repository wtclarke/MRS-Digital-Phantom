"""Synthetic supplied BASIS shared by retained simulation tests."""
import numpy as np
import pytest


@pytest.fixture
def basis(tmp_path):
    path=tmp_path/'test.BASIS';n=256;cf=100.;bw=2000;t=np.arange(n)/bw
    signals={m:a*np.exp(-np.pi*3*t+1j*(2*np.pi*f*t+p)) for m,a,f,p in [
        ('NAA',3,-267,.2),('Cr',3,-165,-.1),('GPC',9,-147,.05),('PCh',9,-148,.07),
        ('Glu',2,-233,-.2),('H2O',2,2,.4)]}
    text=" $SEQPAR\n HZPPPM=100, ECHOT=35, SEQ='PRESS'\n $END\n $BASIS1\n BADELT=0.0005, NDATAB=256, FMTBAS='(2E25.16)'\n $END\n"
    for name,x in signals.items():
        text+=f" $BASIS\n METABO='{name}', CONC=1, TRAMP=1, VOLUME=1, ISHIFT=0\n $END\n"
        text+=''.join(f'{z.real:25.16E} {z.imag:25.16E}\n' for z in np.fft.fft(x.conj()))
    path.write_text(text)
    return path,signals

