from spectrum import *
from spectrum.mtm import dpss, pmtm
from spectrum import data_cosine

def test_dpss():
    dpss(64, 2.5, 4)

def test_pmtm():
    data = data_cosine(N=64, A=0.1, sampling=1024, freq=200)
    res = pmtm(data, 2.5, 4, show=False)
    res = pmtm(data, 2.5, show=False)
