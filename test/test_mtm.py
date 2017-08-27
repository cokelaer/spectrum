from spectrum import *
from spectrum.mtm import dpss, pmtm
from spectrum import data_cosine
from spectrum import mtm

def test_dpss():
    dpss(64, 2.5, 4)


def test_pmtm():
    data = data_cosine(N=64, A=0.1, sampling=1024, freq=200)
    res = pmtm(data, 2.5, 4, show=False)
    res = pmtm(data, 2.5, show=False)



    res = pmtm(data, 2.5, show=False, method="eigen")
    res = pmtm(data, 2.5, show=False, method="unity")
    res = pmtm(data, 2.5, method="eigen", show=True)
    res = pmtm(data, 2.5, method="adapt", show=True)
    #res = pmtm(data, 2.5, show=False, method="eigen", show=True)

    # e and v must be provided together
    try:
        res = pmtm(data, 2.5, show=False, e=1, v=None)
        assert False
    except:
        assert True

    # provide v and e
    v,e = dpss(64,4,2)
    pmtm(marple_data, NW=4, k=2, v=v, e=e);


    try:
        pmtm(marple_data, NW=None, k=2);
        assert False
    except:
        assert True


def test_fftconvolve():
    from spectrum import mtm
    from pylab import randn
    mtm._fftconvolve(randn(128), randn(128),mode="full")
    mtm._fftconvolve(randn(128), randn(128), mode="same")
    mtm._fftconvolve(randn(128), randn(128), mode="valid")


def test_crosscorr():
    from spectrum import mtm
    from pylab import randn
    mtm._crosscorr(randn(128), randn(128), all_lags=True)
    mtm._crosscorr(randn(128), randn(128), all_lags=False)


def test_mtm():
    mtm._other_dpss_method(64,4,10)


def test_Multitapering():
    p = MultiTapering(data_two_freqs(), 4,2)
    p()
    p = MultiTapering(marple_data, 4,2, method="eigen")
    p()
    print(p)
    p._str_title()
