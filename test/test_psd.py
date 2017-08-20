import numpy
from spectrum import *
from spectrum.psd import *
from spectrum.psd import Range # not exposed in __all__, so need to import manually
from spectrum.errors import *
from numpy.testing import assert_array_almost_equal
import pylab
data = marple_data


def test_psd_module_range():
    N = 10
    fs = 1024.
    r = Range(N, fs)
    assert r.N == N
    r.df == float(fs)/N
    r.sampling == fs

    #test attributes
    r.N = 1024;
    assert r.df == 1.
    r.sampling = 2048;
    assert r.df == 2.

    # test methods
    r.N = 10;
    r.onesided()
    r.twosided()
    r.N = 11;
    r.onesided()
    r.twosided()
    print(r)



def test_spectrum_simple():
    data = [1,2,3,4]
    s = Spectrum(data)
    s.N == 4
    print(s)
    assert s.range.df == 0.25
    s.method

    # test sampling property and detrend
    s.sampling = 1024
    # test detrend property
    s.detrend = 'mean'
    try:
        s.detrend = s.detrend
        s.detrend = 'dummy'
        assert False
    except SpectrumChoiceError:
        assert True
    try:
        s.sides = s.sides
        s.sides = 'dummy'
        assert False
    except SpectrumChoiceError:
        assert True

    # test nextpow2 and float values for NFFT attribute
    data = [1,2,3,4,5]
    s = Spectrum(data)
    s.NFFT = "nextpow2"
    assert s.NFFT == 8
    try:
        s.NFFT = -1.1
        assert False
    except:
        assert True

    # change the data
    assert s.modified is True
    s.modified = False
    s.data = [1,2,3,4]
    assert s.modified is True
   
    # test the data_y attribute
    assert s.data_y is None
    s.data_y = [1,2,3,4]
    assert s.data_y == [1,2,3,4]
    
    try:
        s.frequencies("dummy") 
        assert False
    except:
        assert True

    # test the one/two sided conversion
    p = parma(data_two_freqs(), 8,4,15)
    p()
    aa = p.psd.copy()
    p.sides = "twosided"
    p.sides = "onesided"
    p.sides = "centerdc"
    p.sides = "twosided"
    p.sides = "centerdc"
    p.sides = "onesided"
    assert all(p.psd == aa)

    try:
        p.sides = "dummy"
        assert False
    except:
        assert True

    p = parma(data_two_freqs(), 8,4,15)
    try:
        p.plot()
        assert False
    except:
        assert True

    # power method
    p = parma(marple_data, 8,4,10)
    p()
    assert p.power() < 7381 and p.power() > 7380
    p.scale_by_freq = True
    assert p.power() > 0.286 and p.power() < 0.287



class test_spectrum():
    def __init__(self):
        #create a spectrum
        self.data = [1,2,3,4]
        self.init()
        assert self.s.N == 4

    def init(self):
        self.s = Spectrum(self.data)
        print(self.s) #__str__ when psd is not yet computed

    def test_sides(self):
        """Functional tests on p.sides = ..."""
        # test all sides cases by starting from a side, converting to anothre
        # and convert back to the original
        self.init()
        assert self.s.NFFT == 4
        assert self.s.sides == 'onesided'
        self.s.psd = [10, 2, 3, 4, 6]
        for x in ['onesided', 'twosided', 'centerdc']:
            for y in ['onesided', 'twosided', 'centerdc']:
                if x != y:
                    self.s.sides = x
                    p0 = self.s.psd.copy()
                    self.s.sides = y
                    self.s.sides = x
                    assert_array_almost_equal(self.s.psd, p0)

    def test_frequencies(self):
        """Functional tests one s.frequencies(...)"""
        self.init()
        self.s.frequencies('onesided')
        self.s.frequencies('twosided')
        self.s.frequencies('centerdc')
        try:
            self.s.frequencies("dummy")
            assert False
        except:
            assert True

    def test_set_psd(self):
        self.init()
        self.s.psd = [1,2,3]
        assert self.s.NFFT == 4
        self.s.psd = numpy.array([1,2,3])
        assert self.s.NFFT == 4

def test_fourier_spectrum():

    from spectrum import datasets
    from spectrum import FourierSpectrum
    s = FourierSpectrum(datasets.data_cosine(), sampling=1024, NFFT=512)
    s.periodogram()
    s.plot()


    s = FourierSpectrum(datasets.marple_data, sampling=1024, NFFT=512)
    s.periodogram()
    s.plot()

    print(s)

def test_arma_spectrum():
    from spectrum import datasets
    from spectrum import FourierSpectrum
    s = ParametricSpectrum(datasets.data_cosine(), sampling=1024, ar_order=30, ma_order=15)
    #s.ma()
    #s.arma()


def _test_create_all_psd():
    f = pylab.linspace(0, 1, 4096)
    pylab.clf()
    #MA 15 order
    b, rho = MA(data, 15, 30)
    psd = arma2psd(B=b, rho=rho)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='MA 15')
    pylab.hold(True)

    #ARMA 15 order
    a, b, rho = ARMA(data, 15,15, 30)
    psd = arma2psd(A=a,B=b, rho=rho)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='ARMA 15,15')
    pylab.hold(True)

    #yulewalker
    ar, P,c = aryule(data, 15, norm='biased')
    psd = arma2psd(A=ar, rho=P)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq

    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='YuleWalker 15')

    #burg method
    ar, P,k = arburg(data, order=15)
    psd = arma2psd(A=ar, rho=P)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq

    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='Burg 15')

    #covar method
    af, pf, ab, pb, pv = arcovar(data, 15)
    psd = arma2psd(A=af, B=ab, rho=pf)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='covar 15')

    #modcovar method
    a, p, pv = MODCOVAR(data, 15)
    psd = arma2psd(A=a)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='modcovar 15')

    #correlogram
    psd = CORRELOGRAMPSD(data, data, lag=15)
    newpsd = cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='correlogram 15')

    #minvar
    psd = MINVAR(data, 15, 1.)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='MINVAR 15')

    #music
    psd,db = EIGENFRE(data, 15, 11, method='music')
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), '--',label='MUSIC 15')

    #ev music
    psd,db = EIGENFRE(data, 15, 11, method='ev')
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), '--',label='EV 15')


    pylab.legend(loc='upper left', prop={'size':10}, ncol=2)
    pylab.ylim([-80,10])
    pylab.savefig('psd_all.png')


#test_create_all_psd()
