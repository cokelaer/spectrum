import numpy
from spectrum import *
from spectrum.psd import *
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
    print r


class test_spectrum():
    def __init__(self):
        #create a spectrum
        self.data = [1,2,3,4]
        self.init()
        assert self.s.N == 4

    def init(self):
        self.s = Spectrum(self.data)
        print self.s #__str__ when psd is not yet computed
        
    def test_attributes(self):
        self.init()
        # test sampling property
        self.s.sampling = 1024
        # test detrend property
        self.s.detrend = 'mean'
        try:
            self.s.detrend = self.s.detrend
            self.s.detrend = 'dummy'
            assert False
        except SpectrumChoiceError:
            assert True
        try:
            self.s.sides = self.s.sides
            self.s.sides = 'dummy'
            assert False
        except SpectrumChoiceError:
            assert True
        
    
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
    
    print s

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
