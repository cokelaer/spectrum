import numpy
import spectrum
from spectrum import tools
#from spectrum.psd import *
#from spectrum.errors import *
from numpy.testing import assert_array_almost_equal
import pylab
data = spectrum.marple_data


def create_all_psd():


    f = pylab.linspace(0, 1, 4096)
    pylab.clf()

    pylab.figure(figsize=(12,8))

    #MA 15 order
    b, rho = spectrum.ma(data, 15, 30)
    psd = spectrum.arma2psd(B=b, rho=rho)
    newpsd = tools.cshift(psd, len(psd)//2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='MA 15')

    #ARMA 15 order
    a, b, rho = spectrum.arma_estimate(data, 15,15, 30)
    psd = spectrum.arma2psd(A=a,B=b, rho=rho)
    newpsd = tools.cshift(psd, len(psd)//2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='ARMA 15,15')

    #yulewalker
    ar, P,c = spectrum.aryule(data, 15, norm='biased')
    psd = spectrum.arma2psd(A=ar, rho=P)
    newpsd = tools.cshift(psd, len(psd)//2) # switch positive and negative freq

    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='YuleWalker 15')

    #burg method
    ar, P,k = spectrum.arburg(data, order=15)
    psd = spectrum.arma2psd(A=ar, rho=P)
    newpsd = tools.cshift(psd, len(psd)//2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='Burg 15')

    #covar method
    af, pf, ab, pb, pv = spectrum.arcovar_marple(data, 15)
    psd = spectrum.arma2psd(A=af, B=ab, rho=pf)
    newpsd = tools.cshift(psd, len(psd)//2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='covar 15')

    #modcovar method
    a, p, pv = spectrum.modcovar_marple(data, 15)
    psd = spectrum.arma2psd(A=a)
    newpsd = tools.cshift(psd, len(psd)//2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='modcovar 15')

    #correlogram
    psd = spectrum.CORRELOGRAMPSD(data, data, lag=15)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='correlogram 15')

    #minvar
    psd = spectrum.minvar(data, 15)
    #newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='MINVAR 15')

    #music
    psd,db = spectrum.music(data, 15, 11)
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), '--',label='MUSIC 15')

    #ev music
    psd,db = spectrum.ev(data, 15, 11)
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), '--',label='EV 15')


    pylab.legend(loc='upper left', prop={'size':10}, ncol=2)
    pylab.ylim([-80,10])
    pylab.savefig('psd_all.png')


create_all_psd()
