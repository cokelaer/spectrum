import numpy
from spectrum import *
from spectrum import ma, arma_estimate
from spectrum.psd import *
from spectrum.errors import *
from numpy.testing import assert_array_almost_equal
import pylab
data = marple_data


def create_all_psd():
    f = pylab.linspace(0, 1, 4096)
    pylab.clf()
    #MA 15 order
    b, rho = ma(data, 15, 30)
    psd = arma2psd(B=b, rho=rho)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='MA 15')
    pylab.hold(True)
    
    #ARMA 15 order
    a, b, rho = arma_estimate(data, 15,15, 30)
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
    af, pf, ab, pb, pv = arcovar_marple(data, 15)
    psd = arma2psd(A=af, B=ab, rho=pf)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='covar 15')

    #modcovar method
    a, p, pv = modcovar_marple(data, 15)
    psd = arma2psd(A=a)
    newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='modcovar 15')

    #correlogram
    psd = CORRELOGRAMPSD(data, data, lag=15)
    newpsd = cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='correlogram 15')

    #minvar
    psd = minvar(data, 15)
    #newpsd = tools.cshift(psd, len(psd)/2) # switch positive and negative freq
    pylab.plot(f, 10 * pylab.log10(newpsd/max(newpsd)), label='MINVAR 15')

    #music
    psd,db = music(data, 15, 11)
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), '--',label='MUSIC 15')

    #ev music
    psd,db = ev(data, 15, 11)
    pylab.plot(f, 10 * pylab.log10(psd/max(psd)), '--',label='EV 15')


    pylab.legend(loc='upper left', prop={'size':10}, ncol=2)
    pylab.ylim([-80,10])
    pylab.savefig('psd_all.png')


create_all_psd()
