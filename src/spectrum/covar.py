import numpy
from psd import ParametricSpectrum
debug = False


__all__ = ["arcovar", "arcovar_marple", "pcovar", 'arcovar_marple']

def arcovar_marple(x, order):
    r"""Estimate AR model parameters using covariance method
    
    This implementation is based on [Marple]_. This code is far more 
    complicated and slower than :func:`arcovar` function, which is now the official version.
    See :func:`arcovar` for a detailed description of Covariance method.
    
    This function should be used in place of arcovar only if order<=4, for 
    which :func:`arcovar` does not work. 

    Fast algorithm for the solution of the covariance least squares normal 
    equations from Marple.

    :param array X:  Array of complex data samples
    :param int oder: Order of linear prediction model

    :return:
        * AF   - Array of complex forward linear prediction coefficients
        * PF   - Real forward linear prediction variance at order IP
        * AB   - Array of complex backward linear prediction coefficients
        * PB   - Real backward linear prediction variance at order IP
        * PV   - store linear prediction coefficients

    .. note:: this code and the original code in Marple diverge for ip>10.
        it seems that this is related to single precision used with 
        complex type in fortran whereas numpy uses double precision for 
        complex type.

    :validation: the AR parameters are the same as those returned by 
        a completely different function :func:`arcovar`.
        
    :References: [Marple]_
    """
    assert len(x) >= order, "X must be dimensioned >=N"
    
    #   ----------------------------------------------------- Initialization
    x = numpy.array(x)
    N = len(x)
    
    
    # Equations 8.C.42
    r0 = sum(abs(x)**2.)
    r1 = abs(x[0])**2
    rN = abs(x[N-1])**2
    
    pf = r0 - r1
    pb = r0 - rN
    
    delta = 1. - r1 / r0
    gamma = 1. - rN / r0 

    c = numpy.zeros(N, dtype=complex)
    d = numpy.zeros(N, dtype=complex)
    r = numpy.zeros(N, dtype=complex)
    af = numpy.zeros(N, dtype=complex)
    ab = numpy.zeros(N, dtype=complex)
    
    c[0] = x[N-1].conjugate() / r0
    d[0] = x[0].conjugate() / r0
    
    # special case
    if order == 0:
        pf = r0 / float(N)
        pb = pf
        return af, pf, ab, pb, 0
    
    # ----------------------------------------------------------  MAIN LOOP

    #ip +1 because we want to enter in the loop to run the first part of the code. 
    pbv = []
    for m in range(0, order+1):
        if debug:
            print '----------------------------m=', m
            print c[0:2]
            print d[0:2]
        r1 = 1./pf
        r2 = 1./pb
        r3 = 1./delta
        r4 = 1./gamma
        if debug:
            print 'starting r1r2r3r4=', r1, r2, r3, r4, pf, pb, delta, gamma
        #Order update: AF and AB vectors ; time update: C and D vectors
        temp = 0.+0.j
        for k in range(m+1, N):
            temp = temp + x[k]*x[k-m-1].conjugate()

        if debug:
            print 'temp=', temp
        r[m] =  temp.conjugate()
        theta = x[0] * c[m]
        if debug:
            print 'theta', theta
            print 'cccccccccc', c[0:2]
            print 'dddddd', d[0:2]
        if m == 0:
            pass
        else:
            for k in range(0, m):
                theta = theta + x[m-k] * c[k]                   # Eq. (8.C.39)
                r[k] = r[k] - x[N-m-1] * x[N-m+k].conjugate() # Eq. (8.C.32)
                temp = temp + af[m-k-1] * r[k].conjugate()
                #print 'loop1 k=', k
                #print '         theta=',theta, 'r[k]=',r[k], 'temp=', temp
                #print '         c=',c[k], 'af=',af[m-k-1]
        if m>0:
            if debug:
                print m, N-m
                print 'Xk=0',x[m-0],x[N-m-1], x[N-m+0]
            if m>1:
                if debug:
                    print 'Xk=1',x[m-1],x[N-m-1], x[N-m+1]
        c1 = -temp * r2
        c2 = -r1 * temp.conjugate()
        c3 = theta * r3
        c4 = r4 *theta.conjugate()
        if debug:
            print 'c1c2c3c4 before af=',c1 ,c2 ,c3 ,c4
        af[m] = c1                    #             ! Eq. (8.C.19)
        ab[m] = c2                    #             ! Eq. (8.C.22)
        save = c[m]
        c[m] = save + c3*d[m]
        d[m] = d[m] + c4*save
        if debug:
            print 'res',m,'af[m]=',af[m], ab[m], save, 'temp=',temp

        if m == 0:
            pass
        else:
            if debug:print 'af before', af[0:2]
            for k in range(0, m):
                save = af[k]
                af[k] = save + c1 * ab[m-k-1] # Eq. (8.C.18)
                ab[m-k-1] = ab[m-k-1] + c2 * save   # Eq. (8.C.21)
                
                save = c[k]
                c[k] = save + c3*d[k]       # Eq. (8.C.37)
                d[k] = d[k] + c4*save       # Eq. (8.C.38)
                if debug:
                    print 'loop2 k=', k
                    print '      af[k]=', af[k]
                    print '      ab[m-k-1]=', ab[m-k-1]
                    print '      c[k]=', c[k]
                    print '      d[k]=', d[k]
            if debug:
                if debug:
                    print 'af after=', af[0:2]
                    print 'ab=', ab[0:2]

        r5 = temp.real**2 + temp.imag**2
        pf = pf - r5*r2         # Eq. (8.C.20)
        pb = pb - r5*r1         # Eq. (8.C.23)
        r5 = theta.real**2 + theta.imag**2
        delta = delta - r5*r4               # Eq. (8.C.39)
        gamma = gamma - r5*r3               # Eq. (8.C.40)
        if debug:
            print 'r5r2r1deltagamma', r5, r2, r1 , delta, gamma
            print 'pf before norm', pf, pb, N-m-1
        if m != order-1:
            pass
        else:
            pf = pf / float(N-m-1)
            pb = pb / float(N-m-1)
            if debug:
                print 'ENDING', N-m-1
            break
        if debug:
            print 'pf and pb', pf, pb
        if pf > 0 and pb > 0:
            pass
        else: 
            ValueError("Negative PF or PB value")
        if (delta > 0. and delta <=1 and gamma > 0. and gamma <=1):
            pass
        else:
            ValueError("Invalid delta or gamma value")

        #C   Time update:  AF and AB vectors; order update:  C and D vectors

        r1 = 1./pf
        r2 = 1./pb
        r3 = 1./delta
        r4 = 1./gamma
        if debug:
            print '--------time update', r1, r2, r3, r4, m+1, N-m-1, x[m+1], x[N-m-2]
        ef = x[m+1]
        eb = x[(N-1)-m-1]
        if debug:
            print 'delta, gamma=', delta, gamma
        if debug:
            print 'before eb=', eb, ' ef=', ef
            
            
        for k in range(0,m+1):
            #print 'k=', k, 'ef=', ef, ' eb=',eb,' af=',af[k], ab[k]
            #print x[m-k],x[N-m+k-1]
            ef = ef + af[k] * x[m-k]             # Eq. (8.C.1)
            eb = eb + ab[k] * x[N-m+k-1]                   # Eq. (8.C.2)
           
        #ef = sum(af)
            
        if debug:
            print 'efweb', ef , eb
        c1 = ef*r3
        c2 = eb*r4
        c3 = eb.conjugate() * r2
        c4 = ef.conjugate() * r1
        if debug:
            print 'c1c2c3c4', c1, c2, c3, c4
            print 'af before', af[0:2]
        for k in range(m, -1, -1):
            save = af[k]
            af[k] = save + c1 * d[k]                    #  Eq. (8.C.33)
            d[k+1] = d[k] + c4 * save                    # Eq. (8.C.25)
            save = ab[k]
            ab[k] = save + c2 * c[m-k]                 # Eq. (8.C.35)
            c[m-k] = c[m-k] + c3 * save              # Eq. (8.C.24)
        if debug:
            print 'af after', af[0:2]
            print 'd', d[0:2]
            print 'ab', ab[0:2]
            print 'c', c[0:2]
        if debug:print 'Pb before', pf, pb
        c[m+1] = c3
        d[0] = c4
        #r5 = abs(ef)**2
        r5 = ef.real**2 + ef.imag**2
        pf = pf - r5 * r3                              # Eq. (8.C.34)
        
        delta = delta-r5 * r1                        # Eq. (8.C.30)
        #r5 = abs(eb)**2
        r5 = eb.real**2 + eb.imag**2
        pb = pb - r5 * r4                              # Eq. (8.C.36)
        if debug:
            print 'Pb---------------------', m, pb, r5, r4
        gamma = gamma-r5*r2                        # Eq. (8.C.31)
        if debug:
            print 'Gamma', gamma
            print 'r5 r3,r1,and delta', r5, r3, r1, delta
            print 'eb=', eb
            print 'ef=', ef
            print 'pf=', pf
            print 'pb=', pb
        pbv.append(pb)

        if (pf > 0. and pb > 0.):
            pass
        else:
            ValueError("Negative PF or PB value")
        if debug:
            print delta, gamma
        if (delta > 0. and delta <= 1.) and (gamma > 0. and gamma <= 1.):
            pass
        else:
            ValueError("Invalid delta or gamma value")


    #af=array of forward coeff   
    #ab=array of barward coeff   
    #pb=backward variance
    #pf=forward variance
    return af, pf, ab, pb, pbv


def arcovar(x, order):
    r"""Simple and fast implementation of the covariance AR estimate
    
    This code is 10 times faster than :func:`arcovar_marple` and more importantly
    only 10 lines of code, compared to a 200 loc for :func:`arcovar_marple`


    :param array X:  Array of complex data samples
    :param int oder: Order of linear prediction model

    :return:
        * a - Array of complex forward linear prediction coefficients
        * e - error

    The covariance method fits a Pth order autoregressive (AR) model to the 
    input signal, which is assumed to be the output of 
    an AR system driven by white noise. This method minimizes the forward 
    prediction error in the least-squares sense. The output vector 
    contains the normalized estimate of the AR system parameters
    
    The white noise input variance estimate is also returned.
    
    If is the power spectral density of y(n), then:
    
    .. math:: \frac{e}{\left| A(e^{jw}) \right|^2} = \frac{e}{\left| 1+\sum_{k-1}^P a(k)e^{-jwk}\right|^2}
    
    Because the method characterizes the input data using an all-pole model, 
    the correct choice of the model order p is important.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import plot, log10, linspace, axis
        ar_values, error = arcovar(marple_data, 15)
        psd = arma2psd(ar_values, sides='centerdc')
        plot(linspace(-0.5, 0.5, len(psd)), 10*log10(psd/max(psd)))
        axis([-0.5, 0.5, -60, 0])

    .. seealso:: :class:`pcovar`
    
    :validation: the AR parameters are the same as those returned by 
        a completely different function :func:`arcovar_marple`.
        
    :References: [Mathworks]_ 
    """
    
    from spectrum import corrmtx
    import scipy.linalg
    
    X = corrmtx(x, order, 'covariance')
    Xc = numpy.matrix(X[:, 1:]) 
    X1 = numpy.array(X[:, 0])
    
    # Coefficients estimated via the covariance method
    # Here we use lstsq rathre than solve function because Xc is not square 
    # matrix

    a, _residues, _rank, _singular_values = scipy.linalg.lstsq(-Xc, X1)
    
    # Estimate the input white noise variance
    Cz = numpy.dot(X1.conj().transpose(), Xc)
    e = numpy.dot(X1.conj().transpose(), X1) + numpy.dot(Cz, a)
    assert e.imag < 1e-4, 'wierd behaviour'
    e = float(e.real) # ignore imag part that should be small
    
    return a, e


class pcovar(ParametricSpectrum):
    """Class to create PSD based on covariance algorithm 
    
    See :func:`arcovar` for description.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        p = pcovar(marple_data, 15, NFFT=4096)
        p()
        p.plot(sides='centerdc')

    .. seealso:: :class:`arcovar`
    """
    def __init__(self, data, order, NFFT=None, sampling=1., 
                 scale_by_freq=False):
        """**Constructor**

        For a detailled description of the parameters, see :func:`arcovar`.

        :param array data:     input data (list or numpy.array)
        :param int order:  
        :param int NFFT:       total length of the final data sets (padded with zero if needed; default is 4096)
        :param float sampling: sampling frequency of the input :attr:`data`.


        """
        super(pcovar, self).__init__(data, ar_order=order, 
                                            NFFT=NFFT, sampling=sampling, 
                                            scale_by_freq=scale_by_freq)

    def __call__(self):
        from spectrum import arma2psd
        ar, _e = arcovar(self.data, self.ar_order)
        self.ar = ar
        psd = arma2psd(A=ar, T=self.sampling, NFFT=self.NFFT)
        
        if self.datatype == 'real':
            from tools import twosided_2_onesided
            newpsd  = twosided_2_onesided(psd)
            #\psd[0:self.NFFT/2]*2
            #newpsd[0] /= 2.
            #newpsd = numpy.append(newpsd, psd[-1])
            self.psd = newpsd
        else:
            self.psd = psd
        if self.scale_by_freq is True:
            self.scale()
        
    def _str_title(self):
        return "Covariance PSD estimate\n"
    
    def __str__(self):
        return super(pcovar, self).__str__()


