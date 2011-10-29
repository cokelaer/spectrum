"""ARMA and MA estimates, ARMA and MA PSD estimates.

.. topic:: ARMA model and Power Spectral Densities. 

    .. autosummary::
       :nosignatures:

       arma_estimate
       ma
       arma2psd
       parma
       pma

    .. codeauthor:: Thomas Cokelaer 2011

    :References: See [Marple]_
"""
from numpy import zeros, append, insert, real
from numpy.fft import fft
from spectrum.correlation import CORRELATION
from spectrum.covar import arcovar, arcovar_marple
import spectrum.yulewalker as yulewalker
from spectrum.psd import ParametricSpectrum


__all__ = ["arma2psd", "arma_estimate", "ma", "pma", "parma"]


def arma2psd(A=None, B=None, rho=1., T=1., NPSD=4096, sides='default', norm=False):
    r"""Computes power spectral density given ARMA values.

    This function computes the power spectral density values
    given the ARMA parameters of an ARMA model. It is suppose that
    the driving sequence is a white noise process of zero mean and 
    variance :math:`\rho_w`. The sampling frequency and noise variance are
    used to scale the PSD output, which length is set by the user with the 
    `NPSD` parameter. 
    
    :param array A:   Array of AR parameters (complex or real)
    :param array B:   Array of MA parameters (complex or real)
    :param float rho: White noise variance to scale the returned PSD
    :param float T:   Sample interval in seconds to scale the returned PSD
    :param int NPSD:  Final size of the PSD
    :param str sides: Default PSD is two-sided, but sides can be set to centerdc.
    
    .. warning:: By convention, the AR or MA arrays does not contain the
        A0=1 value.
    
    If :attr:`B` is None, the model is a pure AR model. If :attr:`A` is None, 
    the model is a pure MA model.

    :return: two-sided PSD 

    .. rubric:: Details:
    
    AR case: the power spectral density is:
    
    .. math:: P_{ARMA}(f) = T \rho_w \left|\frac{B(f)}{A(f)}\right|^2
    
    where:
    
    .. math:: A(f) = 1 + \sum_{k=1}^q b(k) e^{-j2\pi fkT}
    .. math:: B(f) = 1 + \sum_{k=1}^p a(k) e^{-j2\pi fkT}
            
    .. rubric:: **Example:**

    .. plot::
        :width: 80%
        :include-source:

        import spectrum.arma
        from pylab import plot, log10, hold, legend
        plot(10*log10(spectrum.arma.arma2psd([1,0.5],[0.5,0.5])), label='ARMA(2,2)')
        plot(10*log10(spectrum.arma.arma2psd([1,0.5],None)), label='AR(2)')
        plot(10*log10(spectrum.arma.arma2psd(None,[0.5,0.5])), label='MA(2)')
        legend()
        
    :References: [Marple]_
    """
    if NPSD == None:
        NPSD = 4096
    
    if A == None and B == None:
        raise ValueError("Either AR or MA model must be provided")
    
    psd = zeros(NPSD, dtype=complex)
    
    if A != None:
        ip = len(A)
        den = zeros(NPSD, dtype=complex)
        den[0] = 1.+0j
        for k in range(0, ip):
            den[k+1] = A[k]
        denf = fft(den, NPSD)

    if B != None:
        iq = len(B)
        num = zeros(NPSD, dtype=complex)
        num[0] = 1.+0j
        for k in range(0, iq):
            num[k+1] = B[k]
        numf = fft(num, NPSD)

    if A != None and B != None:
        psd = rho * T * abs(numf)**2. /  abs(denf)**2.
    elif A != None:
        psd = rho * T / abs(denf)**2.
    elif B != None:
        psd = rho * T * abs(numf)**2.
        
    psd = real(psd)
    # The PSD is a twosided PSD.
    # to obtain the centerdc
    if sides != 'default':
        import tools
        assert sides in ['centerdc']
        if sides == 'centerdc':
            psd = tools.twosided_2_centerdc(psd)
        
    if norm == True:
        psd /= max(psd)
        
    return psd


def arma_estimate(X, P, Q, lag):
    """Autoregressive and moving average estimators.

    This function provides an estimate of the autoregressive 
    parameters, the moving average parameters, and the driving
    white noise variance of  an ARMA(P,Q) for a complex or real data sequence.
    
    The parameters are estimated using three steps:    
    
        * Estimate the AR parameters from the original data based on a least 
          squares modified Yule-Walker technique,  
        * Produce a residual time sequence by filtering the original data 
          with a filter based on the AR parameters, 
        * Estimate the MA parameters from the residual time sequence.

    :param array X: Array of data samples (length N)
    :param int P: Desired number of AR parameters
    :param int Q: Desired number of MA parameters
    :param int lag: Maximum lag to use for autocorrelation estimates

    :return:
        * A     - Array of complex P AR parameter estimates
        * B     - Array of complex Q MA parameter estimates
        * RHO   - White noise variance estimate

    .. note::
      *  lag must be >= Q (MA order)

    **dependencies**: 
        * :meth:`spectrum.correlation.CORRELATION`
        * :meth:`spectrum.covar.arcovar`
        * :meth:`spectrum.arma.ma`

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import *
        a,b, rho = arma_estimate(marple_data, 15, 15, 30)
        psd = arma2psd(A=a, B=b, rho=rho, sides='centerdc', norm=True)
        plot(10 * log10(psd))
        ylim([-50,0])

    :reference: [Marple]_
    """
    R = CORRELATION(X, maxlags=lag, norm='unbiased')
    R0 = R[0]
    #C   Estimate the AR parameters (no error weighting is used).
    #C   Number of equation errors is M-Q .
    MPQ = lag - Q + P
    
    N = len(X)
    Y = zeros(N-P, dtype=complex)

    for K in range(0, MPQ):
        KPQ = K + Q - P+1
        if KPQ < 0:
            Y[K] = R[-KPQ].conjugate()
        if KPQ == 0:
            Y[K] = R0
        if KPQ > 0:
            Y[K] = R[KPQ]

    # The resize is very important for the normalissation.
    Y.resize(lag)
    if P <= 4:
        res = arcovar_marple(Y.copy(), P)    #! Eq. (10.12)
        ar_params = res[0]
    else:
        res = arcovar(Y.copy(), P)    #! Eq. (10.12)
        ar_params = res[0]
        
    # the .copy is used to prevent a reference somewhere. this is a bug
    # to be tracked down.
    Y.resize(N-P)
    
    #C   Filter the original time series
    for k in range(P, N):
        SUM = X[k]
        #SUM += sum([ar_params[j]*X[k-j-1] for j in range(0,P)])
        for j in range(0, P):
            SUM = SUM + ar_params[j] * X[k-j-1]   #! Eq. (10.17)
        Y[k-P] = SUM

    #  Estimate the MA parameters (a "long" AR of order at least 2*IQ
    #C   is suggested)
    #Y.resize(N-P)
    ma_params, rho = ma(Y, Q, 2*Q)     #! Eq. (10.3)

    return ar_params, ma_params, rho


class parma(ParametricSpectrum):
    """Class to create PSD using ARMA estimator.
    
    See :func:`arma_estimate` for description.
    
    .. plot::
        :width: 80%
        :include-source:
    
        from spectrum import *
        p = parma(marple_data, 4, 4, 30, NFFT=4096)
        p()
        p.plot(sides='centerdc')
    
    
    """
    def __init__(self, data, P, Q, lag, NFFT=None, sampling=1., scale_by_freq=False):
        """**Constructor:**
        
        For a detailed description of the parameters, see :func:`arma_estimate`.
        
        :param array data:
        :param int P:
        :param int Q:
        :param int lag:
        :param int NFFT:
        :param float sampling:
                
        """
        super(parma, self).__init__(data, ma_order=Q, ar_order=P, lag=lag, 
                                    NFFT=NFFT, sampling=sampling, 
                                    scale_by_freq=scale_by_freq)
        self.lag = lag

    def __call__(self):
        ar_params, ma_params, rho = arma_estimate(self.data, self.ar_order, 
                                         self.ma_order, self.lag)
        self.ma = ma_params
        self.ar = ar_params
        self.rho = rho
        psd = arma2psd(A=self.ar, B=self.ma, rho=self.rho, 
                      T=self.sampling, NPSD=self.NFFT)
        #self.psd = psd
        if self.datatype == 'real':
            newpsd  = psd[0:self.NFFT/2]*2
            newpsd[0] /= 2.
            newpsd = append(newpsd, psd[-1])
            self.psd = newpsd
        else:
            self.psd = psd
        if self.scale_by_freq is True:
            self.scale()

    def _str_title(self):
        return "Periodogram PSD estimate\n"

    def __str__(self):
        return super(parma, self).__str__()



class pma(ParametricSpectrum):
    """Class to create PSD using MA estimator.

    See :func:`ma` for description.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        p = pma(marple_data, 15, 30, NFFT=4096)
        p()
        p.plot(sides='centerdc')
    
    
    """
    def __init__(self, data, Q, M, NFFT=None, sampling=1., scale_by_freq=False):
        """**Constructor:**
        
        For a detailed description of the parameters, see :func:`ma`.
        
        :param array data:
        :param int Q: MA order
        :param int M: AR model used to estimate the MA parameters
        :param int NFFT: 
        :param float sampling:
                
        """
        super(pma, self).__init__(data, ma_order=Q, ar_order=M, 
                                  NFFT=NFFT, sampling=sampling,
                                  scale_by_freq=scale_by_freq)

    def __call__(self):
        ma_params, rho = ma(self.data, self.ma_order, self.ar_order)
        self.ma = ma_params
        self.rho = rho
        psd = arma2psd(A=None, B=self.ma, rho=self.rho, 
                      T=self.sampling, NPSD=self.NFFT)
        #self.psd = psd
        if self.datatype == 'real':
            import tools
            self.psd = tools.twosided_2_onesided(psd)
            #newpsd  = psd[self.NFFT/2:]*2
            #newpsd[0] /= 2.
            #newpsd = numpy.append(newpsd, psd[0])
            #self.psd = newpsd
        else:
            self.psd = psd
        if self.scale_by_freq is True:
            self.scale()
        self.modified = False
        
    def _str_title(self):
        return "Periodogram PSD estimate\n"
    
    def __str__(self):
        return super(pma, self).__str__()
    
        


def ma(X, Q, M):
    """Moving average estimator.

    This program provides an estimate of the moving average parameters
    and driving noise variance for a data sequence based on a
    long AR model and a least squares fit.

    :param array X: The input data array
    :param int Q: Desired MA model order (must be >0 and <M)
    :param int M: Order of "long" AR model (suggest at least 2*Q )

    :return:
        * MA    - Array of Q complex MA parameter estimates 
        * RHO   - Real scalar of white noise variance estimate

    .. plot::
        :width: 80%
        :include-source:

        from pylab import *
        from spectrum import *
        # Estimate 15 Ma parameters
        b, rho = ma(marple_data, 15, 30)
        # Create the PSD from those MA parameters
        psd = arma2psd(B=b, rho=rho, sides='centerdc')
        # and finally plot the PSD
        plot(linspace(-0.5, 0.5, 4096), 10 * log10(psd/max(psd)))
        axis([-0.5, 0.5, -30, 0])

    :reference: [Marple]_
    """
    if Q <= 0 or Q >= M:
        raise ValueError('Q(MA) must be in ]0,lag[')

    #C   Fit a high-order AR to the data
    a, rho, _c = yulewalker.aryule(X, M, 'biased')   #! Eq. (10.5)

    #add an element unity to the AR parameter array
    a = insert(a, 0, 1)

    #C   Find MA parameters from autocorrelations by Yule-Walker method
    ma_params, _p, _c = yulewalker.aryule(a, Q, 'biased')    #! Eq. (10.7)

    return ma_params, rho
