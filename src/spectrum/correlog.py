"""Correlogram PSD estimates

.. topic:: This module provides Correlograms methods 

    .. autosummary::

        CORRELOGRAMPSD
        Correlogram
        
    .. codeauthor:: Thomas Cokelaer 2011
 
    :References: See [Marple]_
"""


import numpy
from correlation import CORRELATION, xcorr
from window import Window
from numpy.fft import fft
from psd import FourierSpectrum


__all__ = ["CORRELOGRAMPSD", "Correlogram"]



def CORRELOGRAMPSD(X, Y=None, lag=-1, window='hamming', 
                    norm='unbiased', NFFT=4096, window_params={}, 
                    correlation_method='xcorr'):
    """PSD estimate using correlogram method.


    :param array X: complex or real data samples X(1) to X(N)
    :param array Y: complex data samples Y(1) to Y(N). If provided, computes 
        the cross PSD, otherwise the PSD is returned
    :param int lag:  highest lag index to compute. Must be less than N
    :param str window_name: see :mod:`window` for list of valid names
    :param str norm: one of the valid normalisation of :func:`xcorr` (biased, unbiased, coeff, None)
    :param int NPSD: the length of the final PSD 
    :param str correlation_method: either `xcorr` or `CORRELATION`. 
        CORRELATION should be removed in the future.

    :return:
        * Array of real (cross) power spectral density estimate values. This is
          a two sided array with negative values following the positive ones
          whatever is the input data (real or complex).
          
    .. rubric:: Description: 
    
    The exact power spectral density is the Fourier transform of the 
    autocorrelation sequence:
    
    .. math:: P_{xx}(f) = T \sum_{m=-\infty}^{\infty} r_{xx}[m] exp^{-j2\pi fmT}
    
    The correlogram method of PSD estimation substitutes a finite sequence of
    autocorrelation estimates :math:`\hat{r}_{xx}` in place of :math:`r_{xx}`.    
    This estimation can be computed with :func:`xcorr` or :func:`CORRELATION` by
    chosing a proprer lag `L`. The estimated PSD is then
    
    .. math:: \hat{P}_{xx}(f) = T \sum_{m=-L}^{L} \hat{r}_{xx}[m] exp^{-j2\pi fmT}
    
    The lag index must be less than the number of data samples `N`. Ideally, it
    should be around `L/10` [Marple]_ so as to avoid greater statistical
    variance associated with higher lags.
    
    To reduce the leakage of the implicit rectangular window and therefore to 
    reduce the bias in the estimate, a tapering window is normally used and lead
    to the so-called Blackman and Tukey correlogram:
    
    .. math:: \hat{P}_{BT}(f) = T \sum_{m=-L}^{L} w[m] \hat{r}_{xx}[m] exp^{-j2\pi fmT}

    The correlogram for the cross power spectral estimate is 
    
    .. math:: \hat{P}_{xx}(f) = T \sum_{m=-L}^{L} \hat{r}_{xx}[m] exp^{-j2\pi fmT}
    
    which is computed if :attr:`Y` is not provide. In such case, 
    :math:`r_{yx} = r_{xy}` so we compute the correlation only once.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import *

        psd = CORRELOGRAMPSD(marple_data, marple_data, lag=15)
        f = linspace(-0.5, 0.5, len(psd))
        psd = cshift(psd, len(psd)/2)
        plot(f, 10*log10(psd/max(psd)))
        axis([-0.5,0.5,-50,0])
        grid(True)

    .. seealso:: :func:`create_window`, :func:`CORRELATION`, :func:`xcorr`,
        :class:`Correlogram`.
    """
    N = len(X)
    assert lag<N, 'lag must be < size of input data'
    assert correlation_method in ['CORRELATION', 'xcorr']
    if Y == None:
        Y = numpy.array(X)
        crosscorrelation = False
    else:
        crosscorrelation = True
    
    if NFFT == None:
        NFFT = N
    psd = numpy.zeros(NFFT, dtype=complex)
    
    # Window should be centered around zero. Moreover, we want only the
    # positive values. So, we need to use 2*lag + 1 window and keep values on 
    # the right side.
    w = Window(2.*lag+1, window, **window_params)
    w = w.data[lag+1:]
   
    # compute the cross correlation
    if correlation_method == 'CORRELATION':
        rxy = CORRELATION (X, Y, maxlags=lag, norm=norm)
    elif correlation_method == 'xcorr':
        rxy, _l = xcorr (X, Y, maxlags=lag, norm=norm)
        rxy = rxy[lag:]
    
    # keep track of the first elt.
    psd[0] = rxy[0]
    
    # create the first part of the PSD
    psd[1:lag+1] = rxy[1:] * w
    
    
    # create the second part. 
    # First, we need to compute the auto or cross correlation ryx
    if crosscorrelation is True:
        # compute the cross correlation
        if correlation_method == 'CORRELATION':
            ryx = CORRELATION(Y, X, maxlags=lag, norm=norm)
        elif correlation_method == 'xcorr':
            ryx, _l = xcorr(Y, X, maxlags=lag, norm=norm)
            ryx = ryx[lag:]
        #print len(ryx), len(psd[-1:NPSD-lag-1:-1])
        
        psd[-1:NFFT-lag-1:-1] = ryx[1:].conjugate() * w
    else: #autocorrelation no additional correlation call required
        psd[-1:NFFT-lag-1:-1] = rxy[1:].conjugate() * w
    
    
    psd = numpy.real(fft(psd))

    return psd



class Correlogram(FourierSpectrum):
    """The Correlogram class provides an interface to :func:`CORRELOGRAMPSD`.
    
    It returns an object that inherits from :class:`FourierSpectrum` and 
    therefore ease the manipulation of PSDs.
       
    .. plot::
        :width: 80%
        :include-source:
        
        from spectrum import *
        p = Correlogram(data_cosine(N=1024), lag=15)
        p()
        p.plot()
        p.plot(sides='twosided')
    
    
    """
    def __init__(self, data, sampling=1., lag=-1,
                 window='hamming', NFFT=None, scale_by_freq=True, 
                 detrend=None):
        """**Correlogram Constructor**
        
        :param array data:
        :param float sampling:
        :param int lag:
        :param str window:
        :param int NFFT:
        :param bool scale_by_freq:
        :param str detrend:
        
        """
        super(Correlogram, self).__init__(data, 
                                          window=window, 
                                          sampling=sampling,
                                          NFFT=NFFT, 
                                          scale_by_freq=scale_by_freq,
                                          lag=lag,
                                          detrend=detrend)

    def __call__(self):
        psd = CORRELOGRAMPSD(self.data, self.data_y,
                             lag=self.lag,
                             window=self.window, 
                             NFFT=self.NFFT, 
                             #   scale_by_freq=self.scale_by_freq,
                             )
        if self.datatype == 'real':
            newpsd  = psd[self.NFFT/2:]*2
            newpsd[0] /= 2.
            newpsd = numpy.append(newpsd, psd[0])
            self.psd = newpsd
        else:
            self.psd = psd
        self.scale()
        #self.sides = 'twosided'
        #self.modified = False
        
    def _str_title(self):
        return "Correlogram PSD estimate\n"
    
    def __str__(self):
        return super(Correlogram, self).__str__()

