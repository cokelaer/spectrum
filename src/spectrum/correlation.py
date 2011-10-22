"""
.. topic:: Correlation module


    Provides two correlation functions. :func:`CORRELATION` is slower than 
    :func:`xcorr`. However, the output is as expected by some other functions. 
    Ultimately, it should be replaced by :func:`xcorr`.
    
    For real data, the behaviour of the 2 functions is identical. However, for
    complex data, xcorr returns a 2-sides correlation.
 

    .. autosummary:: 

        ~spectrum.correlation.CORRELATION
        ~spectrum.correlation.xcorr
        
    .. codeauthor: Thomas Cokelaer, 2011



"""#from numpy.fft import fft, ifft
import numpy
from numpy import  arange, isrealobj
from pylab import rms_flat

__all__ = ['CORRELATION', 'xcorr']


def CORRELATION(x, y=None, maxlags=None, norm='unbiased'):
    r"""Correlation function

    This function should give the same results as :func:`xcorr` but it 
    returns the positive lags only. Moreover the algorithm does not use
    FFT as compared to other algorithms. 
     
    :param array x: first data array of length N
    :param array y: second data array of length N. If not specified, computes the 
        autocorrelation. 
    :param int maxlags: compute cross correlation between [0:maxlags]
        when maxlags is not specified, the range of lags is [0:maxlags].
    :param str norm: normalisation in ['biased', 'unbiased', None, 'coeff']
     
        * *biased*   correlation=raw/N, 
        * *unbiased* correlation=raw/(N-`|lag|`)
        * *coeff*    correlation=raw/(rms(x).rms(y))/N
        * None       correlation=raw

    :return: 
        * a numpy.array correlation sequence,  r[1,N]
        * a float for the zero-lag correlation,  r[0]
    
    The *unbiased* correlation has the form:
    
    .. math::

        \hat{r}_{xx} = \frac{1}{N-m}T \sum_{n=0}^{N-m-1} x[n+m]x^*[n] T 

    The *biased* correlation differs by the front factor only:

    .. math::

        \check{r}_{xx} = \frac{1}{N}T \sum_{n=0}^{N-m-1} x[n+m]x^*[n] T 

    with :math:`0\leq m\leq N-1`.
    
    .. doctest::
    
        >>> from spectrum import *
        >>> x = [1,2,3,4,5]
        >>> res = CORRELATION(x,x, maxlags=0, norm='biased')
        >>> res[0]
        11.0
        
    .. note:: this function should be replaced by :func:`xcorr`.
    
    .. seealso:: :func:`xcorr`
    """
    assert norm in ['unbiased','biased', 'coeff', None]
    #transform lag into list if it is an integer
    if y == None:
        y = x
    
    # N is the max of x and y
    N = max(len(x), len(y))
    if len(x)<N:
        y = y.copy()
        y.resize(N)
    if len(y)<N:
        y = y.copy()
        y.resize(N)
            
    #default lag is N-1
    if maxlags == None:
        maxlags = N - 1
    assert maxlags < N, 'lag must be less than len(x)'
    
    realdata = isrealobj(x) and isrealobj(y)
    #create an autocorrelation array with same length as lag
    if realdata == True:
        r = numpy.zeros(maxlags, dtype=float)
    else:
        r = numpy.zeros(maxlags, dtype=complex)

    if norm == 'coeff':
        rmsx = rms_flat(x)
        rmsy = rms_flat(y)
        
    for k in range(0, maxlags+1):
        nk = N - k - 1
        
        if realdata == True:
            sum = 0
            for j in range(0, nk+1):
                sum = sum + x[j+k] * y[j]
        else:
            sum = 0. + 0j
            for j in range(0, nk+1):
                sum = sum + x[j+k] * y[j].conjugate()
        if k == 0:
            if norm in ['biased', 'unbiased']:
                r0 = sum/float(N)
            elif norm == None:
                r0 = sum
            else:
                r0 =  1.
        else:
            if norm == 'unbiased':
                r[k-1] = sum / float(N-k)
            elif norm == 'biased':
                r[k-1] = sum / float(N)
            elif norm == None:
                r[k-1] = sum
            elif norm == 'coeff':
                r[k-1] =  sum/(rmsx*rmsy)/float(N)

    r = numpy.insert(r, 0, r0)
    return r
 

def xcorr(x, y=None, maxlags=None, norm='biased'):
    """Cross-correlation using numpy.correlate
    
    Estimates the cross-correlation (and autocorrelation) sequence of a random
    process of length N. By default, there is no normalisation and the output
    sequence of the cross-correlation has a length 2*N+1. 
    
    :param array x: first data array of length N
    :param array y: second data array of length N. If not specified, computes the 
        autocorrelation. 
    :param int maxlags: compute cross correlation between [-maxlags:maxlags]
        when maxlags is not specified, the range of lags is [-N+1:N-1].
    :param str option: normalisation in ['biased', 'unbiased', None, 'coeff']
     
    The true cross-correlation sequence is
    
    .. math:: r_{xy}[m] = E(x[n+m].y^*[n]) = E(x[n].y^*[n-m])

    However, in practice, only a finite segment of one realization of the 
    infinite-length random process is available.
    
    The correlation is estimated using numpy.correlate(x,y,'full'). 
    Normalisation is handled by this function using the following cases:

        * 'biased': Biased estimate of the cross-correlation function
        * 'unbiased': Unbiased estimate of the cross-correlation function
        * 'coeff': Normalizes the sequence so the autocorrelations at zero 
           lag is 1.0.

    :return:
        * a numpy.array containing the cross-correlation sequence (length 2*N-1)
        * lags vector
        
    .. note:: If x and y are not the same length, the shorter vector is 
        zero-padded to the length of the longer vector.
               
    .. rubric:: Examples
    
    .. doctest::
    
        >>> from spectrum import *
        >>> x = [1,2,3,4,5]
        >>> c, l = xcorr(x,x, maxlags=0, norm='biased')
        >>> c
        array([ 11.])
    
    .. seealso:: :func:`CORRELATION`.  
    """
    N = len(x)
    if y == None:
        y = x
    assert len(x) == len(y), 'x and y must have the same length. Add zeros if needed'
    assert maxlags <= N, 'maxlags must be less than data length'
    
    if maxlags == None:
        maxlags = N-1
        lags = arange(0, 2*N-1)
    else:
        assert maxlags < N
        lags = arange(N-maxlags-1, N+maxlags)
              
    res = numpy.correlate(x, y, mode='full')
    
    if norm == 'biased':
        Nf = float(N)
        res = res[lags] / float(N)    # do not use /= !! 
    elif norm == 'unbiased':
        res = res[lags] / (float(N)-abs(arange(-N+1, N)))[lags]
    elif norm == 'coeff':        
        Nf = float(N)
        rms = rms_flat(x) * rms_flat(y)
        res = res[lags] / rms / Nf
    else:
        res = res[lags]

    lags = arange(-maxlags, maxlags+1)        
    return res, lags
