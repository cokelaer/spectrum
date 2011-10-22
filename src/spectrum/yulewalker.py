"""Yule Walker method to estimate AR values.

.. topic:: Estimation of AR values using Yule-Walker method 

    .. autosummary::

        aryule
        pyule

    .. codeauthor:: Thomas Cokelaer 2011

"""

from correlation import CORRELATION
from levinson import LEVINSON
from psd import ParametricSpectrum

__all__ = ['aryule', 'pyule']

def aryule(X, order, norm='biased', allow_singularity=True):
    r"""Compute AR coefficients using Yule-Walker method

    :param X: Array of complex data values, X(1) to X(N)
    :param int order: Order of autoregressive process to be fitted (integer)
    :param str norm: Use a biased or unbiased correlation.
    :param bool allow_singularity: 

    :return:
        * AR coefficients (complex)
        * variance of white noise (Real)
        * reflection coefficients for use in lattice filter 
 
    .. rubric:: Description:
    
    The Yule-Walker method returns the polynomial A corresponding to the 
    AR parametric signal model estimate of vector X using the Yule-Walker
    (autocorrelation) method. The autocorrelation may be computed using a 
    **biased** or **unbiased** estimation. In practice, the biased estimate of
    the autocorrelation is used for the unknown true autocorrelation. Indeed, an 
    unbiased estimate may result in nonpositive-definite autocorrelation matrix.
    So, a biased estimate leads to a stable AR filter. 
    The following matrix form represents the Yule-Walker equations. The are 
    solved by means of the Levinson-Durbin recursion:
    
     .. math::
    
        \left( \begin{array}{cccc}
        r(1) & r(2)^* & \dots & r(n)^*\\
        r(2) & r(1)^* & \dots & r(n-1)^*\\
        \dots & \dots & \dots & \dots\\
        r(n) & \dots & r(2) & r(1) \end{array} \right)    
        \left( \begin{array}{cccc}
        a(2)\\
        a(3) \\
        \dots \\
        a(n+1)  \end{array} \right)
        =
        \left( \begin{array}{cccc}
        -r(2)\\
        -r(3) \\
        \dots \\
        -r(n+1)  \end{array} \right)

    The outputs consists of the AR coefficients, the estimated variance of the
    white noise process, and the reflection coefficients. These outputs can be
    used to estimate the optimal order by using :mod:`~spectrum.criteria`.
         
    .. rubric:: Examples:

    From a known AR process or order 4, we estimate those AR parameters using
    the aryule function.
    
    .. doctest::
    
        >>> from scipy.signal import lfilter
        >>> from spectrum import *
        >>> from numpy.random import randn
        >>> A  =[1, -2.7607, 3.8106, -2.6535, 0.9238]
        >>> noise = randn(1, 1024)
        >>> y = lfilter([1], A, noise); 
        >>> #filter a white noise input to create AR(4) process
        >>> [ar, var, reflec] = aryule(y[0], 4)
        >>> # ar should contains values similar to A

    The PSD estimate of a data samples is computed and plotted as follows:
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import *

        ar, P, k = aryule(marple_data, 15, norm='biased')
        psd = arma2psd(ar)
        plot(linspace(-0.5, 0.5, 4096), 10 * log10(psd/max(psd)))
        axis([-0.5, 0.5, -60, 0])

    .. note:: The outputs have been double checked against (1) octave outputs 
        (octave has norm='biased' by default) and (2) Marple test code.
    
    .. seealso:: This function uses :func:`~spectrum.levinson.LEVINSON` and 
        :func:`~spectrum.correlation.CORRELATION`. See the :mod:`~spectrum.criteria`
        module for criteria to automatically select the AR order.
        
    :References: [Marple]_
    
    """
    assert norm in ['biased', 'unbiased']
    r = CORRELATION(X, maxlags=order, norm=norm)
    A, P, k = LEVINSON(r, allow_singularity=allow_singularity)
    return A, P, k




class pyule(ParametricSpectrum):
    """Class to create PSD based on the Yule Walker method
    
    See :func:`aryule` for description.
    
    .. plot::
        :width: 80%
        :include-source:
    
        from spectrum import *
        p = pyule(marple_data, 15, NFFT=4096)
        p()
        p.plot(sides='centerdc')
    
    
    """
    def __init__(self, data, order, norm='biased', NFFT=None, sampling=1., 
                 scale_by_freq=False):
        """**Constructor**
        
        For a detailled description of the parameters, see :func:`aryule`.
         
        :param array data:
        :param int order:
        :param int NFFT:
        :param float sampling:
        :param str norm: don't change if you do not know
                
        """
        super(pyule, self).__init__(data, ar_order=order,  NFFT=NFFT, 
                                    scale_by_freq=scale_by_freq)
        self.sampling = sampling
        self._norm_aryule = norm

    def __call__(self):
        import arma
        ar, rho, k = aryule(self.data, self.ar_order, norm=self._norm_aryule)
        psd = arma.arma2psd(ar, NPSD=self.NFFT, rho=rho, T=self.sampling)
        # save the AR and reflection coefficients.
        self.ar = ar
        self.reflection = k
        
        # save the PSD
        if self.datatype == 'real':
            import tools
            self.psd = tools.twosided_2_onesided(psd)
            #psd = res[0]
            #newpsd  = psd[0:self.NFFT/2]*2
            #newpsd[0] /= 2.
            #newpsd = numpy.append(newpsd, res[-1])
            #self.psd = newpsd
        else:
            self.psd = psd
        if self.scale_by_freq is True:
            self.scale()
        
    def _str_title(self):
        return "PYule PSD estimate\n"
    
    def __str__(self):
        return super(pyule, self).__str__()
    


