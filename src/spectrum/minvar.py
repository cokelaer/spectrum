"""
.. topic:: Minimum Variance Spectral Estimators

    
    .. autosummary::

        :mod:`spectrum.minvar`
        pminvar
        
    .. codeauthor:: Thomas Cokelaer, 2011
"""
import numpy
from burg import arburg
from numpy.fft import fft
from psd import ParametricSpectrum
from spectrum import default_NPSD
import errors

__all__ = ["minvar", "pminvar"]


def minvar(X, order, sampling=1., NFFT=default_NPSD):
    r"""Minimum Variance Spectral Estimation (MV)

    This function computes the minimum variance spectral estimate using
    the Musicus procedure.  The Burg algorithm from :func:`~spectrum.burg.arburg`
    is used for the estimation of the autoregressive parameters. 
    The MV spectral estimator is given by:
    
    .. math:: P_{MV}(f) = \frac{T}{e^H(f) R^{-1}_p e(f)}
    
    
    where :math:`R^{-1}_p` is the inverse of the estimated autocorrelation
    matrix  (Toeplitz) and :math:`e(f)` is the complex sinusoid vector. 
    
    :param X: Array of complex or real data samples (length N)
    :param int order: Dimension of correlation matrix (AR order = order - 1 )
    :param float T: Sample interval (PSD scaling)
    :param int NFFT: length of the final PSD

    :return:
        * PSD  - Power spectral density values (two-sided)
        * AR   - AR coefficients (Burg algorithm)
        * k    - Reflection coefficients (Burg algorithm)

    .. note:: The MV spectral estimator is not a true PSD function because the
        area under the MV estimate does not represent the total power in the 
        measured process. MV minimises the variance of the output of a narrowband
        filter and adpats itself to the spectral content of the input data 
        at each frequency.  
         
    :Example: The following example computes a PSD estimate using :func:`minvar`  
        The output PSD is transformed to a ``centerdc`` PSD and plotted.
        
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import plot, log10, linspace, xlim
        psd, A, k = minvar(marple_data, 15)
        psd = twosided_2_centerdc(psd) # switch positive and negative freq
        f = linspace(-0.5, 0.5, len(psd))
        plot(f, 10 * log10(psd/max(psd)))
        xlim(-0.5, 0.5 )

    .. seealso::
        
        * External functions used are :meth:`~spectrum.burg.arburg` 
          and  numpy.fft.fft
        * :class:`pminvar`, a Class dedicated to MV method.

    :Reference: [Marple]_
    
    """
    errors.is_positive_integer(order)
    errors.is_positive_integer(NFFT)
    
    psi = numpy.zeros(NFFT, dtype=complex)
    
    # First, we need to compute the AR values (note that order-1)
    A, P, k = arburg (X, order - 1)
    # add the order 0 
    A = numpy.insert(A, 0, 1.+0j)

    # We cannot compare the output with those of MARPLE in a precise way.
    # Indeed the burg algorithm is only single precision in fortram code
    # So, the AR values are slightly differnt. 
    # The followign values are those from Marple     
    """A[1] = 2.62284255-0.701703191j
    A[2] = 4.97930574-2.32781982j 
    A[3] = 6.78445101-5.02477741j
    A[4] =7.85207081-8.01284409j
    A[5] =7.39412165-10.7684202j 
    A[6] =6.03175116-12.7067814j
    A[7] =3.80106878-13.6808891j
    A[8] =1.48207295-13.2265558j
    A[9] =-0.644280195-11.4574194j
    A[10] =-2.02386642-8.53268814j
    A[11] =-2.32437634-5.25636244j 
    A[12] =-1.75356281-2.46820402j
    A[13] =-0.888899028-0.781434655j 
    A[14] =-0.287197977-0.0918145925j
    P = 0.00636525545
    """
    
    # if we use exqtly the same AR coeff and P from MArple BUrg output, then 
    # we can compare the following code. This has been done and reveals that
    # the FFT in marple is also slightly different (precision) from this one.
    # However, the results are sufficiently close (when NPSD is small) that
    # we are confident the following code is correct.

    # Compute the psi coefficients
    for K in range(0, order):
        SUM = 0.
        MK = order-K

        #  Correlate the autoregressive parameters
        for I in range(0, order - K):
            SUM = SUM + float(MK-2*I) * A[I].conjugate()*A[I+K]  # Eq. (12.25)
        
        SUM = SUM/P
        if K != 0:
            psi[NFFT-K] = SUM.conjugate()
        psi[K] = SUM

    # Compute FFT of denominator
    psi = fft(psi, NFFT)          
    
    #  Invert the psi terms at this point to get PSD values
    PSD = sampling / numpy.real(psi)
    
    return PSD, A, k



class pminvar(ParametricSpectrum):
    """Class to create PSD based on the Minimum variance spectral estimation
    
    See :func:`minvar` for description.
    
    .. plot::
        :width: 80%
        :include-source:
    
        from spectrum import *
        p = pminvar(marple_data, 15, NFFT=4096)
        p()
        p.plot(sides='centerdc')
    
    
    """
    def __init__(self, data, order, NFFT=None, sampling=1.):
        """**Constructor**
        
        For a detailled description of the parameters, see :func:`minvar`.
         
        :param array data:
        :param int order:
        :param int NFFT:
        :param float sampling:
                
        """
        super(pminvar, self).__init__(data, ar_order=order, sampling=sampling, 
                                            NFFT=NFFT)
        

    def __call__(self):
        res = minvar(self.data, self.ar_order, sampling=self.sampling, 
                     NFFT=self.NFFT)
        
        # save the AR and reflection coefficients.
        self.ar = res[1]
        self.reflection = res[2]
        
        # save the PSD
        if self.datatype == 'real':
            import tools
            self.psd = tools.twosided_2_onesided(res[0])
            #psd = res[0]
            #newpsd  = psd[0:self.NFFT/2]*2
            #newpsd[0] /= 2.
            #newpsd = numpy.append(newpsd, res[-1])
            #self.psd = newpsd
        else:
            self.psd = res[0]
        self.scale()
        
    def _str_title(self):
        return "Minimum Variance spectral estimation\n"
    
    def __str__(self):
        return super(pminvar, self).__str__()
    


