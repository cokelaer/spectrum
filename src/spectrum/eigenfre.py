import numpy
from numpy.fft import fft
from numpy.linalg import svd
from spectrum import default_NPSD
from psd import ParametricSpectrum
#import spectrum
from tools import twosided_2_onesided, centerdc_2_twosided, twosided_2_centerdc

verbose = False

__all__ = ["music", "ev", "pmusic", "eigen", "pev"]


class pmusic(ParametricSpectrum):
    """Class to create PSD using ARMA estimator.

    See :func:`pmusic` and :func:`eigenfre` for description.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        p = pmusic(marple_data, 15, NFFT=4096)
        p()
        p.plot()


    """
    def __init__(self, data, IP, NSIG=None, NFFT=None, sampling=1.):
        """**Constructor:**

        For a detailed description of the parameters, see :func:`arma_estimate`.

        :param array data:
        :param int IP:
        :param int NFFT:
        :param float sampling:

        """
        super(pmusic, self).__init__(data, ar_order=IP,  
                                            NFFT=NFFT, sampling=sampling)
        self.NSIG = NSIG

    def __call__(self):
        psd, _eigenvalues = eigen(self.data, self.ar_order, NSIG=self.NSIG, 
            NPSD=self.NFFT, threshold=None, criteria='aic', verbose=verbose, method='music')
        if self.datatype == 'real':
            self.psd = twosided_2_onesided(psd)
        else:
            self.psd = centerdc_2_twosided(psd)
        self.scale()
        

    def _str_title(self):
        return "Music PSD estimate\n"

    def __str__(self):
        return super(pmusic, self).__str__()

class pev(ParametricSpectrum):
    """Class to create PSD using ARMA estimator.

    See :func:`eigenfre` for description.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        p = pev(marple_data, 15, NFFT=4096)
        p()
        p.plot()


    """
    def __init__(self, data, IP, NSIG=None, NFFT=None, sampling=1.):
        """**Constructor:**

        For a detailed description of the parameters, see :func:`arma_estimate`.

        :param array data:
        :param int IP:
        :param int NFFT:
        :param float sampling:

        """
        super(pev, self).__init__(data, ar_order=IP,  
                                            NFFT=NFFT, sampling=sampling)
        self.NSIG = NSIG

    def __call__(self):
        psd, _eigenvalues = eigen(self.data, self.ar_order, NSIG=self.NSIG, 
            NPSD=self.NFFT, threshold=None, criteria='aic', verbose=verbose, 
            method='ev')
        if self.datatype == 'real':
            self.psd = twosided_2_onesided(psd)
        else:
            self.psd = centerdc_2_twosided(psd)
        self.scale()
        

    def _str_title(self):
        return "EV PSD estimate\n"

    def __str__(self):
        return super(pev, self).__str__()



def music(X, IP, NSIG=None, NPSD=default_NPSD, threshold=None, criteria='aic',
        verbose=False):
    """Eigen value pseudo spectrum estimate. See :func:`eigenfre`"""
    return eigen(X, IP, NSIG=NSIG, method='music', NPSD=NPSD, 
                 threshold=threshold, criteria=criteria, verbose=verbose)

def ev(X, IP, NSIG=None, NPSD=default_NPSD, threshold=None, criteria='aic',
        verbose=False):
    """Eigen value pseudo spectrum estimate. See :func:`eigenfre`"""
    return eigen(X, IP, NSIG=NSIG, method='ev', NPSD=NPSD, 
                 threshold=threshold, criteria=criteria, verbose=verbose)
#ev.__doc__ += "\n\n"+eigen.__doc__

def eigen(X, P, NSIG=None, method='music', threshold=None, NPSD=default_NPSD,
          criteria='aic', verbose=False):
    r"""Pseudo spectrum using eigenvector method (EV or Music)

    This function computes either the Music or EigenValue (EV) noise
    subspace frequency estimator.
    
    First, an autocorrelation matrix of order `P` is computed from 
    the data. Second, this matrix is separated into vector subspaces, 
    one a signal subspace and the other a noise
    subspace using a SVD method to obtain the eigen values and vectors. 
    From the eigen values :math:`\lambda_i`, and eigen vectors :math:`v_k`, 
    the **pseudo spectrum** (see note below) is computed as follows:
    
    .. math:: P_{ev}(f) = \frac{1}{e^H(f)\left(\sum\limits_{k=M+1}^{p} \frac{1}{\lambda_k}v_kv_k^H\right)e(f)}
    
    The separation of the noise and signal subspaces requires expertise
    of the signal. However, AIC and MDL criteria may be used to automatically
    perform this task.
    
    You still need to provide the parameter `P` to indicate the maximum number
    of eigen values to be computed. The criteria will just select a subset 
    to estimate the pseudo spectrum (see :func:`~spectrum.criteria.aic_eigen`
    and :func:`~spectrum.criteria.mdl_eigen` for details.
    
    .. note:: **pseudo spectrum**. func:`eigen` does not compute a PSD estimate.
        Indeed, the method does not preserve the measured process power.
    
    :param X: Array data samples
    :param int P: maximum number of eigen values to compute. NSIG (if 
        specified) must therefore be less than P. 
    :param str method: 'music' or 'ev'.
    :param int NSIG: If specified, the signal sub space uses NSIG eigen values.
    :param float threshold: If specified, the signal sub space is made of the
        eigen values larger than :math:`\rm{threshold} \times \lambda_{min}`,
        where :math:`\lambda_{min}` is the minimum eigen values.
    :param NPSD: :func:`spectrum.default_NPSD`

    :return:
        * PSD: Array of real frequency estimator values (two sided for 
                complex data and one sided for real data)
        * S, the eigen values

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import plot, log10, linspace, legend, axis

        psd, ev = pev(marple_data, 15, NSIG=11)
        f = linspace(-0.5, 0.5, len(psd))
        plot(f, 10 * log10(psd/max(psd)), label='User defined')
        psd, ev = pev(marple_data, 15, threshold=2)
        plot(f, 10 * log10(psd/max(psd)), label='threshold method (100)')
        psd, ev = pev(marple_data, 15)
        plot(f, 10 * log10(psd/max(psd)), label='AIC method (8)')
        legend()
        axis([-0.5, 0.5, -120, 0])


    .. seealso:: 
        :func:`pev`, 
        :func:`pmusic`, 
        :func:`~spectrum.criteria.aic_eigen`

    :References: [Marple]_, Chap 13
    
    .. todo:: for developers:
    
        * what should be the second argument of the criteria N, N-P, P...?
        * what should be the max value of NP
    """
    if method not in ['music', 'ev']:
        raise ValueError("method must be 'music' or 'ev'")
    
    if NSIG!=None and threshold!=None:
        raise ValueError("NSIG and threshold cannot be provided together")
    
    if NSIG != None:
        if NSIG < 0:
            raise ValueError('NSIG must be positive')
        if NSIG >= P:
            raise ValueError("NSIG must be stricly less than IP")
        
    #    
    N = len(X)
    NP = N - P
    
    assert 2 * NP > P-1, 'decrease the second argument' 
    if NP > 100:
        NP = 100
    FB = numpy.zeros((2*NP, P), dtype=complex)
    #FB = numpy.zeros((MAXU, IP), dtype=complex)
    Z = numpy.zeros(NPSD, dtype=complex)
    PSD = numpy.zeros(NPSD)
    
    
    # These loops can surely be replaced by a function that create such matrix
    for I in range(0, NP):
        for K in range(0, P):
            FB[I, K] = X[I-K+P-1]
            FB[I+NP, K] = X[I+K+1].conjugate()
    
    # This commented line produces the correct FB, as the 2 for loops above
    # It is more elegant but slower...corrmtx needs to be optimised (20/4/11)
    #FB2 = spectrum.linalg.corrmtx(X, P-1, method='modified')
    
    #Compute the eigen
    _U, S, V = svd (FB)
    # U and V are not the same as in Marple. Real or Imaginary absolute values
    # are correct but signs are not. This is wierd because the svd function
    # gives the same result as cvsd in Marple. Is FB correct ? it seems so.
    # The following operation has to be done. Otherwise, the resulting PSD is
    # not corect  
    V = -V.transpose()
    
    
    NSIG  = _get_signal_space(S, 2*NP, 
                             verbose=verbose, threshold=threshold, 
                             NSIG=NSIG, criteria=criteria)    


    #C   AI or Expert Knowledge to choose "signal" singular values, or input
    #C   NSIG at this point
    for I in range(NSIG, P):
        Z[0:P] = V[0:P, I]
        Z[P:NPSD] = 0

        Z  = fft(Z, NPSD)

        if method == 'music':  
            PSD = PSD + abs(Z)**2.
        elif method == 'ev' :  
            PSD = PSD + abs(Z)**2. / S[I]

    PSD = 1./PSD

    #for some reasons, we need to rearrange the output. this is related to 
    #the way U and V are order in the routine svd
    nby2 = NPSD/2
    newpsd = numpy.append(PSD[nby2:0:-1], PSD[nby2*2-1:nby2-1:-1])

    return newpsd, S

        
        
    
def _get_signal_space(S, NP, verbose=False, threshold=None, NSIG=None, 
                     criteria='aic'):    
    """todo
    
    
    """
    from criteria import aic_eigen, mdl_eigen
    # This section selects automatically the noise and signal subspaces.
    # NSIG being the number of eigenvalues corresponding to signals.
    if NSIG == None:
        if threshold == None:
            if verbose:
                print 'computing NSIG using AIC method'
            # get the minimum index of the AIC vector
            if criteria == 'aic':
                aic = aic_eigen(S, NP*2)
            elif criteria == 'mdl':
                aic = mdl_eigen(S, NP*2)
            # get the minimum index of the AIC vector, add 1 to get the NSIG
            NSIG = numpy.argmin(aic) + 1
            if verbose:print 'NSIG=', NSIG, ' found as the number of pertinent sinusoids'
        else:
            if verbose:
                print 'computing NSIG using user threshold '
            # following an idea from Matlab, pmusic, we look at the minimum
            # eigen value, and split the eigen values above and below
            # K times min eigen value, where K is >1
            m = threshold * min(S)
            new_s = S[numpy.where(S>m)]
            NSIG = len(new_s)
            if verbose:
                print 'found', NSIG
            if NSIG == 0:
                NSIG = 1
    return NSIG
