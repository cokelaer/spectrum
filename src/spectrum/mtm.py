"""Multitapering method.

This module provides a multi-tapering method for spectral estimation. 
The computation of the tapering windows being computationally expensive, a C 
module is used to compute the tapering window (see :func:`spectrum.mtm.dpss`). 

The estimation itself can be performed with the :func:`spectrum.mtm.pmtm` function


"""
import numpy as np
import ctypes
from ctypes import *
import ctypes as C
import platform
import os
from ctypes import POINTER
import os
from os.path import join as pj
from spectrum.tools import nextpow2
from pylab import semilogy

"""

from ctypes import *
my_sum=CDLL('sum.so')
a=numpy.array(range(10),dtype=float)
my_sum.sum(a.ctypes.data_as(c_void_p),int(10))


which you want to call from python. First make a shared object library by doing (at the command line)

gcc -c sum.c
gcc -shared -o sum.so sum.o

Note that on OSX -shared should be replaced by -dynamiclib and sum.so should be called sum.dylib Then you can do

"""

# Import shared mtspec library depending on the platform.
if platform.system() == 'Windows':
    lib_name = 'mydpss.pyd'
elif platform.system() == 'Darwin':
    lib_name = 'mydpss.so'
else:
    lib_name = 'mydpss.so'

p = os.path.abspath(os.path.dirname(__file__))
mtspeclib = ctypes.cdll.LoadLibrary(pj(p, lib_name))


def pmtm(x, NW=None, k=None, NFFT=None, e=None, v=None, method='adapt', show=True):
    """Multitapering spectral estimation

    :param array x: the data
    :param float NW: The time half bandwidth parameter (typical values are 2.5,3,3.5,4). Must be provided otherwise the tapering windows and eigen values (outputs of dpss) must be provided  
    :param int k: uses the first k Slepian sequences. If *k* is not provided, *k* is set to *NW*2*.
    :param NW
    :param e: the matrix containing the tapering windows
    :param v: the window concentrations (eigenvalues)
    :param str method: set how the eigenvalues are used. Must be in ['unity', 'adapt', 'eigen']
    :param bool show: plot results


    Usually in spectral estimation the mean to reduce bias is to use tapering window. 
    In order to reduce variance we need to average different spectrum. The problem 
    is that we have only one set of data. Thus we need to decompose a set into several 
    segments. Such method are well-known: simple daniell's periodogram, Welch's method
    and so on. The drawback of such methods is a loss of resolution since the segments 
    used to compute the spectrum are smaller than the data set. 
    The interest of multitapering method is to keep a good resolution while reducing 
    bias and variance. 

    How does it work? First we compute different simple periodogram with the whole data 
    set (to keep good resolution) but each periodgram is computed with a different 
    tapering windows. Then, we average all these spectrum. To avoid redundancy and bias
    due to the tapers mtm use special tapers.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *

        data = data_cosine(N=2048, A=0.1, sampling=1024, freq=200)
        # If you already have the DPSS windows
        [tapers, eigen] = dpss(2048, 2.5, 4)
        res = pmtm(data, e=tapers, v=eigen, show=False)
        # You do not need to compute the DPSS before end
        res = pmtm(data, NW=2.5, show=False)
        res = pmtm(data, NW=2.5, k=4, show=True)

    """
    assert method in ['adapt','eigen','unity']

    N = len(x)

    # if dpss not provided, compute them
    if e == None and v == None:
        if NW != None:
            [tapers, eigenvalues] = dpss(N, NW, k=k)
        else:
            raise ValueError("NW must be provided (e.g. 2.5, 3, 3.5, 4")
    elif e != None and v != None:
        eigenvalues = v[:]
        tapers = e[:]
    else:
        raise ValueError("if e provided, v must be provided as well and viceversa.")
    nwin = len(eigenvalues) # length of the eigen values vector to be used later

    # set the NFFT
    if NFFT==None:
        NFFT = max(256, 2**nextpow2(N))
    # si nfft smaller thqn N, cut otherwise add wero.
    # compute 
    if method in ['eigen', 'unity']:
        if method == 'unity':
            weights = np.ones((nwin, 1))
        elif method == 'eigen':
            # The S_k spectrum can be weighted by the eigenvalues, as in Park et al. 
            weights = np.array([_x/float(i+1) for i,_x in enumerate(eigenvalues)])
            weights = weights.reshape(nwin,1)

        Sk = abs(np.fft.fft(np.multiply(tapers.transpose(), x), NFFT))**2
        Sk = np.mean(Sk * weights, axis=0)


    elif method == 'adapt':
        # This version uses the equations from [2] (P&W pp 368-370).

        # Wrap the data modulo nfft if N > nfft
        sig2 = np.dot(x, x) / float(N)
        Sk = abs(np.fft.fft(np.multiply(tapers.transpose(), x), NFFT))**2
        Sk = Sk.transpose()
        S = (Sk[:,0] + Sk[:,1]) / 2    # Initial spectrum estimate
        S = S.reshape(NFFT, 1)
        Stemp = np.zeros((NFFT,1))
        S1 = np.zeros((NFFT,1))
        # Set tolerance for acceptance of spectral estimate:
        tol = 0.0005 * sig2 / float(NFFT)
        i = 0
        a = sig2 * (1 - eigenvalues)

        # converges very quickly but for safety; set i<100
        while sum(np.abs(S-S1))/NFFT > tol and i<100:
            i = i + 1
            # calculate weights
            b1 = np.multiply(S, np.ones((1,nwin))) 
            b2 = np.multiply(S,eigenvalues.transpose()) + np.ones((NFFT,1))*eigenvalues.transpose()
            b = b1/b2
            
            # calculate new spectral estimate
            wk=(b**2)*(np.ones((NFFT,1))*eigenvalues.transpose())
            S1 = sum(wk.transpose()*Sk.transpose())/ sum(wk.transpose())
            S1 = S1.reshape(NFFT, 1)
            Stemp = S1
            S1 = S
            S = Stemp   # swap S and S1
        Sk = S
    #clf(); p.plot(); plot(arange(0,0.5,1./512),20*log10(res[0:256]))
    if show==True:
        semilogy(Sk)
    return Sk



def dpss(N, NW=None, k=None):
    r"""Discrete prolate spheroidal (Slepian) sequences

    Calculation of the Discrete Prolate Spheroidal Sequences also known as the
    slepian sequences, and the corresponding eigenvalues. 

    :param int N: desired window length
    :param float NW: The time half bandwidth parameter (typical values are 2.5,3,3.5,4). 
    :param int k: returns the first k Slepian sequences. If *k* is not provided, *k* is set to *NW*2*.
    :return:
        * tapers, a matrix of tapering windows. Matrix is a N by *k* (k is the number of windows)
        * eigen, a vector of eigenvalues of length *k*

    The discrete prolate spheroidal or Slepian sequences derive from the following 
    time-frequency concentration problem. For all finite-energy sequences index 
    limited to some set , which sequence maximizes the following ratio:

    .. math::

        \lambda = \frac{\int_{-W}^{W}\left| X(f) \right|^2 df} 
            {\int_{-F_s/2}^{F_s/2}\left| X(f) \right|^2 df}
   
    where :math:`F_s` is the sampling frequency and :math:`|W| < F_s/2`. 
    This ratio determines which index-limited sequence has the largest proportion of its 
    energy in the band :math:`[-W,W]` with :math:`0  < \lambda < 1`. 
    The sequence maximizing the ratio is the first 
    discrete prolate spheroidal or Slepian sequence. The second Slepian sequence  
    maximizes the ratio and is orthogonal to the first Slepian sequence. The third 
    Slepian sequence maximizes the ratio of integrals and is orthogonal to both 
    the first and second Slepian sequences and so on.

    .. note:: Note about the implementation. Since the slepian generation is computationally expensive, 
        we use a C implementation based on the C code written by Lees as published in:

            Lees, J. M. and J. Park (1995): Multiple-taper spectral analysis: A stand-alone 
            C-subroutine: Computers & Geology: 21, 199-236.

        However, the original C code has been trimmed. Indeed, we only require the multitap 
        function (that depends on jtridib, jtinvit functions only).

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import *
        N = 512
        [w, eigens] = dpss(N, 2.5, 4)
        plot(w)
        title('Slepian Sequences N=%s, NW=2.5' % N)
        axis([0, N, -0.15, 0.15])
        legend(['1st window','2nd window','3rd window','4th window'])

    Windows are normalised:

    .. math::  \sum_k h_k h_k = 1
     
    :references: [Percival]_

        Slepian, D. Prolate spheroidal wave functions, Fourier analysis, and
        uncertainty V: The discrete case. Bell System Technical Journal,
        Volume 57 (1978), 1371430
    
    .. note:: the C code to create the slepian windows is extracted from original C code 
        from Lees and Park (1995) and uses the conventions of Percival and Walden (1993). 
        Functions that are not used here were removed.

    """
    assert NW<N/2 , "NW must be stricly less than N/2"
    if k == None:
        k = min(round(2*NW),N)
        k = int(max(k,1))
    from numpy import dot, zeros, arange, sqrt
    mtspeclib.multitap.restype = None

    lam = zeros(k, dtype=float)
    tapers = zeros(k*N, dtype=float)
    tapsum = zeros(k, dtype=float)

    res = mtspeclib.multitap(
        c_int(N),
        c_int(k),
        lam.ctypes.data_as(c_void_p), 
        c_float(NW), 
        tapers.ctypes.data_as(c_void_p), 
        tapsum.ctypes.data_as(c_void_p),
        )

    # normalisation by sqtr(N). It is required to have normalised windows 
    tapers = tapers.reshape(k,N).transpose() / sqrt(N)

    for i in range(k):
        # By convention (Percival and Walden, 1993 pg 379)
        # * symmetric tapers (k=0,2,4,...) should have a positive average.
        # * antisymmetric tapers should begin with a positive lobe
        if i%2 == 0:
            if tapsum[i]<0: 
                tapsum[i] *= -1
                tapers[:,i] *= -1
        else:
            if tapers[0,i] < 0:
                tapsum[i] *= -1
                tapers[:,i] *= -1
    
    # Now find the eigenvalues of the original 
    # Use the autocovariance sequence technique from Percival and Walden, 1993
    # pg 390 to get the eigenvalues more precisely (same as matlab output)

    # The values returned in lam are not exacly the same as in the following methods.
    acvs = _autocov(tapers.transpose(), debias=False) * N
    nidx = arange(N)
    W = float(NW)/N
    r = 4*W*np.sinc(2*W*nidx)
    r[0] = 2*W
    eigvals = dot(acvs, r)

    #return (tapers, lam)
    return (tapers, eigvals)




def _other_dpss_method(N, NW, Kmax):
    """Returns the Discrete Prolate Spheroidal Sequences of orders [0,Kmax-1]
    for a given frequency-spacing multiple NW and sequence length N. 

    See dpss function that is the official version. This version is indepedant 
    of the C code and relies on Scipy function. However, it is slower by a factor 3

    Tridiagonal form of DPSS calculation from:

    """
    # here we want to set up an optimization problem to find a sequence
    # whose energy is maximally concentrated within band [-W,W].
    # Thus, the measure lambda(T,W) is the ratio between the energy within
    # that band, and the total energy. This leads to the eigen-system
    # (A - (l1)I)v = 0, where the eigenvector corresponding to the largest
    # eigenvalue is the sequence with maximally concentrated energy. The
    # collection of eigenvectors of this system are called Slepian sequences,
    # or discrete prolate spheroidal sequences (DPSS). Only the first K,
    # K = 2NW/dt orders of DPSS will exhibit good spectral concentration
    # [see http://en.wikipedia.org/wiki/Spectral_concentration_problem]
    
    # Here I set up an alternative symmetric tri-diagonal eigenvalue problem
    # such that
    # (B - (l2)I)v = 0, and v are our DPSS (but eigenvalues l2 != l1)
    # the main diagonal = ([N-1-2*t]/2)**2 cos(2PIW), t=[0,1,2,...,N-1]
    # and the first off-diangonal = t(N-t)/2, t=[1,2,...,N-1]
    # [see Percival and Walden, 1993]
    from scipy import linalg as la
    Kmax = int(Kmax)
    W = float(NW)/N
    ab = np.zeros((2,N), 'd')
    nidx = np.arange(N)
    ab[0,1:] = nidx[1:]*(N-nidx[1:])/2.
    ab[1] = ((N-1-2*nidx)/2.)**2 * np.cos(2*np.pi*W)
    # only calculate the highest Kmax-1 eigenvectors
    l,v = la.eig_banded(ab, select='i', select_range=(N-Kmax, N-1))
    dpss = v.transpose()[::-1]

    # By convention (Percival and Walden, 1993 pg 379)
    # * symmetric tapers (k=0,2,4,...) should have a positive average.
    # * antisymmetric tapers should begin with a positive lobe
    fix_symmetric = (dpss[0::2].sum(axis=1) < 0)
    for i, f in enumerate(fix_symmetric):
        if f:
            dpss[2*i] *= -1
    fix_skew = (dpss[1::2,1] < 0)
    for i, f in enumerate(fix_skew):
        if f:
            dpss[2*i+1] *= -1

    # Now find the eigenvalues of the original 
    # Use the autocovariance sequence technique from Percival and Walden, 1993
    # pg 390
    # XXX : why debias false? it's all messed up o.w., even with means
    # on the order of 1e-2
    acvs = _autocov(dpss, debias=False) * N
    r = 4*W*np.sinc(2*W*nidx)
    r[0] = 2*W
    eigvals = np.dot(acvs, r)
    return dpss, eigvals


def _autocov(s, **kwargs):
    """Returns the autocovariance of signal s at all lags.

    Adheres to the definition
    sxx[k] = E{S[n]S[n+k]} = cov{S[n],S[n+k]}
    where E{} is the expectation operator, and S is a zero mean process
    """
    # only remove the mean once, if needed
    debias = kwargs.pop('debias', True)
    axis = kwargs.get('axis', -1)
    if debias:
        s = remove_bias(s, axis)
    kwargs['debias'] = False
    return _crosscov(s, s, **kwargs)


def _crosscov(x, y, axis=-1, all_lags=False, debias=True):
    """Returns the crosscovariance sequence between two ndarrays.
    This is performed by calling fftconvolve on x, y[::-1]

    Parameters


    x: ndarray
    y: ndarray
    axis: time axis

    all_lags: {True/False}
    whether to return all nonzero lags, or to clip the length of s_xy
    to be the length of x and y. If False, then the zero lag covariance
    is at index 0. Otherwise, it is found at (len(x) + len(y) - 1)/2

    debias: {True/False}
    Always removes an estimate of the mean along the axis, unless
    told not to.


    cross covariance is defined as
    sxy[k] := E{X[t]*Y[t+k]}, where X,Y are zero mean random processes
    """
    if x.shape[axis] != y.shape[axis]:
        raise ValueError(
            'crosscov() only works on same-length sequences for now'
            )
    if debias:
        x = _remove_bias(x, axis)
        y = _remove_bias(y, axis)
    slicing = [slice(d) for d in x.shape]
    slicing[axis] = slice(None,None,-1)
    sxy = _fftconvolve(x, y[tuple(slicing)], axis=axis, mode='full')
    N = x.shape[axis]
    sxy /= N
    if all_lags:
        return sxy
    slicing[axis] = slice(N-1,2*N-1)
    return sxy[tuple(slicing)]
    
def _crosscorr(x, y, **kwargs):
    """
    Returns the crosscorrelation sequence between two ndarrays.
    This is performed by calling fftconvolve on x, y[::-1]

    Parameters


    x: ndarray
    y: ndarray
    axis: time axis
    all_lags: {True/False}
    whether to return all nonzero lags, or to clip the length of r_xy
    to be the length of x and y. If False, then the zero lag correlation
    is at index 0. Otherwise, it is found at (len(x) + len(y) - 1)/2

    Notes


    cross correlation is defined as
    rxy[k] := E{X[t]*Y[t+k]}/(E{X*X}E{Y*Y})**.5,
    where X,Y are zero mean random processes. It is the noramlized cross
    covariance.
    """
    sxy = _crosscov(x, y, **kwargs)
    # estimate sigma_x, sigma_y to normalize
    sx = np.std(x)
    sy = np.std(y)
    return sxy/(sx*sy)

def _remove_bias(x, axis):
    "Subtracts an estimate of the mean from signal x at axis"
    padded_slice = [slice(d) for d in x.shape]
    padded_slice[axis] = np.newaxis
    mn = np.mean(x, axis=axis)
    return x - mn[tuple(padded_slice)]


def _fftconvolve(in1, in2, mode="full", axis=None):
    """ Convolve two N-dimensional arrays using FFT. See convolve.

    This is a fix of scipy.signal.fftconvolve, adding an axis argument and
    importing locally the stuff only needed for this function
    
    """
    #Locally import stuff only required for this:
    from scipy.fftpack import fftn, fft, ifftn, ifft
    from scipy.signal.signaltools import _centered
    from numpy import array, product


    s1 = array(in1.shape)
    s2 = array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex) or
                      np.issubdtype(in2.dtype, np.complex))

    if axis is None:
        size = s1+s2-1
        fslice = tuple([slice(0, int(sz)) for sz in size])
    else:
        equal_shapes = s1==s2
        # allow equal_shapes[axis] to be False
        equal_shapes[axis] = True
        assert equal_shapes.all(), 'Shape mismatch on non-convolving axes'
        size = s1[axis]+s2[axis]-1
        fslice = [slice(l) for l in s1]
        fslice[axis] = slice(0, int(size))
        fslice = tuple(fslice)

    # Always use 2**n-sized FFT
    fsize = 2**np.ceil(np.log2(size))
    if axis is None:
        IN1 = fftn(in1,fsize)
        IN1 *= fftn(in2,fsize)
        ret = ifftn(IN1)[fslice].copy()
    else:
        IN1 = fft(in1,fsize,axis=axis)
        IN1 *= fft(in2,fsize,axis=axis)
        ret = ifft(IN1,axis=axis)[fslice].copy()
    if not complex_result:
        del IN1
        ret = ret.real
    if mode == "full":
        return ret
    elif mode == "same":
        if product(s1,axis=0) > product(s2,axis=0):
            osize = s1
        else:
            osize = s2
        return _centered(ret,osize)
    elif mode == "valid":
        return _centered(ret,abs(s2-s1)+1)
