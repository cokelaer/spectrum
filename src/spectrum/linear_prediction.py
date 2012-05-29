"""Linear prediction tools

:References: [Kay]_
"""

from levinson import LEVINSON, rlevinson
import numpy
from scipy.signal import deconvolve
__all__ = ['ac2poly', 'poly2ac', 'ac2rc', 'rc2poly', 'rc2ac', 'is2rc', 'rc2is', 
    'rc2lar', 'lar2rc', 'poly2rc', 'lsf2poly', 'poly2lsf']

"""
lpc Linear prediction filter coefficients
schurrc Compute reflection coefficients from autocorrelation sequence
"""


def ac2poly(data):
    """Convert autocorrelation sequence to prediction polynomial

    :param array data:    input data (list or numpy.array)
    :return: 
        * AR parameters
        * noise variance

    This is an alias to::

        a, e, c = LEVINSON(data)

    :Example:

    .. doctest::

        >>> r = [5, -2, 1]
        >>> ar, e = ac2poly(r)
        >>> ar
        array([ 1. ,  0. , -0.2])
        >>> e
        4.7999999999999998

    """
    a, e, _c = LEVINSON(data)
    a = numpy.insert(a, 0, 1) 
    return a, e


def ac2rc(data):
    """Convert autocorrelation sequence to reflection coefficients

    :param data: an autorrelation vector
    :return: the reflection coefficient and data[0]

    This is an alias to::

        a, e, c = LEVINSON(data)
        c, data[0]

    """
    a, e, _c = LEVINSON(data)
    return a, data[0]
    
    
def poly2ac(poly, efinal):
    """ Convert prediction filter polynomial to autocorrelation sequence

    :param array poly: the AR parameters
    :param efinal: an estimate of the final error
    :return: the autocorrelation  sequence

    .. doctest::

        >>> ar = [ 1. ,  0. , -0.2]
        >>> efinal = 4.8
        >>> r = poly2ac(ar, efinal)

    """
    results = rlevinson(poly, efinal)
    return results[0] 
    
def ar2rc(ar):
    """ Convert autoregressive parameters into reflection coefficients """
    raise NotImplementedError
    
    
def poly2rc(a, efinal):
    """Convert prediction filter polynomial to reflection coefficients

    :param a: AR parameters
    :param efinal: 
    """
    results = rlevinson(a, efinal)
    return results[2]
    
def rc2poly(kr, r0=None):
    """convert reflection coefficients to prediction filter polynomial

    :param k: reflection coefficients




    """
    # Initialize the recursion
    from levinson import levup
    p = len(kr)              #% p is the order of the prediction polynomial.
    a = numpy.array([1, kr[0]])           #% a is a true polynomial.
    e = numpy.zeros(5)

    if r0 == None:
        e0 = 0
    else:
        e0 = r0

    e[0] = e0 * (1. - numpy.conj(numpy.conjugate(kr[0])*kr[0]))

    # Continue the recursion for k=2,3,...,p, where p is the order of the 
    # prediction polynomial.

    for k in range(1, p):
        [a, e[k]] = levup(a, kr[k], e[k-1])

    efinal = e[-1]
    return a, efinal


def rc2ac(k,R0):
    """Convert reflection coefficients to autocorrelation sequence.

    :param k: reflection coefficients
    :param R0: zero-lag autocorrelation
    :returns: the autocorrelation sequence

    .. seealso:: :func:`ac2rc`, :func:`poly2rc`, :func:`ac2poly`, :func:`poly2rc`, :func:`rc2poly`.

    """
    [a,efinal] = rc2poly(k, R0)
    R, u, kr, e = rlevinson(a, efinal)
    return R


def is2rc(inv_sin):
    """Convert inverse sine parameters to reflection coefficients.

    :param inv_sin: inverse sine parameters
    :return: reflection coefficients
    
    .. seealso::  :func:`rc2is`, :func:`poly2rc`, :func:`ac2rc`, :func:`lar2rc`.

    :Reference: J.R. Deller, J.G. Proakis, J.H.L. Hansen, "Discrete-Time Processing of Speech Signals", Prentice Hall, Section 7.4.5.

    """
    return numpy.sin(numpy.array(inv_sin)*numpy.pi/2);



def  rc2is(k):
    """Convert reflection coefficients to inverse sine parameters.

    :param k: reflection coefficients
    :return: inverse sine parameters

    .. seealso:: :func:`is2rc`, :func:`rc2poly`, :func:`rc2acC`, :func:`rc2lar`.

    Reference: J.R. Deller, J.G. Proakis, J.H.L. Hansen, "Discrete-Time 
       Processing of Speech Signals", Prentice Hall, Section 7.4.5.

    """
    assert numpy.isrealobj(k), 'Inverse sine parameters not defined for complex reflection coefficients.'
    if max(numpy.abs(k)) >= 1:
        raise ValueError('All reflection coefficients should have magnitude less than unity.')

    return (2/numpy.pi)*numpy.arcsin(k)

def rc2lar(k):
    """Convert reflection coefficients to log area ratios.

    :param k: reflection coefficients
    :return: inverse sine parameters

    The log area ratio is defined by G = log((1+k)/(1-k)) , where the K is the reflection coefficient.

    .. seealso:: :func:`lar2rc`, :func:`rc2poly`, :func:`rc2ac`, :func:`rc2ic`.

    :References:
       [1] J. Makhoul, "Linear Prediction: A Tutorial Review," Proc. IEEE, Vol.63, No.4, pp.561-580, Apr 1975.

    """
    assert numpy.isrealobj(k), 'Log area ratios not defined for complex reflection coefficients.'
    if max(numpy.abs(k)) >= 1:
        raise ValueError('All reflection coefficients should have magnitude less than unity.')

    # Use the relation, atanh(x) = (1/2)*log((1+k)/(1-k))
    return -2 * numpy.arctanh(-numpy.array(k))



def lar2rc(g):
    """Convert log area ratios to reflection coefficients.

    :param g:  log area ratios
    :returns: the reflection coefficients

    .. seealso: :func:`rc2lar`, :func:`poly2rc`, :func:`ac2rc`, :func:`is2rc`.

    :References:
       [1] J. Makhoul, "Linear Prediction: A Tutorial Review," Proc. IEEE,  Vol.63, No.4, pp.561-580, Apr 1975.

    """
    assert numpy.isrealobj(g), 'Log area ratios not defined for complex reflection coefficients.'
    # Use the relation, tanh(x) = (1-exp(2x))/(1+exp(2x))
    return -numpy.tanh(-numpy.array(g)/2)



def lsf2poly(lsf):
    """Convert line spectral frequencies to prediction filter coefficients

    returns a vector a containing the prediction filter coefficients from a vector lsf of line spectral frequencies.

    .. doctest::

        >>> lsf = [0.7842 ,   1.5605  ,  1.8776 ,   1.8984,    2.3593]
        >>> a = lsf2poly(lsf)
        array([  1.00000000e+00,   6.14837835e-01,   9.89884967e-01,
            9.31594056e-05,   3.13713832e-03,  -8.12002261e-03 ])

    .. seealso:: poly2lsf, rc2poly, ac2poly, rc2is
    """
    #   Reference: A.M. Kondoz, "Digital Speech: Coding for Low Bit Rate Communications
    #   Systems" John Wiley & Sons 1994 ,Chapter 4 

    # Line spectral frequencies must be real.

    lsf = numpy.array(lsf)

    if max(lsf) > numpy.pi or min(lsf) < 0:
        raise ValueError('Line spectral frequencies must be between 0 and pi.')

    p = len(lsf) # model order

    # Form zeros using the LSFs and unit amplitudes
    z  = numpy.exp(1.j * lsf)

    # Separate the zeros to those belonging to P and Q
    rQ = z[0::2]
    rP = z[1::2]
    
    # Include the conjugates as well
    rQ = numpy.concatenate((rQ, rQ.conjugate()))
    rP = numpy.concatenate((rP, rP.conjugate()))
    
    # Form the polynomials P and Q, note that these should be real
    Q  = numpy.poly(rQ);
    P  = numpy.poly(rP);
    
    # Form the sum and difference filters by including known roots at z = 1 and
    # z = -1 
    
    if p%2:
        # Odd order: z = +1 and z = -1 are roots of the difference filter, P1(z)
        P1 = numpy.convolve(P, [1, 0, -1])
        Q1 = Q
    else:
        # Even order: z = -1 is a root of the sum filter, Q1(z) and z = 1 is a
        # root of the difference filter, P1(z)
        P1 = numpy.convolve(P, [1, -1])
        Q1 = numpy.convolve(Q, [1,  1])

    # Prediction polynomial is formed by averaging P1 and Q1

    a = .5 * (P1+Q1)
    return a[0:-1:1] # do not return last element 


def poly2lsf(a):
    """Prediction polynomial to line spectral frequencies.

    converts the prediction polynomial specified by A,
    into the corresponding line spectral frequencies, LSF. 
    normalizes the prediction polynomial by A(1).

    .. doctest::

        >>> a = [1.0000  ,  0.6149   , 0.9899   , 0.0000 ,   0.0031,   -0.0082
        >>> lsf = poly2lsf(a)
        >>> lsf =  array([  0.7842,    1.5605 ,   1.8776 ,   1.8984,    2.3593])

    .. seealso:: lsf2poly, poly2rc, poly2qc, rc2is
    """

    #Line spectral frequencies are not defined for complex polynomials.

    # Normalize the polynomial

    a = numpy.array(a)
    if a[0] != 1:
        a/=a[0]

    if max(numpy.abs(numpy.roots(a))) >= 1.0:
        error('The polynomial must have all roots inside of the unit circle.');


    # Form the sum and differnce filters

    p  = len(a)-1   # The leading one in the polynomial is not used
    a1 = numpy.concatenate((a, numpy.array([0])))        
    a2 = a1[-1::-1]
    P1 = a1 - a2        # Difference filter
    Q1 = a1 + a2        # Sum Filter 

    # If order is even, remove the known root at z = 1 for P1 and z = -1 for Q1
    # If odd, remove both the roots from P1

    if p%2: # Odd order
        P, r = deconvolve(P1,[1, 0 ,-1])
        Q = Q1
    else:          # Even order 
        P, r = deconvolve(P1, [1, -1])
        Q, r = deconvolve(Q1, [1,  1])
    
    rP  = numpy.roots(P)
    rQ  = numpy.roots(Q)

    aP  = numpy.angle(rP[1::2])
    aQ  = numpy.angle(rQ[1::2])

    lsf = sorted(numpy.concatenate((-aP,-aQ)))

    return lsf


def schurrc():
    raise NotImplementedError



