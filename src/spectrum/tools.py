"""
.. topic:: Tools module

    .. autosummary::

        db2mag
        db2pow
        mag2db
        nextpow2
        pow2db
        onesided_2_twosided
        twosided_2_onesided
        centerdc_2_twosided
        twosided_2_centerdc
        
    .. codeauthor:: Thomas Cokelaer, 2011
"""
import numpy
from numpy import ceil, log2
from collections import deque


#__all__ = ["cshift", "pow2dB", "nextpow2", "twosided", "twosided_zerolag", 
#           "fftshift"]


def fftshift(x):
    """wrapper to numpy.fft.fftshift
    
    .. doctest::
    
        >>> from spectrum import fftshift
        >>> x = [100, 2, 3, 4, 5]
        >>> fftshift(x)
        array([  4,   5, 100,   2,   3])
    """
    return numpy.fft.fftshift(x)
    
def _swapsides(data):
    """todo is it really useful ? 
    
    Swap  sides
    
    .. doctest::
    
        >>> from spectrum import swapsides
        >>> x = [-2, -1, 1, 2]
        >>> swapsides(x)
        array([ 2, -2, -1])
    
    """
    N = len(data)
    from numpy import concatenate
    return concatenate((data[N/2+1:], data[0:N/2]))

def twosided_2_onesided(data):
    """Convert a one-sided PSD to a twosided PSD

    In order to keep the power in the onesided PSD the same
    as in the twosided version, the onesided values are twice
    as much as in the input data (except for the zero-lag value). 

    ::

        >>> twosided_2_onesided([10, 2,3,3,2,8])
        array([ 10.,   4.,   6.,   8.])

    """
    assert len(data)%2 == 0
    N = len(data) 
    psd = numpy.array(data[0:N/2+1])*2.
    psd[0]/=2.
    psd[-1] = data[-1]
    return psd

def onesided_2_twosided(data):
    """Convert a two-sided PSD to a one-sided PSD

    In order to keep the power in the twosided PSD the same
    as in the onesided version, the twosided values are 2 times
    lower than the input data (except for the zero-lag and N-lag 
    values).

    ::

        >>> twosided_2_onesided([10, 4, 6, 8])
        array([ 10.,   2.,   3.,   3., 2., 8.])

    """
    psd = numpy.concatenate((data[0:-1], cshift(data[-1:0:-1], -1)))/2.
    psd[0]*=2.
    psd[-1]*=2.
    return psd

def twosided_2_centerdc(data):
    """Convert a two-sided PSD to a center-dc PSD"""
    N = len(data)
    newpsd = numpy.concatenate((cshift(data[N/2:], 1), data[0:N/2]))
    newpsd[0] = data[-1]
    return newpsd

def centerdc_2_twosided(data):
    """Convert a center-dc PSD to a twosided PSD"""
    N = len(data)
    newpsd = numpy.concatenate((data[N/2:], (cshift(data[0:N/2], -1))))
    return newpsd


def twosided(data):
    """return a twosided vector with non-duplication of the first element 

    .. doctest::

        >>> from spectrum import twosided
        >>> a = [1,2,3]
        >>> twosided(a)
        array([3, 2, 1, 2, 3])

    """
    twosided = numpy.concatenate((data[::-1], data[1:]))
    return twosided #remove the first element to have a power of 2 and compatiable with pylab.psd


def _twosided_zerolag(data, zerolag):
    """Build a symmetric vector out of stricly positive lag vector and zero-lag  
    
    .. doctest::
    
        >>> data = [3,2,1]
        >>> zerolag = 4
        >>> twosided_zerolag(data, zerolag)
        array([1, 2, 3, 4, 3, 2, 1])
    
    .. seealso:: Same behaviou as :func:`twosided_zerolag`
    """
    res = twosided(numpy.insert(data, 0, zerolag))
    return res


def cshift(data, offset):
    """Circular shift to the right (within an array) by a given offset

    :param array data: input data (list or numpy.array)
    :param int offset: shift the array with the offset
    
    .. doctest::

        >>> from spectrum import cshift
        >>> cshift([0, 1, 2, 3, -2, -1], 2)
        array([-2, -1,  0,  1,  2,  3])

    """
    # the deque method is suppose to be optimal when using rotate to shift the
    # data that playing with the data to build a new list.
    a = deque(data)
    a.rotate(offset)
    return numpy.array(a)  #convert back to an array. Is it necessary?


def pow2db(x):
    """returns the corresponding decibel (dB) value for a power value x. 

    The relationship between power and decibels is:

    .. math::    X_{dB} = 10 * \log_{10}(x)

    .. doctest::

        >>> x = pow2db(0.1)
        >>> x
        -10.0
    """
    return 10 * numpy.log10(x)


def db2pow(xdb):
    """Convert decibels (dB) to power
    
    .. doctest::
    
        >>> p = db2pow(-10)
        >>> p
        0.1
        
    .. seealso:: :func:`pow2db`
    """
    return 10.**(xdb/10.)
     

def nextpow2(x):
    """returns the smallest power of two that is greater than or equal to the 
    absolute value of x.

    This function is useful for optimizing FFT operations, which are 
    most efficient when sequence length is an exact power of two.

    :Example:

    .. doctest::
    
        >>> x = [255, 256, 257]
        >>> nextpow2(x)
        array([8, 8, 9])

    """
    res = ceil(log2(x))
    return res.astype('int')  #we want integer values only but ceil gives float


def db2mag(xdb):
    """Convert decibels (dB) to magnitude
    
    .. doctest::
    
        >>> db2mag(-20)
        0.1
        
    .. seealso:: :func:`pow2db`
    """
    return 10.**(xdb/20.)
    
    
def mag2db(x):
    """Convert magnitude to decibels (dB)

    The relationship between magnitude and decibels is:

    .. math::    X_{dB} = 20 * \log_{10}(x)
    
    .. doctest::
    
        >>> mag2db(0.1)
        -20.0
    
    .. seealso:: :func:`db2mag`
    """
    
    return 20. * numpy.log10(x)

