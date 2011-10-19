"""
.. topic:: Cholesky methods

    .. autosummary::

        CHOLESKY
        
    .. codeauthor:: Thomas Cokelaer, 2011
"""
import numpy

__all__ = ["CHOLESKY"]

def _numpy_cholesky(A, B):
    """Solve Ax=B using numpy cholesky solver

    A = LU

    in the case where A is square and Hermitian, A = L.L* where L* is
    transpoed and conjugate matrix

    Ly = b

    where

    Ux=y

    so x = U^{-1} y
    where U = L*
    and y = L^{-1} B
    """
    L = numpy.linalg.cholesky(A)
    # A=L*numpy.transpose(L).conjugate()
    # Ly = b
    y = numpy.linalg.solve(L,B)
    # Ux = y
    x = numpy.linalg.solve(L.transpose().conjugate(),y)
    return x, L

def _numpy_solver(A, B):
    """This function solve Ax=B directly without taking care of the input
    matrix properties.
    """
    x = numpy.linalg.solve(A, B)
    return x

def CHOLESKY(A, B, method='scipy'):
    """Solve linear system `AX=B` using CHOLESKY method.
    
    :param A: an input Hermitian matrix
    :param B: an array
    :param str method: a choice of method in [numpy, scipy, numpy_solver]
        
        * `numpy_solver` relies entirely on numpy.solver (no cholesky decomposition)
        * `numpy` relies on the numpy.linalg.cholesky for the decomposition and
          numpy.linalg.solve for the inversion.
        * `scipy` uses scipy.linalg.cholesky for the decomposition and 
          scipy.linalg.cho_solve for the inversion.
    
    .. rubric:: Description
    
    When a matrix is square and Hermitian (symmetric with lower part being 
    the complex conjugate of the upper one), then the usual triangular
    factorization takes on the special form:
    
    .. math:: A = R R^H 
    
    where :math:`R` is a lower triangular matrix with nonzero real principal
    diagonal element. The input matrix can be made of complex data. Then, the 
    inversion to find :math:`x` is made as follows: 
    
    .. math::  Ry = B
    
    and 
    
    .. math::   Rx = y
    
    .. doctest::

        >>> import numpy
        >>> from spectrum import CHOLESKY
        >>> A = numpy.array([[ 2.0+0.j ,  0.5-0.5j, -0.2+0.1j],
        ...    [ 0.5+0.5j,  1.0+0.j ,  0.3-0.2j],
        ...    [-0.2-0.1j,  0.3+0.2j,  0.5+0.j ]])
        >>> B = numpy.array([ 1.0+3.j ,  2.0-1.j ,  0.5+0.8j])
        >>> CHOLESKY(A, B)
        array([ 0.95945946+5.25675676j,  4.41891892-7.04054054j,
               -5.13513514+6.35135135j])

    """
    if method == 'numpy_solver':
        X = _numpy_solver(A,B)
        return X
    elif method == 'numpy':
        X, _L = _numpy_cholesky(A, B)
        return X
    elif method == 'scipy':
        import scipy.linalg
        L = scipy.linalg.cholesky(A)
        X = scipy.linalg.cho_solve((L, False), B)
    else:
        raise ValueError('method must be numpy_solver, numpy_cholesky or cholesky_inplace')
    return X
