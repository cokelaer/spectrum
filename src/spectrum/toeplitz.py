"""
.. topic:: Toeplitz module

    These functions are not yet used by other functions, which explains the
    lack of documentation, test, examples. 
    
    Nevertheless, they can be used for production. 

    .. autosummary::

        HERMTOEP
        
    .. codeauthor:: Thomas Cokelaer, 2011
"""

import numpy

__all__ = [ "HERMTOEP"]


def TOEPLITZ(T0, TC, TR, Z):
    """solve the general toeplitz linear equations

    Solve TX=Z
    
    :param T0: zero lag value
    :param TC: r1 to rN 
    :param TR: r1 to rN

    returns X

    requires 3M^2+M operations instead of M^3 with gaussian elimination
    
    .. warning:: not used right now
    """
    assert len(TC)>0
    assert len(TC)==len(TR)
    M = len(TC)
    X = numpy.zeros(M+1,dtype=complex)
    A = numpy.zeros(M,dtype=complex)
    B = numpy.zeros(M,dtype=complex)
    P = T0
    if P == 0: raise ValueError("P must be different from zero")
    if P == 0: raise ValueError("P must be different from zero")
    X[0] = Z[0]/T0 
    for k in range(0, M):
        save1 = TC[k]
        save2 = TR[k]
        beta = X[0]*TC[k]
        if k == 0:
            temp1 = -save1 / P
            temp2 = -save2 / P
        else:
            for j in range(0, k):
                save1 = save1 + A[j] * TC[k-j-1]
                save2 = save2 + B[j] * TR[k-j-1]
                beta = beta + X[j+1] * TC[k-j-1]
            temp1 = -save1 / P
            temp2 = -save2/P
        P = P * (1. - (temp1*temp2))
        if P <= 0:
            raise ValueError("singular matrix")
        A[k] = temp1
        B[k] = temp2
        alpha = (Z[k+1]-beta)/P
        if k == 0: 
            X[k+1] = alpha
            for j in range(0,k+1):
                X[j] = X[j] + alpha * B[k-j]
            continue
        
        for j in range(0, k):
            kj = k-j-1
            save1 = A[j]
            A[j] = save1 + temp1 * B[kj] 
            B[kj] = B[kj] + temp2*save1
            
        X[k+1] = alpha
        for j in range(0,k+1):
            X[j] = X[j] + alpha*B[k-j]
    return X


def HERMTOEP(T0, T, Z):
    """solve Tx=Z by a variation of Levinson algorithm where T 
    is a complex hermitian toeplitz matrix

    :param T0: zero lag value
    :param T: r1 to rN
     
    :return: X
    
    used by eigen PSD method
    """
    assert len(T)>0
    M = len(T)
    X = numpy.zeros(M+1,dtype=complex)
    A = numpy.zeros(M,dtype=complex)
    P = T0
    if P == 0: raise ValueError("P must be different from zero")
    X[0] = Z[0]/T0 
    for k in range(0, M):
        save = T[k]
        beta = X[0]*T[k]
        if k == 0:
            temp = -save / P
        else:
            for j in range(0, k):
                save = save + A[j] * T[k-j-1]
                beta = beta + X[j+1] * T[k-j-1]
            temp = -save / P
        P = P * (1. - (temp.real**2+temp.imag**2))
        if P <= 0:
            raise ValueError("singular matrix")
        A[k] = temp
        alpha = (Z[k+1]-beta)/P

        if k == 0:
            #print 'skipping code for k=0' 
            X[k+1] = alpha
            for j in range(0,k+1):
                X[j] = X[j] + alpha * A[k-j].conjugate()
            continue
        khalf = (k+1)/2
        for j in range(0, khalf):
            kj = k-j-1
            save=A[j]
            A[j] = save+temp*A[kj].conjugate() 
            if j != kj:
                A[kj] = A[kj] + temp*save.conjugate()
        X[k+1] = alpha
        for j in range(0,k+1):
            X[j] = X[j] + alpha * A[k-j].conjugate()
    return X


