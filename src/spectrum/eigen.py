import numpy
import toeplitz

__all__ = ["MINEIGVAL"]


def MINEIGVAL(T0, T, TOL):
    """Finds the minimum eigenvalue of a Hermitian Toeplitz matrix

    The classical power method is used together with a  fast  Toeplitz
    equation solution routine.  The eigenvector is normalized to unit length.

    :param T0: Scalar corresponding to real matrix element t(0)
    :param T: Array of  M  complex matrix elements t(1),...,t(M) C            from the left column of the Toeplitz matrix
    :param TOL: Real scalar tolerance; routine exits when [ EVAL(k) - EVAL(k-1) ]/EVAL(k-1) < TOL , where the index k denotes the iteration number.

    :return:
        * EVAL - Real scalar denoting the minimum eigenvalue of matrix
        * EVEC - Array of  M  complex eigenvector elements associated


    .. note::
     * External array T must be dimensioned >= M
     * array EVEC must be >= M+1
     *  Internal array E must be dimensioned >= M+1 .  

     * **dependencies**
        *  :meth:`spectrum.toeplitz.HERMTOEP`
    """
    M = len(T)
    eigval = 10
    eigvalold = 1
    eigvec  = numpy.zeros(M+1, dtype=complex)
    for k in range(0,M+1):
        eigvec[k] = 1+0j
    it=0
    #print 'initialisation',T0, T, eigval, eigvec
    maxit = 15
    while abs(eigvalold-eigval)>TOL*eigvalold and it<maxit:
        it=it+1
        eigvalold = eigval
        #print 'iteration ',it, 'eigvalold=',eigvalold, 'eigval=', eigval


        eig = toeplitz.HERMTOEP(T0, T, eigvec)
        SUM = 0
        save =0.+0j
        for k in range(0, M+1):
            SUM = SUM + eig[k].real**2+eig[k].imag**2
            save = save +eig[k]*eigvec[k].conjugate()
        SUM=1./SUM
        eigval = save.real*SUM
        for k in range(0,M+1):
            eigvec[k] = SUM * eig[k]
    if it==maxit:
        print 'warning reached max number of iteration (%s)' % maxit
    return eigval, eigvec


