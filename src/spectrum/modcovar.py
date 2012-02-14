import numpy
from numpy import real
from psd import ParametricSpectrum

__all__ = ['modcovar', 'modcovar_marple', 'pmodcovar']


def modcovar_marple (X,IP):
    """Fast algorithm for the solution of the modified covariance least squares normal equations.

    This implementation is based on [Marple]_. This code is far more 
    complicated and slower than :func:`modcovar` function, which is now the official version.
    See :func:`modcovar` for a detailed description of Modified Covariance method.

    :param X:    - Array of complex data samples X(1) through X(N)
    :param int IP:   - Order of linear prediction model (integer)

    :return:
        * P    - Real linear prediction variance at order IP
        * A    - Array of complex linear prediction coefficients
        * ISTAT - Integer status indicator at time of exit
              0. for normal exit (no numerical ill-conditioning)
              1. if P is not a positive value
              2. if DELTA' and GAMMA' do not lie in the range 0 to 1
              3. if P' is not a positive value
              4. if DELTA and GAMMA do not lie in the range 0 to 1

    
    :validation: the AR parameters are the same as those returned by 
        a completely different function :func:`modcovar`.
        
    .. note:: validation. results similar to test example in Marple but 
        starts to differ for    ip~8. with ratio of 0.975 for ip=15 probably 
        due to precision.

        
    :References: [Marple]_
    """
    Pv = [] 
    N = len(X)
    A = numpy.zeros(N, dtype=complex)
    D = numpy.zeros(N, dtype=complex)
    C = numpy.zeros(N, dtype=complex)
    R = numpy.zeros(N, dtype=complex)
    #   Initialization
    R1=0.
    for K in range(1, N-1):
        R1=R1 + 2.*(X[K].real**2 + X[K].imag**2)
    R2 = X[0].real**2 + X[0].imag**2
    R3 = X[N-1].real**2 + X[N-1].imag**2
    R4 = 1. / (R1 + 2. * (R2 + R3))
    P = R1 + R2 + R3
    DELTA = 1. - R2 * R4
    GAMMA = 1. - R3 * R4
    LAMBDA = (X[0] * X[N-1]).conjugate()*R4
    C[0] = X[N-1] * R4
    D[0] = X[0].conjugate() * R4
      
    M = 0
    if (IP ==0):
        P = (.5*R1+R2+R3)/float(N)
        return [], P, []

    #C   Main loop

    for M in range(0, IP):
        #print '---------------------', M
        SAVE1 = 0+0j
        for K in range(M+1, N):
            SAVE1 = SAVE1 + X[K]*X[K-M-1].conjugate()
        SAVE1 *= 2.
        R[M] = SAVE1.conjugate()
        THETA = X[N-1]*D[0]
        PSI=X[N-1]*C[0]
        XI = X[0].conjugate() * D[0]
        if M==0:
            pass
        else:
            for K in range(0, M):
                THETA=THETA+X[N-K-2]*D[K+1]    # Eq. (8.D.45)
                PSI = PSI + X[N-K-2]*C[K+1]    # Eq. (8.D.45)
                XI = XI + X[K+1].conjugate() * D[K+1] # Eq. (8.D.45)
                R[K] = R[K]-X[N-M-1] * X[N+1-M+K-1].conjugate() - X[M].conjugate() * X[M-K-1] # Eq. (8.D.37)
                SAVE1=SAVE1+R[K].conjugate()*A[M-K-1]             # Eq. (8.D.24)
                #print 'withini loop', K, THETA, PSI, XI, R[K], SAVE1
                #print 'x   -------', C[K], C[K+1]
                #print 'x   -------', D[K], D[K+1]
                #print 'x   -------', N-K-2, K+1,N-M-1,  M, M-K-1 
                #print 'x   -------', X[N-K-2], X[K+1],X[N-M-1],  X[M], X[M-K-1] 


        #  Order update of A vector

        C1 = -SAVE1/P
        A[M]=C1                                   # Eq. (8.D.23)
        P=P*(1.-C1.real**2-C1.imag**2)            # Eq. (8.D.25)
        #print 'C1=' , C1, 'A[M]', A[M], 'P=',P
        if M==0:
            pass
        else:
            #for k in range(0, (M-1)/2+1):
            #print 'BEFORE', A[0:4]
            for K in range(0, (M+1)/2):
                MK = M-K-1
                SAVE1=A[K]
                A[K]=SAVE1+C1*A[MK].conjugate() # Eq. (8.D.22)
                #print 'A uopdate---->',M, K, MK, A[K], A[MK], SAVE1, C1
                if (K != MK):
                    #print 'ENTERI IF', K, MK, A[MK], C1, SAVE1
                    A[MK]=A[MK]+C1*(SAVE1.conjugate()) # Eq. (8.D.22)
            #print 'AFTER', A[0:4]
        #print M+1, IP
        if M+1 == IP:
            #Pv.append(P)
            P=.5*P/float(N-M-1)
            Pv.append(P)
            return A, P, Pv
        else:
            Pv.append(.5*P/float(N-M-1))
            #Pv.append(P)

        #   Time update of C,D vectors and GAMMA,DELTA,LAMBDA scalars

        R1=1./(DELTA*GAMMA-(LAMBDA.real)**2-(LAMBDA.imag)**2)
        C1=(THETA*(LAMBDA.conjugate())+PSI*DELTA)*R1
        C2=(PSI*LAMBDA+THETA*GAMMA)*R1
        C3=(XI*(LAMBDA.conjugate())+THETA*DELTA)*R1
        C4=(THETA*LAMBDA+XI*GAMMA)*R1

        #print 'C1C2C3C4', C1, C2, C3, C4
        #print M, M/2+1, (M-1)/2+1
        for K in range(0, (M)/2+1):
            MK=M-K
            SAVE1=C[K].conjugate()
            SAVE2=D[K].conjugate()
            SAVE3=C[MK].conjugate()
            SAVE4=D[MK].conjugate()
            C[K]=C[K]+C1*SAVE3+C2*SAVE4                # Eq. (8.D.43)
            D[K]=D[K]+C3*SAVE3+C4*SAVE4                # Eq. (8.D.44)
            #print '-------->', K,MK, SAVE1, SAVE2, SAVE3, SAVE4
            if K != MK:
                C[MK]=C[MK]+C1*SAVE1+C2*SAVE2              # Eq. (8.D.43)
                D[MK]=D[MK]+C3*SAVE1+C4*SAVE2              # Eq. (8.D.44)
                #print 'witninloop loop ',K,MK,C[MK], D[MK]
        #print 'after loop'
        #print 'C=',C[0:2]
        #print 'D=',D[0:2]
        R2=PSI.real**2+PSI.imag**2
        R3=THETA.real**2+THETA.imag**2
        R4=XI.real**2+XI.imag**2
        R5=GAMMA-(R2*DELTA+R3*GAMMA+2.*real(PSI*LAMBDA*THETA.conjugate()))*R1
        R2=DELTA-(R3*DELTA+R4*GAMMA+2.*real(THETA*LAMBDA*XI.conjugate()))*R1
        GAMMA=R5                                     # Eq. (8.D.46)
        DELTA=R2                                     # Eq. (8.D.47)
        LAMBDA=LAMBDA+C3*PSI.conjugate()+C4*THETA.conjugate()  # Eq. (8.D.48)
        #print 'before time updaqt',R2, R3, R4, R5, LAMBDA
        if P <= 0.:
            raise ValueError('Found an invalid negative value for P.')
        if (DELTA > 0. and DELTA <= 1. and GAMMA > 0. and GAMMA <=1.):
            pass
        else:
            raise ValueError('Found an invalid DELTA or GAMMA value.')
        #   Time update of A vector; order updates of C,D vectors aand GAMMA,
        #   DELTA,LAMBDA scalars
        #print 'time update--------------'
        R1=1./P
        R2=1./(DELTA*GAMMA-LAMBDA.real**2-LAMBDA.imag**2) # Eq. (8.D.41)
        EF=X[M+1]
        EB=X[N-M-2]
        #print 'efeb', EF, EB
        #print 'A=',A[0:2]
        for K in range(0, M+1):
            #print 'efeb lloop', K,A[K], X[M-K], X[N-M+K-1]
            EF=EF+A[K]*X[M-K]                        # Eq. (8.D.1)
            EB=EB+A[K].conjugate()*X[N-M+K-1]            # Eq. (8.D.2)
        C1=EB*R1                                     # Eq. (8.D.28)
        C2=EF.conjugate()*R1                              # Eq. (8.D.29)
        C3=(EB.conjugate()*DELTA+EF*LAMBDA)*R2
        C4=(EF*GAMMA+(EB*LAMBDA).conjugate())*R2

        #print 'c1c2c3c4', C1 ,C2, C3, C4
        #print 'efeb',EF, EB

        for K in range(M, -1, -1):
            SAVE1=A[K]
            A[K]=SAVE1+C3*C[K]+C4*D[K]                 # Eq. (8.D.38)
            C[K+1]=C[K]+C1*SAVE1                       # Eq. (8.D.26)
            D[K+1]=D[K]+C2*SAVE1                       # Eq. (8.D.27)
            #print 'ACD',K, A[K], C[K+1], D[K+1]
        C[0]=C1
        D[0]=C2
        R3=EB.real**2+EB.imag**2
        R4=EF.real**2+EF.imag**2
        P=P-(R3*DELTA+R4*GAMMA+2.*real(EF*EB*LAMBDA))*R2  # Eq. (8.D.42)
        DELTA=DELTA-R4*R1                            # Eq. (8.D.32)
        GAMMA=GAMMA-R3*R1                            # Eq. (8.D.33)
        LAMBDA=LAMBDA+(EF*EB).conjugate()*R1                # Eq. (8.D.35)

        #print 'final'
        #print R3, R4, P, DELTA, GAMMA, LAMBDA
        #print 'Afinal=',A[0:2]
        #print 'P=',P
        if (P > 0.): 
            pass
        else:
            raise ValueError("Found a invalid negative P value ")
        if (DELTA > 0. and DELTA <= 1. and  GAMMA > 0. and GAMMA <= 1.):
            pass
        else:
            raise ValueError("Found an invalid negative GAMMA or DELTA value ")



 


def modcovar(x, order):
    """Simple and fast implementation of the covariance AR estimate
    
    This code is 10 times faster than :func:`modcovar_marple` and more importantly
    only 10 lines of code, compared to a 200 loc for :func:`modcovar_marple`

    :param X:        Array of complex data samples
    :param int order:   Order of linear prediction model

    :return:
        * P    - Real linear prediction variance at order IP
        * A    - Array of complex linear prediction coefficients
    

    .. plot::
        :include-source:
        :width: 80%
        
        from spectrum import *
        from pylab import *

        a, p = modcovar(marple_data, 15)
        PSD = arma2psd(a)
        PSD = cshift(PSD, len(PSD)/2) # switch positive and negative freq
        plot(linspace(-0.5, 0.5, 4096), 10*log10(PSD/max(PSD)))
        axis([-0.5,0.5,-60,0])
    
    .. seealso:: :class:`~spectrum.modcovar.pmodcovar`
    
    :validation: the AR parameters are the same as those returned by 
        a completely different function :func:`modcovar_marple`.
    
    
    :References: Mathworks 
    """
    from spectrum import corrmtx
    import scipy.linalg
    X = corrmtx(x, order, 'modified')
    Xc = numpy.matrix(X[:,1:]) 
    X1 = numpy.array(X[:,0])
    
    # Coefficients estimated via the covariance method
    # Here we use lstsq rathre than solve function because Xc is not square matrix
    a, residues, rank, singular_values = scipy.linalg.lstsq(-Xc, X1)
    
    # Estimate the input white noise variance
    

    Cz = numpy.dot(X1.conj().transpose(), Xc)
    e = numpy.dot(X1.conj().transpose(), X1) + numpy.dot(Cz, a)
    assert e.imag < 1e-4, 'wierd behaviour'
    e = float(e.real) # ignore imag part that should be small
    
    return a, e




class pmodcovar(ParametricSpectrum):
    """Class to create PSD based on modified covariance algorithm 
    
    See :func:`modcovar` for description.

    .. rubric:: Examples
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        p = pmodcovar(marple_data, 15, NFFT=4096)
        p()
        p.plot(sides='centerdc')

    .. seealso:: :class:`modcovar`

    """
    def __init__(self, data, order, NFFT=None, sampling=1., 
                 scale_by_freq=False):
        """**Constructor**

        For a detailled description of the parameters, see :func:`modcovar`.

        :param array data:     input data (list or numpy.array)
        :param int order:  
        :param int NFFT:       total length of the final data sets (padded with zero if needed; default is 4096)
        :param float sampling: sampling frequency of the input :attr:`data`.


        """
        super(pmodcovar, self).__init__(data, ar_order=order, 
                                            NFFT=NFFT, sampling=sampling, 
                                            scale_by_freq=scale_by_freq)

    def __call__(self):
        from spectrum import arma2psd
        ar, e = modcovar(self.data, self.ar_order)
        self.ar = ar
        psd = arma2psd(A=ar, T=self.sampling, NFFT=self.NFFT)
        
        if self.datatype == 'real':
            from tools import twosided_2_onesided
            newpsd  = twosided_2_onesided(psd)
            #\psd[0:self.NFFT/2]*2
            #newpsd[0] /= 2.
            #newpsd = numpy.append(newpsd, psd[-1])
            self.psd = newpsd
        else:
            self.psd = psd
        if self.scale_by_freq is True:
            self.scale()
        
    def _str_title(self):
        return "Modified covariance PSd estimate\n"
    
    def __str__(self):
        return super(pmodcovar, self).__str__()
