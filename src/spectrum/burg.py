"""BURG method of AR model estimate

.. topic:: This module provides BURG method and BURG PSD estimate 

    .. autosummary::

        arburg
        pburg

    .. codeauthor:: Thomas Cokelaer 2011

"""


#TODO: convert arburg into arburg2 to get a nicer and faster algorithm.
 
import numpy
from psd import ParametricSpectrum

__all__ = ["arburg", 'pburg']


def _arburg2(X, order):
    """This version is 10 times faster than arburg, but the output rho is not correct.
    
    
    returns [1 a0,a1, an-1]
    
    """
    x = numpy.array(X)
    N = len(x)

    if order == 0.: 
        raise ValueError("order must be > 0")

    # Initialisation
    # ------ rho, den
    rho = sum(abs(x)**2.) / N  # Eq 8.21 [Marple]_
    den = rho * 2. * N 

    # ------ backward and forward errors
    ef = numpy.zeros(N, dtype=complex)
    eb = numpy.zeros(N, dtype=complex)    
    for j in range(0, N):  #eq 8.11
        ef[j] = x[j]
        eb[j] = x[j]

    # AR order to be stored
    a = numpy.zeros(1, dtype=complex)
    a[0] = 1
    # ---- rflection coeff to be stored
    ref = numpy.zeros(order, dtype=complex)

    temp = 1.
    E = numpy.zeros(order+1)
    E[0] = rho

    for m in range(0, order):
        #print m
        # Calculate the next order reflection (parcor) coefficient
        efp = ef[1:]
        ebp = eb[0:-1]
        #print efp, ebp
        num = -2.* numpy.dot(ebp.conj().transpose(),  efp)
        den = numpy.dot(efp.conj().transpose(),  efp)
        den += numpy.dot(ebp,  ebp.conj().transpose())
        ref[m] = num / den

        # Update the forward and backward prediction errors
        ef = efp + ref[m] * ebp
        eb = ebp + ref[m].conj().transpose() * efp

        # Update the AR coeff.
        a.resize(len(a)+1)
        a = a + ref[m] * numpy.flipud(a).conjugate()

        # Update the prediction error
        E[m+1] = (1 - ref[m].conj().transpose()*ref[m]) * E[m]
        #print 'REF', ref, num, den
    return a, E[-1], ref




class pburg(ParametricSpectrum):
    """Class to create PSD based on Burg algorithm 
    
    See :func:`arburg` for description.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        p = pburg(marple_data, 15, NFFT=4096)
        p()
        p.plot(sides='centerdc')


    """
    def __init__(self, data, order, criteria=None, NFFT=None, sampling=1.):
        """**Constructor**

        For a detailled description of the parameters, see :func:`burg`.

        :param array data:
        :param int order:  
        :param str criteria:
        :param int NFFT:
        :param float sampling: 

        """
        super(pburg, self).__init__(data, ar_order=order, 
                                            NFFT=NFFT, sampling=sampling)
        self.criteria = criteria

    def __call__(self):
        from arma import arma2psd
        ar, rho, ref = arburg(self.data, self.ar_order, self.criteria)
        self.ar = ar
        self.rho = rho
        self.reflection = ref
        psd = arma2psd(A=self.ar, B=self.ma, rho=self.rho, 
                      T=self.sampling, NPSD=self.NFFT)
        #self.psd = psd
        if self.datatype == 'real':
            newpsd  = psd[0:self.NFFT/2]*2
            newpsd[0] /= 2.
            newpsd = numpy.append(newpsd, psd[-1])
            self.psd = newpsd
        else:
            self.psd = psd
        self.scale()
        
    def _str_title(self):
        return "Periodogram PSD estimate\n"
    
    def __str__(self):
        return super(pburg, self).__str__()




def arburg(X, order, criteria=None):
    r"""Estimate the complex autoregressive parameters by the Burg algorithm.

    .. math:: x(n) = \sqrt{(v}) e(n) + \sum_{k=1}^{P+1} a(k) x(n-k)

    :param x:  Array of complex data samples (length N)
    :param order: Order of autoregressive process (0<order<N)
    :param criteria: select a criteria to automatically select the order 

    :return:
        * A Array of complex autoregressive parameters A(1) to A(order). First 
          value (unity) is not included !! 
        * P Real variable representing driving noise variance (mean square 
          of residual noise) from the whitening operation of the Burg 
          filter.
        * reflection coefficients defining the filter of the model.

    .. plot::
        :width: 80%
        :include-source:

        from pylab import plot, log10, linspace, axis
        from spectrum import *

        AR, P, k = arburg(marple_data, 15)
        PSD = arma2psd(AR, sides='centerdc')
        plot(linspace(-0.5, 0.5, len(PSD)), 10*log10(PSD/max(PSD)))
        axis([-0.5,0.5,-60,0])


    .. note::
        1. no detrend. Should remove the mean trend to get PSD. Be careful if
           presence of large mean.
        2. If you don't know what the order value should be, choose the
           criterion='AKICc', which has the least bias and best 
           resolution of model-selection criteria.
      
    .. note:: real and complex results double-checked versus octave using 
        complex 64 samples stored in marple_data. It does not agree with Marple
        fortran routine but this is due to the simplex precision of complex
        data in fortran.
      
    :reference: [Marple]_ [octave]_   
    """
    if order == 0.: 
        raise ValueError("order must be > 0")

    x = numpy.array(X)
    N = len(x)

    # Initialisation
    # ------ rho, den
    rho = sum(abs(x)**2.) / float(N)  # Eq 8.21 [Marple]_
    den = rho * 2. * N 

    # ---- criteria
    if criteria:
        from spectrum import Criteria
        crit = Criteria(name=criteria, N=N)
        crit.data = rho
        print 0, 'old criteria=',crit.old_data, 'new criteria=',crit.data, 'new_rho=', rho

    #p =0
    a = numpy.zeros(0, dtype=complex)
    ref = numpy.zeros(0, dtype=complex)
    ef = x.astype(complex)
    eb = x.astype(complex)
    temp = 1.
    #   Main recursion
    
    for k in range(0, order):
        
        # calculate the next order reflection coefficient Eq 8.14 Marple
        num = sum([ef[j]*eb[j-1].conjugate() for j in range(k+1, N)])
        den = temp * den - abs(ef[k])**2 - abs(eb[N-1])**2  
        kp = -2. * num / den #eq 8.14
        
        temp = 1. - abs(kp)**2.
        new_rho = temp * rho 

        if criteria:
            print k+1, 'old criteria=',crit.old_data, 'new criteria=',crit.data, 'new_rho=',new_rho
            #k+1 because order goes from 1 to P whereas k starts at 0.
            status = crit(rho=temp*rho, k=k+1)
            if status is False:
                #print 'should stop here-----------------', crit.data, crit.old_data 
                break
        # this should be after the criteria
        rho = new_rho
        if rho <= 0:
            raise ValueError("Found an negative value (expected positive stricly) %s" % rho)

        a.resize(a.size+1)
        a[k] = kp
        if k == 0:
            for j in range(N-1, k, -1):     
                save2 = ef[j]
                ef[j] = save2 + kp * eb[j-1]          # Eq. (8.7)
                eb[j] = eb[j-1] + kp.conjugate() *  save2
            
        else:
            # update the AR coeff
            khalf = (k+1)/2
            for j in range(0, khalf):   
                ap = a[j] # previous value
                a[j] = ap + kp * a[k-j-1].conjugate()      # Eq. (8.2)     
                if j != k-j-1:
                    a[k-j-1] = a[k-j-1] + kp * ap.conjugate()    # Eq. (8.2)

            # update the prediction error
            for j in range(N-1, k, -1):
                save2 = ef[j]
                ef[j] = save2 + kp * eb[j-1]          # Eq. (8.7)
                eb[j] = eb[j-1] + kp.conjugate() *  save2
            
        # save the reflection coefficient
        ref.resize(ref.size+1)        
        ref[k] = kp
    return a, rho, ref
