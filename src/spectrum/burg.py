"""BURG method of AR model estimate

.. topic:: This module provides BURG method and BURG PSD estimate

    .. autosummary::

        arburg
        pburg

    .. codeauthor:: Thomas Cokelaer 2011

"""
#TODO: convert arburg into arburg2 to get a nicer and faster algorithm.
import logging

import numpy as np
from spectrum.psd import ParametricSpectrum

__all__ = ["arburg", 'pburg']


def _arburg2(X, order):
    """This version is 10 times faster than arburg, but the output rho is not correct.


    returns [1 a0,a1, an-1]

    """
    x = np.array(X)
    N = len(x)

    if order <= 0.:
        raise ValueError("order must be > 0")

    # Initialisation
    # ------ rho, den
    rho = sum(abs(x)**2.) / N  # Eq 8.21 [Marple]_
    den = rho * 2. * N

    # ------ backward and forward errors
    ef = np.zeros(N, dtype=complex)
    eb = np.zeros(N, dtype=complex)
    for j in range(0, N):  #eq 8.11
        ef[j] = x[j]
        eb[j] = x[j]

    # AR order to be stored
    a = np.zeros(1, dtype=complex)
    a[0] = 1
    # ---- rflection coeff to be stored
    ref = np.zeros(order, dtype=complex)

    temp = 1.
    E = np.zeros(order+1)
    E[0] = rho

    for m in range(0, order):
        #print m
        # Calculate the next order reflection (parcor) coefficient
        efp = ef[1:]
        ebp = eb[0:-1]
        #print efp, ebp
        num = -2.* np.dot(ebp.conj().transpose(),  efp)
        den = np.dot(efp.conj().transpose(),  efp)
        den += np.dot(ebp,  ebp.conj().transpose())
        ref[m] = num / den

        # Update the forward and backward prediction errors
        ef = efp + ref[m] * ebp
        eb = ebp + ref[m].conj().transpose() * efp

        # Update the AR coeff.
        a.resize(len(a)+1, refcheck=False)
        a = a + ref[m] * np.flipud(a).conjugate()

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
        p.plot(sides='centerdc')

    Another example based on a real data set is shown here below. Note here
    that we set the scale_by_freq value to False and True. False should give
    results equivalent to octave or matlab convention while setting to True
    implies that the data is multiplied by  :math:`2\pi df` where 
    :math:`df = sampling / N`.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import data_two_freqs, pburg
        p = pburg(data_two_freqs(), 7, NFFT=4096)
        p.plot()
        p = pburg(data_two_freqs(), 7, NFFT=4096, scale_by_freq=True)
        p.plot()
        from pylab import legend
        legend(["un-scaled", "scaled by 2 pi df"])

    """
    def __init__(self, data, order, criteria=None, NFFT=None, sampling=1.,
        scale_by_freq=False):
        """**Constructor**

        For a detailled description of the parameters, see :func:`burg`.

        :param array data:     input data (list or np.array)
        :param int order:
        :param str criteria:
        :param int NFFT: total length of the final data sets (padded with zero if 
            needed; default is 4096)

        :param float sampling: sampling frequency of the input :attr:`data`.

        """
        super(pburg, self).__init__(data, ar_order=order,
                sampling=sampling, NFFT=NFFT,
                scale_by_freq=scale_by_freq)
        self.criteria = criteria

    def __call__(self):
        from .arma import arma2psd
        ar, rho, ref = arburg(self.data, self.ar_order, self.criteria)
        self.ar = ar
        self.rho = rho
        self.reflection = ref
        psd = arma2psd(A=self.ar, B=self.ma, rho=self.rho,
                      T=self.sampling, NFFT=self.NFFT)

        if self.datatype == 'real':
            #  see doc/concepts.rst for details
            if self.NFFT % 2 == 0:
                newpsd  = psd[0:int(self.NFFT/2 +1)] * 2
            else:
                newpsd  = psd[0:int((self.NFFT+1)/2)] * 2
            self.psd = newpsd
        else:
            self.psd = psd
        self.scale()
        return self

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
    if order <= 0.:
        raise ValueError("order must be > 0")

    if order > len(X):
        raise ValueError("order must be less than length input - 2")

    x = np.array(X)
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
        logging.debug('Step {}. old criteria={} new one={}.  rho={}'.format(
                0, crit.old_data, crit.data, rho))

    #p =0
    a = np.zeros(0, dtype=complex)
    ref = np.zeros(0, dtype=complex)
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
            logging.debug('Step {}. old criteria={} new one={}. rho={}'.format(
                k+1, crit.old_data, crit.data, new_rho))
            #k+1 because order goes from 1 to P whereas k starts at 0.
            status = crit(rho=temp*rho, k=k+1)
            if status is False:
                logging.debug('Stop criteria reached %s %s ' % (crit.data, crit.old_data))
                break
        # this should be after the criteria
        rho = new_rho
        if rho <= 0:
            raise ValueError("Found a negative value (expected positive strictly) %s. Decrease the order" % rho)

        a.resize(a.size+1,refcheck=False)
        a[k] = kp
        if k == 0:
            for j in range(N-1, k, -1):
                save2 = ef[j]
                ef[j] = save2 + kp * eb[j-1]          # Eq. (8.7)
                eb[j] = eb[j-1] + kp.conjugate() *  save2

        else:
            # update the AR coeff
            khalf = (k+1)//2  # FIXME here khalf must be an integer
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
        ref.resize(ref.size+1, refcheck=False)
        ref[k] = kp

    return a, rho, ref
