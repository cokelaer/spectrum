"""Criteria for parametric methods.

.. topic:: This module provides criteria to automatically select order in 
    parametric PSD estimate or pseudo spectrum estimates (e.g, music). 

    Some criteria such as the AIC criterion helps to chose the order of PSD 
    models such as the ARMA model. Nevertheless, it is difficult to estimate
    correctly the order of an ARMA model even by using these criteria. The 
    reason being that even the Akaike criteria (AIC) does not provide the 
    proper order with a probability of 1 with infinite samples. 

    The order choice is related to an expertise of the signal. There is no 
    exact criteria. However, they may provide useful information.

    AIC, AICc, KIC and AKICc are based on information theory.  They  attempt
    to balance the complexity (or length) of the model against how well the
    model fits the data.  AIC and KIC are biased estimates of the asymmetric
    and the symmetric Kullback-Leibler divergence respectively.  AICc and
    AKICc attempt to correct the bias. 
    
    There are also criteria related to eigen analysis, which takes as input
    the eigen values of any PSD estimate method.

.. rubric:: Example

.. plot::
    :width: 80%
    :include-source:

    from spectrum import *
    from pylab import *
    order = arange(1, 25)
    rho = [aryule(marple_data, i, norm='biased')[1] for i in order]
    plot(order, AIC(len(marple_data), rho, order), label='AIC')



:References: bd-Krim Seghouane and Maiza Bekara
   "A small sample model selection criterion based on Kullback's symmetric
   divergence", IEEE Transactions on Signal Processing,
   Vol. 52(12), pp 3314-3323, Dec. 2004

"""


class Criteria(object):
    """Criteria class for an automatic selection of ARMA order.

    Available criteria are

    ======= =====================
    ======= =====================
    AIC     see :func:`AIC`
    AICc    see :func:`AICc`
    KIC     see :func:`KIC`
    AKICc   see :func:`AKICc`
    FPE     see :func:`FPE`
    MDL     see :func:`MDL`
    CAT     see :func:`_CAT`
    ======= =====================


    """
    valid_criteria_names =  ['AIC', 'AICc', 'KIC', 'FPE', 'AKICc', 'MDL']
    error_incorrect_name = 'Invalid name provided. Correct names are %s ' \
        % valid_criteria_names
    error_no_criteria_found = 'No names match the valid criteria names (%s)' \
        % valid_criteria_names 
    def __init__(self, name, N):
        """Create a criteria object

        :param name: a string or list of strings containing valid criteria 
            method's name
        :param int N: size of the data sample.
            
        """
        #valid attributes
        self.__name = name
        self.__N = N
        self.__rho = 0
        self.__k = None
        self.__old_data = None
        self.__data = None
        self.__norm = True
        
    def _getName(self):
        return self.__name
    def _setName(self, name):
        assert isinstance(name, str), 'name must be a string'
        if name in self.valid_criteria_names:
            self.__name = name
        else:
            raise ValueError(self.error_no_criteria_found) 
    name = property(fget=_getName, fset=_setName, doc="Getter/Setter for the criteria name")

    def _getData(self):
        return self.__data
    def _setData(self, data):
        # save the data value in old_data is there is something to save
        if self.data == None:
            self.__data = data
            self.__old_data = 2.*data
        else:
            self.__old_data = self.data
            self.__data = data
    data = property(fget=_getData, fset=_setData, doc="Getter/Setter for the criteria output")
    
    def _getOldData(self):
        return self.__old_data
    old_data = property(fget=_getOldData, doc="Getter/Setter for the previous value")
    
    def _getK(self):
        return self.__k
    k = property(fget=_getK, doc="Getter for k the order of evaluation")
    
    def _getN(self):
        return self.__N
    def _setN(self, N):
        assert N>0, 'N must be positive'
        self.__N = N
    N = property(fget=_getN, fset=_setN, doc="Getter/Setter for N")
    
    def _getRho(self):
        return self.__rho
    def _setRho(self, rho):
        self.__rho = rho
    rho = property(fget=_getRho, fset=_setRho, doc="Getter/Setter for rho")


    def __call__(self, rho=None, k=None, N=None, norm=True):
        """Call the criteria function correspondign to :attr:`name`."""
        self.__norm = norm
        if N != None:
            self.N = N
        
        # we update rho only if it is needed (input different from self.rho)
        # if such case, we also update k
        if rho != None: 
            self.rho = rho
        if k != None:
            self.__k = k        
        self.__norm = norm
        #used to check if the criteria is reached or not
            
        f = eval(self.name)
        self.data = f(self.N, self.rho, self.k)
        # compare the new data with the previous one and return
        # False if the new value is larger so as to stop the iteration
        if self.old_data != None and self.data != None:
            if self.data > self.old_data:
                return False
            else:
                return True
        return True
        
    def plot(self):
        from pylab import plot
        plot(self.data)

    def plot_all(self):
        _ar_criteria(self.__rho, self.N)                    
        

def AIC(N, rho, k):
    r"""Akaike Information Criterion

    :param rho: rho at order k
    :param N: sample size
    :param k: AR order. 

    If k is the AR order and N the size of the sample, then Akaike criterion is

    .. math:: AIC(k) = \log(\rho_k) + 2\frac{k+1}{N}
    
    ::
    
        AIC(64, [0.5,0.3,0.2], [1,2,3])
    
    :validation: double checked versus octave.
    """
    from numpy import log, array
    #k+1 #todo check convention. agrees with octave
    res = N * log(array(rho)) + 2.* (array(k)+1)
    return res

def AICc(N, rho, k, norm=True):
    r"""corrected Akaike information criterion
    
    .. math:: AICc(k) = log(\rho_k) + 2 \frac{k+1}{N-k-2}
    
    
    :validation: double checked versus octave.
    """
    from numpy import log, array
    p = k  #todo check convention. agrees with octave 
    res = log(rho) + 2. * (p+1) / (N-p-2)
    return res
        
def KIC(N, rho, k):
    r"""Kullback information criterion
    
    .. math:: KIC(k) = log(\rho_k) + 3 \frac{k+1}{N}
    
    :validation: double checked versus octave.
    """
    from numpy import log, array
    res = log(rho) + 3. * (k+1.) /float(N)
    return res

def AKICc(N, rho, k):
    r"""approximate corrected Kullback information
    
    .. math:: AKICc(k) = log(rho_k) + \frac{p}{N*(N-k)} + (3-\frac{k+2}{N})*\frac{k+1}{N-k-2}
    
    """
    from numpy import log, array
    p = k
    res = log(rho) + +float(p)/N/(N-p) + (3.-(p+2.)/N) * (p+1.) / (N-p-2.)
    return res




def FPE(N,rho, k=None):
    r"""Final prediction error criterion

    .. math:: FPE(k) = \frac{N + k + 1}{N - k - 1} \rho_k

    :validation: double checked versus octave.

    """
    #k #todo check convention. agrees with octave
    fpe = rho * (N + k + 1.) / (N- k -1)
    return fpe

    

def MDL(N, rho, k):
    r"""Minimum Description Length

    .. math:: MDL(k) = N log \rho_k + p \log N

    :validation: results 
    """
    from numpy import log
    #p = arange(1, len(rho)+1)
    mdl = N* log(rho) + k * log(N)
    return mdl

def CAT(N, rho, k):
    r"""Criterion Autoregressive Transfer Function :

    .. math::  CAT(k) = \frac{1}{N} \sum_{i=1}^k \frac{1}{\rho_i} - \frac{\rho_i}{\rho_k} 
    
    .. todo:: validation
    """
    from numpy import zeros, arange
    cat = zeros(len(rho))
    for p in arange(1, len(rho)+1):
        rho_p = float(N)/(N-p)*rho[p-1]
        s = 0
        for j in range(1, p+1):
            rho_j = float(N)/(N-j)*rho[j-1]
            s = s + 1./rho_j
        print s, s/float(N), 1./rho_p
        cat[p-1] = s/float(N) - 1./rho_p
    return cat
    

def _ar_criteria(rho, N):
    pass
    #from pylab import plot, legend, hold, clf
    #fpe = FPE(N, rho)
    #aic = AIC(N, rho, )
    #mdl = MDL(rho, N)
    #cat = CAT(rho, N)
    #print fpe/max(fpe)
    #print aic/max(aic)
    #print mdl/max(mdl)
    #print cat/max(cat)
    #figure()
    #clf()
    #plot(fpe/fpe[0], 'o-', label='FPE')
    #hold(True)
    #plot(aic/aic[0], 'x-', label='AIC')
    #plot(mdl/mdl[0], '.-', label='MDL')
    #plot(cat/cat[-1], 's-', label='CAT')
    #legend()
    



def aic_eigen(s, N):
    r"""AIC order-selection using eigen values
    
    :param s: a list of `p` sorted eigen values
    :param N: the size of the input data. To be defined precisely. 
    
    :return: 
        * an array containing the AIC values
    
    Given :math:`n` sorted eigen values :math:`\lambda_i` with 
    :math:`0 <= i < n`, the proposed criterion from Wax and Kailath (1985) 
    is:
      
    .. math:: AIC(k) = -2(n-k)N \ln \frac{g(k)}{a(k)} + 2k(2n-k)
    
    where the arithmetic sum :math:`a(k)` is: 
    
    .. math:: a(k) = \sum_{i=k+1}^{n}\lambda_i
    
    and the geometric sum :math:`g(k)` is:
    
    .. math:: g(k) = \prod_{i=k+1}^{n} \lambda_i^{-(n-k)}
    
    The number of relevant sinusoids in the signal subspace is determined by 
    selecting the minimum of `AIC`.
    
    .. seealso:: :func:`~spectrum.eigenfreq.eigen`
    .. todo:: define precisely the input parameter N. Should be the input
       data length but when using correlation matrix (SVD), I suspect it 
       should be the length of the correlation matrix rather than the 
       original data.
       
    :References: 
        * [Marple]_ Chap 13, 
        * [Wax]_
    """
    import numpy as np
    
    kaic = []
    n = len(s)
    for k in range(0, n-1):         
        ak = 1./(n-k) * np.sum(s[k+1:])
        gk = np.prod(s[k+1:]**(1./(n-k))) 
        kaic.append( -2.*(n-k)*N * np.log(gk/ak) + 2.*k*(2.*n-k))
         
    return kaic

def mdl_eigen(s, N):
    r"""MDL order-selection using eigen values
    
    :param s: a list of `p` sorted eigen values
    :param N: the size of the input data. To be defined precisely. 
    
    :return: 
        * an array containing the AIC values
        
    .. math:: MDL(k) = (n-k)N \ln \frac{g(k)}{a(k)} + 0.5k(2n-k) log(N)
    
    .. seealso:: :func:`aic_eigen` for details
    
    :References: 
        * [Marple]_ Chap 13, 
        * [Wax]_
    """
    import numpy as np
    kmdl = []
    n = len(s)
    for k in range(0, n-1):         
        ak = 1./(n-k) * np.sum(s[k+1:])
        gk = np.prod(s[k+1:]**(1./(n-k))) 
        kmdl.append( -(n-k)*N * np.log(gk/ak) + 0.5*k*(2.*n-k)*np.log(N))
    return kmdl
        