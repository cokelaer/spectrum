"""
.. topic:: This module provides Periodograms (classics, daniell, bartlett)
        
        
    .. autosummary::
    
        Periodogram
        DaniellPeriodogram
        speriodogram
    
        
    .. codeauthor:: Thomas Cokelaer 2011
     
    :References: See [Marple]_
    
.. rubric:: Usage

You can compute a periodogram using :func:`speriodogram`:: 
    
    from spectrum import *
    from pylab import plot
    p = speriodogram(marple_data)
    plot(p)

However, the output is not always easy to manipulate or plot, therefore
it is advised to use the class :class:`Periodogram` instead::
    
    from spectrum import *
    p = Periodogram(marple_data)
    p()
    p.plot()
        
This class will take care of the plotting and internal state of 
the computation. For instance, if you can change the output easily::

    p.plot(sides='twosided')
    
"""
from window import Window
from psd import Spectrum, FourierSpectrum
from pylab import pi, fft, mean, rfft
from numpy import array, ceil
import numpy
import pylab as plt 
    
    
    
__all__ = ['pdaniell', 'speriodogram', 'Periodogram', 'pylab_periodogram', 
           'DaniellPeriodogram']

 
def speriodogram(x, NFFT=None, detrend=True, sampling=1., 
                   scale_by_freq=True, window='hamming', axis=0):
    """Simple periodogram, but matrices accepted.

    :param x: an array or matrix of data samples.
    :param NFFT: length of the data before FFT is computed (zero padding)
    :param bool detrend: detrend the data before co,puteing the FFT
    :param float sampling:
    :param scale_by_freq:
    :param str window:

    :return: 2-sided PSD if complex data, 1-sided if real.
    
    if a matrix is provided (using numpy.matrix), then a periodogram 
    is computed for each row. The returned matrix has the same shape as the input
    matrix.
    
    The mean of the input data is also removed from the data before computing
    the psd. 

    .. plot::
        :width: 80%
        :include-source:

        from pylab import *
        from spectrum import *
        data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)
        semilogy(speriodogram(data, detrend=False, sampling=1024), marker='o')
        grid(True)

    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        from pylab import *
        # create N data sets and make the frequency dependent on the time
        N = 100
        m = numpy.concatenate([data_cosine(N=1024, A=0.1, sampling=1024, freq=x) for x in range(1, N)]); 
        m.resize(N, 1024)
        res = speriodogram(m)
        figure(1)
        semilogy(res)
        figure(2)
        imshow(res.transpose(), aspect='auto')

    .. todo:: a proper spectrogram class/function 
    """
    x = array(x)
    # array with 1 dimension case
    if x.ndim == 1:
        print 'array case'
        axis = 0
        r = x.shape[0]
        w = Window(r, window)   #same size as input data
    # matrix case
    elif x.ndim == 2:
        print '2D array. each row is a 1D array'
        [r, c] = x.shape
        if r==1:
            r = c
            w = Window(x.size, window)   #same size as input data
        else:
            w = Window(c, window)   #same size as input data


    if NFFT == None:
        NFFT = len(x)
    
    isreal = numpy.isrealobj(x)
    
    if detrend == True:
        m = mean(x, axis=axis)
    else:
        m = 0
    if isreal == True:
        res =  (abs (rfft (x*w.data - m, NFFT))) ** 2. / r       
    else:
        res =  (abs (fft (x*w.data - m, NFFT ))) ** 2. / r
            
    if scale_by_freq is True:
        df = sampling / float(NFFT)
        res*= 2*pi/df

    return res.transpose()
    



def pylab_periodogram(data, NFFT=None,  sampling=1.):
    r"""Simple periodogram wrapper of numpy.psd with only nfft option

    This function is a simplified version of numpy.psd. The main difference
    is that you do not have the possibility to average over several segments, 
    as a Welch Periodogram would allow.

    :param A: the input data
    :param NFFT: NFFT padding
    :param window: the type of tapering


    :Technical documentation:

    When we calculate the periodogram of a set of data we get an estimation of the spectral density. In fact as we use a Fourier transform and a truncated segments the spectrum is the convolution of the data with a rectangular window which Fourier transform is

    .. math::

        W(s)= \frac{1}{N^2} \left[ \frac{\sin(\pi s)}{\sin(\pi s/N)} \right]^2

    Thus oscillations and sidelobes appears around the main frequency. One aim of t he tapering is to reduced this effects. We multiply data by a window whose  sidelobes are much smaller than the main lobe. Classical window is hanning window.  But other windows are available. However we must take into account this energy and divide the spectrum by energy of taper used. Thus periodogram becomes :

    .. math::

        D_k \equiv \sum_{j=0}^{N-1}c_jw_j \; e^{2\pi ijk/N}  \qquad k=0,...,N-1

    .. math::

        P(0)=P(f_0)=\frac{1}{2\pi W_{ss}}\arrowvert{D_0}\arrowvert^2

    .. math::

        P(f_k)=\frac{1}{2\pi W_{ss}} \left[\arrowvert{D_k}\arrowvert^2+\arrowvert{D_{N-k}}\arrowvert^2\right]        \qquad k=0,1,...,     \left( \frac{1}{2}-1 \right)

    .. math::

        P(f_c)=P(f_{N/2})= \frac{1}{2\pi W_{ss}} \arrowvert{D_{N/2}}\arrowvert^2

    with

    .. math::

        {W_{ss}} \equiv N\sum_{j=0}^{N-1}w_j^2


    .. plot::


        from spectrum import *
        psd = pylab_periodogram(marple_data, 256)


    """
    spectrum = Spectrum(data, sampling=1.)
    
    P = plt.psd(data, NFFT, sampling=sampling)
    spectrum.psd = P[0]
    #spectrum.__Spectrum_sides = 'twosided'
    
    return P, spectrum
    
    
class Periodogram(FourierSpectrum):
    """The Periodogram class provides an interface to periodogram PSDs
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import *
        data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)
        p = Periodogram(data, sampling=1024)
        p()
        p.plot(marker='o')
    
    
    """
    def __init__(self, data, sampling=1.,
                 window='hann', NFFT=None, scale_by_freq=False, 
                 detrend=None):
        """**Periodogram Constructor**
        
        :param array data:
        :param float sampling:
        :param str window:
        :param int NFFT:
        :param bool scale_by_freq:
        :param str detrend:
        
        """
        super(Periodogram, self).__init__(data, 
                                          window=window, 
                                          sampling=sampling,
                                          NFFT=NFFT, 
                                          scale_by_freq=scale_by_freq,
                                          detrend=detrend)

    def __call__(self):
        psd = speriodogram(self.data, window=self.window, sampling=self.sampling, 
                             NFFT=self.NFFT, scale_by_freq=self.scale_by_freq,
                             detrend=self.detrend)
        self.psd = psd
        if self.scale_by_freq is True:
            self.scale()
        
    def _str_title(self):
        return "Periodogram PSD estimate\n"
    
    def __str__(self):
        return super(Periodogram, self).__str__()
    
        
def DaniellPeriodogram(data, P, NFFT=None, detrend='mean', sampling=1., 
                       scale_by_freq=True, window='hamming'):
    r"""Return Daniell's periodogram.
    
    To reduce fast fluctuations of the spectrum one idea proposed by daniell 
    is to average each value with points in its neighboorhood. It's like 
    a low filter. 
    
    .. math:: \hat{P}_D[f_i]= \frac{1}{2P+1} \sum_{n=i-P}^{i+P} \tilde{P}_{xx}[f_n]  
    
    where P is the number of points to average.
    
    Daniell's periodogram is the convolution of the spectrum with a low filter:
    
    .. math:: \hat{P}_D(f)=   \hat{P}_{xx}(f)*H(f)
    
    Example::
    
        >>> DaniellPeriodogram(data,8)
        
    if N/P is not integer, the final values of the original PSD are not used.
    
    using DaniellPeriodogram(data, 0) should give the original PSD.
    
    """
    psd = speriodogram(data, NFFT=NFFT, detrend=detrend, sampling=sampling, 
                   scale_by_freq=scale_by_freq, window=window)
    
    if len(psd) %2 == 1:
        datatype = 'real'
    else:
        datatype = 'complex'
        
    N = len(psd)
    slice = 2 * P + 1
    if datatype== 'real': #must get odd value
        newN = ceil(psd.size/float(slice))
        if newN%2 == 0:
            newN = psd.size/slice
    else:
        newN = ceil(psd.size/float(slice))
        if newN%2 == 1:
            newN = psd.size/slice
    
    newpsd = numpy.zeros(newN) # keep integer division
    for i in range(0, newpsd.size):
        #print i,
        count = 0 #needed to know the number of valid averaged values
        for n in range(i*slice-P, i*slice+P+1): #+1 to have P values on each sides
            if n > 0 and n<N: #needed to start the average 
                count += 1
                newpsd[i] += psd[n]
        #print count
        newpsd[i] /= float(count)
    
    #todo: check this
    if datatype == 'complex':
        freq = plt.linspace(0,sampling, len(newpsd))
    else:
        df = 1./sampling
        freq = plt.linspace(0,sampling/2., len(newpsd))
    #psd.refreq(2*psd.size()/A.freq());  
    #psd.retime(-1./psd.freq()+1./A.size());

    return newpsd, freq




class pdaniell(FourierSpectrum):
    """The pdaniell class provides an interface to DaniellPeriodogram
    
    ::
    
    
        from spectrum import *
        data = data_cosine(N=4096, sampling=4096)
        p = pdaniell(data, 8, NFFT=4096)
        p()
        p.plot()
    
    
    """
    def __init__(self, data, P, sampling=1.,
                 window='hann', NFFT=None, scale_by_freq=True, 
                 detrend=None):
        """**pdaniell Constructor**
        
        :param array data:
        :param int P: number of neighbour to average over.
        :param float sampling:
        :param str window:
        :param int NFFT:
        :param bool scale_by_freq:
        :param str detrend:
        
        """
        super(pdaniell, self).__init__(data, 
                                       window=window, 
                                       sampling=sampling,
                                       NFFT=NFFT, 
                                       scale_by_freq=scale_by_freq,
                                       detrend=detrend)
        self.P = P
    def __call__(self):
        res = DaniellPeriodogram(self.data, self.P, window=self.window, sampling=self.sampling, 
                             NFFT=self.NFFT, scale_by_freq=self.scale_by_freq,
                             detrend=self.detrend)
        
        self.psd = res[0]
        
        
    def _str_title(self):
        return "Daniell Periodogram PSD estimate\n"
    
    def __str__(self):
        return super(pdaniell, self).__str__()
    
