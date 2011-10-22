"""
.. topic:: This module contains tapering windows utilities.

    .. autosummary::

        Window
        window_visu
        create_window
        window_hann
        window_hamming
        window_bartlett
        ...
        
    .. codeauthor:: Thomas Cokelaer 2011
 
    :References: See [Nuttall]_, [Marple]_, [Harris]_
"""
from numpy import pi, cos, arange, array, sin, exp, sinc, linspace, \
    sqrt, ones, sum

window_names = {
                 'bartlett_hann': 'window_bartlett_hann',
                 'blackman_harris': 'window_blackman_harris',
                 'blackman_nuttall': 'window_blackman_nuttall', 
                 'bohman':'window_bohman',                
                 'blackman': 'window_blackman',
                 'chebwin':'window_chebwin',
                 'gaussian': 'window_gaussian',
                 'hamming':'window_hamming',
                 'kaiser': 'window_kaiser',
                 'lanczos': 'window_lanczos',
                 'sinc': 'window_lanczos',
                 'poisson': 'window_poisson',
                 'tukey':'window_tukey',
                 'nuttall':'window_nuttall',
                 'parzen':'window_parzen',
                 'flattop': 'window_flattop',
                  'riesz':'window_riesz',
                 'riemann':'window_riemann',
                 'hann':'window_hann',
                 'hanning':'window_hann',
                 'poisson_hanning': 'window_poisson_hanning',
                 'rectangular':'window_rectangle',
                 'rectangle': 'window_rectangle',
                 'bartlett':'window_bartlett',
                 'triangular':     'window_bartlett',
                 'cosine': 'window_cosine', 
                 'sine':   'window_cosine',
                 'cauchy':'window_cauchy',
                 }

class Window(object):
    r"""Window tapering object

    This class provides utilities to manipulate tapering windows. 
    Plotting functions allows to visualise the time and frequency response.
    It is also possible to retrieve relevant quantities such as the
    equivalent noise band width. 

    The following examples illustrates the usage. First, we create the window
    by providing a name and a size:: 
    
        from spectrum import *
        w = Window(64, 'hamming')
        
    The window has been computed and the data is stored in::
    
        w.data
    
    This object contains plotting methods so that you can see the time or 
    frequency response. 
    
    
    .. plot::
        :include-source:
        :width: 80%

        from spectrum.window import Window
        w = Window(64, 'hamming')
        w.plot_frequencies()

    Some windows may accept optional arguments. For instance, :func:`window_blackman` 
    accepts an optional argument called :math:`\alpha` as well as :class:`Window`.
    Indeed, we use the factory :func:`~spectrum.window.create_window`, which 
    manage all the optional arguments. So you can write::
    
        w = Window(64, 'blackman', alpha=1)
    
    .. seealso:: :func:`create_window`.
    
    
    .. rubric:: Constructor:
    """
    def __init__(self, N, name=None, norm=True, **kargs):
        """Create a tapering window object

        :param N: the window length
        :param name: the type of window, e.g., 'Hann'
        :param norm: normalise the window in frequency domain (for plotting)
        :param kargs: any of :func:`create_window` valid optional arguments.
        
        
        .. rubric:: Attributes:
        
        * data: time series data
        * frequencies: getter to the frequency series
        * response: getter to the PSD
        * enbw:  getter to the Equivalent noise band width.
        
        """
        #input attributse
        assert N>0 , "First argument  must be positive"
        if name == None or name not in window_names.keys():
            raise ValueError("second argument must be a valid window name %s" %
                             window_names.keys())
        
        self.__N = N
        self.__name = name
        self.__norm = norm

        #compute the window data and Fourier domain data
        self.__data = create_window(N, name, **kargs)
        self.__frequencies = None
        self.__response = None
        # equivalent noise BW and other quantities.
        self.__enbw = enbw(self.data)

    def _getNorm(self):
        return self.__norm
    norm = property(fget=_getNorm, doc="Getter of the normalisation flag (True by default)")

    def _getN(self):
        return self.__N
    N = property(fget=_getN, doc="Getter for the window length")

    def _getENBW(self):
        return self.__enbw
    enbw = property(fget=_getENBW, doc="getter for the equivalent noise band\
        width. See :func:`enbw` function")

    def _getMeanSquare(self):
        return sum(self.data**2)/self.N
    mean_square = property(fget=_getMeanSquare, doc="returns :math:`\frac{w^2}{N}`")


    def _getName(self):
        return self.__name
    name = property(fget=_getName, doc="Getter for the window name")

    def _getData(self):
        return self.__data
    data = property(fget=_getData, doc="Getter for the window values (in time)")

    def _getF(self):
        if self.__response == None:
            self.compute_response()
        self.__frequencies = linspace(-0.5, 0.5, len(self.__response))
        return self.__frequencies
    frequencies = property(fget=_getF, doc="Getter for the frequency array")

    def _getResponse(self):
        if self.__response == None:
            self.compute_response()
        return self.__response
    response = property(fget=_getResponse, doc="Getter for the frequency \
        response. See :meth:`compute_response`")


    def compute_response(self, **kargs):
        """Compute the window data frequency response

        :param norm: True by default. normalised the frequency data.
        :param NFFT: 2048 by default. if less than data length, then 
            NFFT is set to the data length*2.
            
        The response is stored in :attr:`response`.
        
        .. note:: Units are dB (20 log10) since we plot the frequency response)
        
        
        """
        from pylab import fft, fftshift, log10

        norm = kargs.get('norm', self.norm)
        
        # do some padding. Default is max(2048, data.len*2) 
        NFFT = kargs.get('NFFT', 2048)
        if NFFT < len(self.data): 
            NFFT = self.data.size * 2

        # compute the fft modulus
        A = fft(self.data, NFFT)
        mag = abs(fftshift(A))
        
        # do we want to normalise the data
        if norm is True:
            mag = mag / max(mag)
        response = 20. * log10(mag) # factor 20 we are looking at the response
                                    # not the powe
        #response = clip(response,mindB,100)
        self.__response = response

    def plot_frequencies(self, mindB=None, maxdB=None, norm=True):
        """Plot the window in the frequency domain

        :param mindB: change the default lower y bound
        :param maxdB: change the default upper lower bound 
        :param bool norm: if True, normalise the frequency response.
            
        .. plot::
            :width: 80%
            :include-source:

            from spectrum.window import Window
            w = Window(64, name='hamming')
            w.plot_frequencies()

        """
        from pylab import plot, title, xlim, grid, ylim, xlabel, ylabel
        # recompute the response
        self.compute_response(norm=norm)
        
        plot(self.frequencies, self.response)
        title("ENBW=%2.1f" % (self.enbw))
        ylabel('Frequency response (dB)')
        xlabel('Fraction of sampling frequency')
        # define the plot limits 
        xlim(-0.5, 0.5)
        y0, y1 = ylim()
        if mindB:
            y0 = mindB
        if maxdB != None:
            y1 = maxdB
        else:
            y1 = max(self.response)
            
        ylim(y0, y1)
        
        grid(True)

    def plot_window(self):
        """Plot the window in the time domain

        .. plot::
            :width: 80%
            :include-source:

            from spectrum.window import Window
            w = Window(64, name='hamming')
            w.plot_window()

        """
        from pylab import plot, xlim, grid, title, ylabel, axis
        x = linspace(0, 1, self.N)
        xlim(0, 1)
        plot(x, self.data)
        grid(True)
        title('%s Window (%s points)' % (self.name.capitalize(), self.N))
        ylabel('Amplitude')
        axis([0, 1, 0, 1.1])

    def plot_time_freq(self, mindB=-100, maxdB=None, norm=True):
        """Plotting method to plot both time and frequency domain results.

        See :meth:`plot_frequencies` for the optional arguments.
        
        .. plot::
            :width: 80%
            :include-source:

            from spectrum.window import Window
            w = Window(64, name='hamming')
            w.plot_time_freq()

        """
        from pylab import subplot

        subplot(1, 2, 1)
        self.plot_window()
        
        subplot(1, 2, 2)
        self.plot_frequencies(mindB=mindB, maxdB=maxdB, norm=norm)
        
    def info(self):
        """Print object information such as length and name"""
        print self
        
    def __str__(self):
        msg = 'Window object:\n'
        msg += 'Name: %s\n' % self.name
        msg += 'Length: %s\n' % self.N
        msg += 'NFFT (for the frequency response): %s\n' % self.response.size
        msg += 'ENBW=%s\n' % self.enbw
        return msg
        

def create_window(N, name=None, **kargs):
    r"""Returns the N-point window given a valid name

    :param int N: window size
    :param str name: window name (default is *rectangular*). Valid names 
        are stored in :func:`~spectrum.window.window_names`.
    :param kargs: optional arguments are:
    
        * *beta*: argument of the :func:`window_kaiser` function (default is 8.6)
        * *attenuation*: argument of the :func:`window_chebwin` function (default is 50dB)
        * *alpha*: argument of the 
            1. :func:`window_gaussian` function (default is 2.5)
            2. :func:`window_blackman` function (default is 0.16)
            3. :func:`window_poisson` function (default is 2)
            4. :func:`window_cauchy` function (default is 3) 
        * *mode*: argument :func:`window_flattop` function (default is *symmetric*, can be *periodic*)
        * *r*: argument of the :func:`window_tukey` function (default is 0.5).
        
    The following windows have been simply wrapped from existing librairies like
    NumPy:

        * **Rectangular**: :func:`window_rectangle`,
        * **Bartlett** or Triangular: see :func:`window_bartlett`,
        * **Hanning** or Hann: see :func:`window_hann`,
        * **Hamming**: see :func:`window_hamming`,
        * **Kaiser**: see :func:`window_kaiser`,
        * **chebwin**: see :func:`window_chebwin`.    
        
    The following windows have been implemented from scratch:
        
        * **Blackman**: See :func:`window_blackman`
        * **Bartlett-Hann** : see :func:`window_bartlett_hann`
        * **cosine or sine**: see :func:`window_cosine`
        * **gaussian**: see :func:`window_gaussian`
        * **Bohman**: see :func:`window_bohman`
        * **Lanczos or sinc**: see :func:`window_lanczos`
        * **Blackman Harris**: see :func:`window_blackman_harris`
        * **Blackman Nuttall**: see :func:`window_blackman_nuttall`
        * **Nuttall**: see :func:`window_nuttall`
        * **Tukey**: see :func:`window_tukey`
        * **Parzen**: see :func:`window_parzen`
        * **Flattop**: see :func:`window_flattop`
        * **Riesz**: see :func:`window_riesz`
        * **Riemann**: see :func:`window_riemann`
        * **Poisson**: see :func:`window_poisson`
        * **Poisson-Hanning**: see :func:`window_poisson_hanning`
    
    .. todo:: on request taylor, potter, Bessel, expo, 
        rife-vincent, Kaiser-Bessel derived (KBD)    

    .. plot::
        :width: 80%
        :include-source:

        from pylab import plot, hold, legend
        from spectrum import create_window

        data = create_window(51, 'hamming')
        plot(data, label='hamming')
        hold(True)
        data = create_window(51, 'kaiser')
        plot(data, label='kaiser')
        legend()

    .. plot::
        :width: 80%
        :include-source:

        from pylab import *
        from spectrum import create_window

        A = fft(create_window(51, 'hamming'), 2048) / 25.5
        mag = abs(fftshift(A))
        freq = linspace(-0.5,0.5,len(A))
        response = 20*log10(mag)
        mindB = -60
        response = clip(response,mindB,100)
        plot(freq, response)

    .. seealso:: :func:`window_visu`, :func:`Window`, :mod:`spectrum.dpss`
    """
    if name == None:
        name = 'rectangle'
    name = name.lower()
    assert name in window_names.keys(), \
        """window name %s not implemented or incorrect. Try to use one of %s"""\
        % (name, window_names)
    
        
    # create the function name
    f = eval(window_names[name])

    windows_with_parameters = \
    {'kaiser': {'beta': eval(window_names['kaiser']).__defaults__[0]}, 
     'blackman': {'alpha': eval(window_names['blackman']).__defaults__[0]},
     'cauchy': {'alpha': eval(window_names['cauchy']).__defaults__[0]},
     'flattop': {'mode': eval(window_names['flattop']).__defaults__[0]},
     'gaussian': {'alpha': eval(window_names['gaussian']).__defaults__[0]},
     'chebwin': {'attenuation':eval(window_names['chebwin']).__defaults__[0]},
     'tukey': {'r':eval(window_names['tukey']).__defaults__[0]},
     'poisson': {'alpha': eval(window_names['poisson']).__defaults__[0]},
     'poisson_hanning': {'alpha': 
                         eval(window_names['poisson_hanning']).__defaults__[0]},
         }
    
    if name not in windows_with_parameters.keys():
        if len(kargs) == 0:
            # not parameters, so we directly call the function    
            w = f(N)
        else:
            raise ValueError("""
            Parameters do not match any of the window. The window provided 
            do not expect any parameters. Try to remove the parameters""")
    elif name in windows_with_parameters.keys():
        # user optional parameters are provided, scan them:
        dargs = {}
        for arg in kargs.keys(): 
            # check that the parameters are valid, and obtain the default value
            try:
                default = windows_with_parameters[name][arg]
            except:
                raise ValueError("""
                    Invalid optional argument (%s) for %s window.
                    Valid optional arguments are (%s)""" % \
                    (arg, name, windows_with_parameters[name].keys()))
            # add the user parameter to the list of parameters
            dargs[arg] = kargs.get(arg, default)
        # call the proper function with the optional arguments
        w = f(N, **dargs)
        
    return w


def enbw(data):
    r"""Computes the equivalent noise bandwidth 

    .. math:: ENBW = N \frac{\sum_{n=1}^{N} w_n^2}{\left(\sum_{n=1}^{N} w_n \right)^2}

    .. doctest::

        >>> from spectrum import *
        >>> w = create_window(64, 'rectangular')
        
        >>> enbw(w)
        1.0

    The following table contains the ENBW values for some of the 
    implemented windows in this module (with N=16384). They have been 
    double checked against litterature (Source: [Harris]_, [Marple]_). 
    
    If not present, it means that it has not been checked.
    
    =================== ============ =============
    name                 ENBW        litterature
    =================== ============ =============
    rectangular         1.           1. 
    triangle            1.3334       1.33
    Hann                1.5001       1.5
    Hamming             1.3629       1.36
    blackman            1.7268       1.73
    kaiser              1.7
    blackmanharris,4    2.004        2.
    riesz               1.2000       1.2
    riemann             1.32         1.3
    parzen              1.917        1.92
    tukey 0.25          1.102        1.1
    bohman              1.7858       1.79
    poisson 2           1.3130       1.3
    hanningpoisson 0.5  1.609        1.61
    cauchy              1.489        1.48
    lanczos             1.3
    =================== ============ =============
    
    
    """
    N = len(data)
    return N * sum(data**2) / sum(data)**2


def _kaiser(n, beta):
    """Independant Kaiser window

    For the definition of the Kaiser window, see A. V. Oppenheim & R. W. Schafer, "Discrete-Time Signal Processing".

    The continuous version of width n centered about x=0 is:

    .. note:: 2 times slower than scipy.kaiser
    """
    from scipy.special import iv as besselI
    m = n - 1
    k = arange(0, m)
    k = 2. * beta / m * sqrt (k * (m - k))
    w = besselI (0, k) / besselI (0, beta)
    return w


def window_visu(N=51, name='hamming', **kargs):
    """A Window visualisation tool

    :param N: length of the window
    :param name: name of the window
    :param NFFT: padding used by the FFT
    :param mindB: the minimum frequency power in dB 
    :param maxdB: the maximum frequency power in dB
    :param kargs: optional arguments passed to :func:`create_window`

    This function plot the window shape and its equivalent in the Fourier domain.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'kaiser', beta=8.)

    """
    # get the default parameters
    mindB = kargs.pop('mindB', -100)
    maxdB = kargs.pop('maxdB', None)
    norm = kargs.pop('norm', True)
     
    # create a window object
    w = Window(N, name, **kargs)
    
    # plot the time and frequency windows
    w.plot_time_freq(mindB=mindB, maxdB=maxdB, norm=norm)


def window_rectangle(N):
    r"""Kaiser window

    :param N: window length
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'rectangle')
        
    """
    return ones(N)

    
def window_kaiser(N, beta=8.6, method='numpy'):
    r"""Kaiser window

    :param N: window length
    :param beta: kaiser parameter (default is 8.6)

    To obtain a Kaiser window that designs an FIR filter with 
    sidelobe attenuation of :math:`\alpha` dB, use the following :math:`\beta` where
    :math:`\beta = \pi \alpha`.

    .. math::

        w_n = \frac{I_0\left(\pi\alpha\sqrt{1-\left(\frac{2n}{M}-1\right)^2}\right)} {I_0(\pi \alpha)}

    where

      * :math:`I_0` is the zeroth order Modified Bessel function of the first kind.
      * :math:`\alpha` is a real number that determines the shape of the window. It determines 
        the trade-off between main-lobe width and side lobe level.
      * the length of the sequence is N=M+1.

    The Kaiser window can approximate many other windows by varying the :math:`\beta` parameter

    ===== ========================
    beta  Window shape
    ===== ========================
    0     Rectangular
    5     Similar to a Hamming
    6     Similar to a Hanning
    8.6   Similar to a Blackman
    ===== ========================

    .. plot::
        :width: 80%
        :include-source:

        from pylab import plot, legend, hold, xlim
        from spectrum import window_kaiser
        N = 64
        for beta in [1,2,4,8,16]:
            plot(window_kaiser(N, beta), label='beta='+str(beta))
            hold(True)
        xlim(0,N)
        legend()


    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'kaiser', beta=8.)

    .. seealso:: numpy.kaiser, :func:`spectrum.window.create_window`
    """
    if N == 1:
        return ones(1)
    if method == 'numpy':
        from numpy import kaiser 
        return kaiser(N, beta)
    else:
        return _kaiser(N, beta)
    

def window_blackman(N, alpha=0.16):
    r"""Blackman window

    :param N: window length

    .. math:: a_0 - a_1 \cos(\frac{2\pi n}{N-1}) +a_2 \cos(\frac{4\pi n }{N-1})

    with 

    .. math::
    
        a_0 = (1-\alpha)/2, a_1=0.5, a_2=\alpha/2 \rm{\;and\; \alpha}=0.16

    When :math:`\alpha=0.16`, this is the unqualified Blackman window with
    :math:`a_0=0.48`  and :math:`a_2=0.08`.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'blackman')
    
    .. note:: Although Numpy implements a blackman window for :math:`\alpha=0.16`, 
        this implementation is valid for any :math:`\alpha`.
            
    .. seealso:: numpy.blackman, :func:`create_window`, :class:`Window`
    
    """    
    a0 = (1. - alpha)/2.
    a1 = 0.5
    a2 = alpha/2.
    
    if (N == 1):
        win = array([1.])
    else:
        k = arange(0, N)/float(N-1.)
        win =  a0 - a1 * cos (2 * pi * k) + a2 * cos (4 * pi * k)
    return win

    
    #from numpy import blackman
    #return blackman(N)


def window_bartlett(N):
    r"""Bartlett window (wrapping of numpy.bartlett) also known as Fejer

    :param int N: window length

    The Bartlett window is defined as
    
    .. math:: w(n) = \frac{2}{N-1} \left(
              \frac{N-1}{2} - \left|n - \frac{N-1}{2}\right|
              \right)

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'bartlett')
        
    .. seealso:: numpy.bartlett, :func:`create_window`, :class:`Window`.
    """
    from numpy import bartlett
    return bartlett(N)


def window_hamming(N):
    r"""Hamming window

    :param N: window length
    

    The Hamming window is defined as

    .. math:: 0.54 -0.46 \cos\left(\frac{2\pi n}{N-1}\right)
               \qquad 0 \leq n \leq M-1

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'hamming')
        
    .. seealso:: numpy.hamming, :func:`create_window`, :class:`Window`.
    """
    from numpy import hamming
    return hamming(N)


def window_hann(N):
    r"""Hann window (or Hanning). (wrapping of numpy.bartlett)

    :param int N: window length

    The Hanning window is also known as the Cosine Bell. Usually, it is called
    Hann window, to avoid confusion with the Hamming window.

    .. math:: w(n) =  0.5\left(1- \cos\left(\frac{2\pi n}{N-1}\right)\right)
               \qquad 0 \leq n \leq M-1

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'hanning')
        
    .. seealso:: numpy.hanning, :func:`create_window`, :class:`Window`.
    """
    from numpy import hanning
    return hanning(N)


def window_gaussian(N, alpha=2.5):
    r"""Gaussian window

    :param N: window length

    .. math:: \exp^{-0.5 \left( \sigma\frac{n}{N/2} \right)^2}
    
    with :math:`\frac{N-1}{2}\leq n \leq \frac{N-1}{2}`.
    
    .. note:: N-1 is used to be in agreement with octave convention. The ENBW of
         1.4 is also in agreement with [Harris]_

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'gaussian', alpha=2.5)

    

    .. seealso:: scipy.signal.gaussian, :func:`create_window`
    """
    t = linspace(-(N-1)/2., (N-1)/2., N)
    #t = linspace(-(N)/2., (N)/2., N)
    w = exp(-0.5*(alpha * t/(N/2.))**2.)
    return w


def window_chebwin(N, attenuation=50):
    """Cheb window

    :param N: window length

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'chebwin', attenuation=50)

    .. seealso:: scipy.signal.chebwin, :func:`create_window`, :class:`Window`
    """
    import scipy.signal
    return scipy.signal.chebwin(N, attenuation)


def window_cosine(N):
    r"""Cosine tapering window also known as sine window.

    :param N: window length
    
    .. math:: w(n) = \cos\left(\frac{\pi n}{N-1} - \frac{\pi}{2}\right) = \sin \left(\frac{\pi n}{N-1}\right)

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'cosine')

    .. seealso:: :func:`create_window`, :class:`Window`
    """
    if N ==1: 
        return ones(1)
    n = arange(0, N)
    win = sin(pi*n/(N-1.))
    return win

def window_lanczos(N):
    r"""Lanczos window also known as sinc window.

    :param N: window length
    
    .. math:: w(n) = sinc \left(  \frac{2n}{N-1} - 1 \right)

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'lanczos')
    
    .. seealso:: :func:`create_window`, :class:`Window`
    """
    if N ==1: 
        return ones(1)
    
    n = linspace(-N/2., N/2., N)
    win = sinc(2*n/(N-1.))
    return win


def window_bartlett_hann(N):
    r"""Bartlett-Hann window

    :param N: window length
    
    .. math:: w(n) = a_0 + a_1 \left| \frac{n}{N-1} -\frac{1}{2}\right| - a_2 \cos \left( \frac{2\pi n}{N-1} \right)

    with :math:`a_0 = 0.62`, :math:`a_1 = 0.48` and :math:`a_2=0.38`

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'bartlett_hann')

    .. seealso:: :func:`create_window`, :class:`Window`
    """
    if N == 1:
        return ones(1)
    n = arange(0, N)

    a0 = 0.62
    a1 = 0.48
    a2 = 0.38

    win = a0 -  a1 *abs(n/(N-1.)-0.5) -a2 * cos(2*pi*n/(N-1.))
    return win



def _coeff4(N, a0, a1, a2, a3):
    """a common internal function to some window functions with 4 coeffs


    For the blackmna harris for instance, the results are identical to octave if N is odd
    but not for even values...if n =0 whatever N is, the w(0) must be equal to a0-a1+a2-a3, which 
    is the case here, but not in octave..."""
    if N == 1:
        return ones(1)
    
    n = arange(0, N)
    N1 = N - 1.

    w = a0 -a1*cos(2.*pi*n / N1) + a2*cos(4.*pi*n / N1) - a3*cos(6.*pi*n / N1)

    return w

def window_nuttall(N):
    r"""Nuttall tapering window

    :param N: window length
    
    .. math:: w(n) = a_0 - a_1 \cos\left(\frac{2\pi n}{N-1}\right)+ a_2 \cos\left(\frac{4\pi n}{N-1}\right)- a_3 \cos\left(\frac{6\pi n}{N-1}\right)

    with :math:`a_0 = 0.355768`, :math:`a_1 = 0.487396`, :math:`a_2=0.144232` and :math:`a_3=0.012604`

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'nuttall', mindB=-80)


    .. seealso:: :func:`create_window`, :class:`Window`
    """
    a0 = 0.355768
    a1 = 0.487396
    a2 = 0.144232
    a3 = 0.012604

    return _coeff4(N, a0, a1, a2, a3)


def window_blackman_nuttall(N):
    r"""Blackman Nuttall window

    returns a minimum, 4-term Blackman-Harris window. The window is minimum in the sense that its maximum sidelobes are minimized. 
    The coefficients for this window differ from the Blackman-Harris window coefficients and produce slightly lower sidelobes.

    :param N: window length
    
    .. math:: w(n) = a_0 - a_1 \cos\left(\frac{2\pi n}{N-1}\right)+ a_2 \cos\left(\frac{4\pi n}{N-1}\right)- a_3 \cos\left(\frac{6\pi n}{N-1}\right)

    with :math:`a_0 = 0.3635819`, :math:`a_1 = 0.4891775`, :math:`a_2=0.1365995` and :math:`0_3=.0106411`

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'blackman_nuttall', mindB=-80)

    .. seealso:: :func:`spectrum.window.create_window`
    .. seealso:: :func:`create_window`, :class:`Window`

    """
    a0 = 0.3635819
    a1 = 0.4891775
    a2 = 0.1365995
    a3 = 0.0106411
    return _coeff4(N, a0, a1, a2, a3)


def window_blackman_harris(N):
    r"""Blackman Harris window

    :param N: window length
    
    .. math:: w(n) = a_0 - a_1 \cos\left(\frac{2\pi n}{N-1}\right)+ a_2 \cos\left(\frac{4\pi n}{N-1}\right)- a_3 \cos\left(\frac{6\pi n}{N-1}\right)

    =============== =========
    coeff            value
    =============== =========
    :math:`a_0`     0.35875
    :math:`a_1`     0.48829
    :math:`a_2`     0.14128
    :math:`a_3`     0.01168
    =============== =========

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'blackman_harris', mindB=-80)

    .. seealso:: :func:`spectrum.window.create_window`
    .. seealso:: :func:`create_window`, :class:`Window`
    """
    a0 = 0.35875
    a1 = 0.48829
    a2 = 0.14128
    a3 = 0.01168
    return _coeff4(N, a0, a1, a2, a3)


def window_bohman(N):
    r"""Bohman tapering window

    :param N: window length
    
    .. math:: w(n) = (1-|x|) \cos (\pi |x|) + \frac{1}{\pi} \sin(\pi |x|)

    where x is a length N vector of linearly spaced values between 
    -1 and 1. 
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'bohman')
        
    .. seealso:: :func:`create_window`, :class:`Window`
    """
    x = linspace(-1, 1, N)
    w = (1.-abs(x)) * cos(pi*abs(x)) + 1./pi * sin(pi*abs(x))
    return w



def window_tukey(N, r=0.5):
    """Tukey tapering window (or cosine-tapered window) 
    
    :param N: window length
    :param r: defines the ratio between the constant section and the cosine 
      section. It has to be between 0 and 1. 
      
    The function returns a Hanning window for `r=0` and a full box for `r=1`.
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'tukey')
        window_visu(64, 'tukey', r=1)

    .. math:: 0.5 (1+cos(2pi/r (x-r/2))) for 0<=x<r/2
    
    .. math:: 0.5 (1+cos(2pi/r (x-1+r/2))) for x>r/2
    
    
    .. seealso:: :func:`create_window`, :class:`Window`
    """
    assert r>=0 and r<=1 , "r must be in [0,1]"
    if N==1:
        return ones(1)

    if r == 0:
        return ones(N)
    elif r == 1:
        return window_hann(N)
    else:
        from numpy import flipud, concatenate, where
        
        ## cosine-tapered window
        x = linspace(0, 1, N)
        x1 = where(x<r/2.)
        w = 0.5*(1+cos(2*pi/r*(x[x1[0]]-r/2)))
        w = concatenate((w, ones(N-len(w)*2), flipud(w)))


        return w


def window_parzen(N):
    r"""Parsen tapering window (also known as de la Valle-Poussin)
    
    :param N: window length

    Parzen windows are piecewise cubic approximations 
    of Gaussian windows. Parzen window sidelobes fall off as :math:`1/\omega^4`.
    
    if :math:`0\leq|x|\leq (N-1)/4`:
    
    .. math:: w(n) = 1-6 \left( \frac{|n|}{N/2} \right)^2 +6 \left( \frac{|n|}{N/2}\right)^3
    
    if :math:`(N-1)/4\leq|x|\leq (N-1)/2`
    
    .. math:: w(n) = 2 \left(1- \frac{|n|}{N/2}\right)^3
    
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'parzen')


    .. seealso:: :func:`create_window`, :class:`Window`
    """
    from numpy import  where, concatenate

    n = linspace(-(N-1)/2., (N-1)/2., N)
    n1 = n[where(abs(n)<=(N-1)/4.)[0]]
    n2 = n[where(n>(N-1)/4.)[0]]
    n3 = n[where(n<-(N-1)/4.)[0]]
    

    w1 = 1. -6.*(abs(n1)/(N/2.))**2 + 6*(abs(n1)/(N/2.))**3
    w2 = 2.*(1-abs(n2)/(N/2.))**3
    w3 = 2.*(1-abs(n3)/(N/2.))**3

    w = concatenate((w3, w1, w2))
    return w


def window_flattop(N, mode='symmetric',precision=None):
    r"""Flat-top tapering window
    
    Returns symmetric or periodic flat top window.
      
    :param N: window length
    :param mode: way the data are normalised. If mode is *symmetric*, then
        divide n by N-1. IF mode is *periodic*, divide by N, 
        to be consistent with octave code.
    
    When using windows for filter design, the *symmetric* mode 
    should be used (default). When using windows for spectral analysis, the *periodic* 
    mode should be used. The mathematical form of the flat-top window in the symmetric
    case is:

    .. math:: w(n) = a_0 
        - a_1 \cos\left(\frac{2\pi n}{N-1}\right)
        + a_2 \cos\left(\frac{4\pi n}{N-1}\right)
        - a_3 \cos\left(\frac{6\pi n}{N-1}\right)
        + a_4 \cos\left(\frac{8\pi n}{N-1}\right)
        
    =====  =============
    coeff  value
    =====  =============
    a0     0.21557895
    a1     0.41663158
    a2     0.277263158
    a3     0.083578947
    a4     0.006947368
    =====  =============
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'bohman')

    
    .. seealso:: :func:`create_window`, :class:`Window`
    """
    assert mode in ['periodic', 'symmetric']
    t = arange(0, N)
    if mode == 'periodic':
        x = 2*pi*t/float(N)
    else:
        if N ==1: 
            return ones(1)
        x = 2*pi*t/float(N-1)
    a0 = 0.21557895
    a1 = 0.41663158
    a2 = 0.277263158
    a3 = 0.083578947
    a4 = 0.006947368

    if precision == 'octave':
        #to compare with octave, same as above but less precise
        d  = 4.6402
        a0 = 1./d
        a1 = 1.93/d
        a2 = 1.29/d
        a3 = 0.388/d
        a4 = 0.0322/d
    w = a0-a1*cos(x)+a2*cos(2*x)-a3*cos(3*x)+a4*cos(4*x)
    return w


def _window_taylor(N, nbar=4, sll=-30):
    """Taylor tapering window
    
    Taylor windows allows you to make tradeoffs between the 
    mainlobe width and sidelobe level (sll). 
    
    :param N: window length
    :param float nbar:
    :param float sll:
    
    The default values gives equal height 
    sidelobes (nbar) and maximum sidelobe level (sll). 
    
    .. warning:: not implemented
    
    .. seealso:: :func:`create_window`, :class:`Window`
    """
    raise NotImplementedError


def window_riesz(N):
    r"""Riesz tapering window
    
    :param N: window length
    
    .. math:: w(n) = 1 - \left| \frac{n}{N/2}  \right|^2
    
    with :math:`-N/2 \leq n \leq N/2`.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'riesz')

    .. seealso:: :func:`create_window`, :class:`Window`
    """
    n = linspace(-N/2., (N)/2., N)
    
    w = 1 - abs(n/(N/2.))**2.
    return w


def window_riemann(N):
    r"""Riemann tapering window
    
    :param int N: window length
    
    .. math:: w(n) = 1 - \left| \frac{n}{N/2}  \right|^2
    
    with :math:`-N/2 \leq n \leq N/2`.

    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'riesz')

    .. seealso:: :func:`create_window`, :class:`Window`
    """
    n = linspace(-N/2., (N)/2., N)
    w = sin(n/N*2.*pi)/(n/N*2.*pi)
    return w



def window_poisson(N, alpha=2):
    r"""Poisson tapering window

    :param int N: window length
    
    .. math:: w(n) = \exp^{-\alpha \frac{|n|}{N/2} }
    
    with :math:`-N/2 \leq n \leq N/2`.
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'poisson')
        window_visu(64, 'poisson', alpha=3)
        window_visu(64, 'poisson', alpha=4)

    .. seealso:: :func:`create_window`, :class:`Window`
    """
    n = linspace(-N/2., (N)/2., N)
    w = exp(-alpha * abs(n)/(N/2.))
    return w


def window_poisson_hanning(N, alpha=2):
    r"""Hann-Poisson tapering window
    
    This window is constructed as the product of the Hanning and Poisson
    windows. The parameter **alpha** is the Poisson parameter.
    
    :param int N: window length
    :param float alpha: parameter of the poisson window
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'poisson_hanning', alpha=0.5)
        window_visu(64, 'poisson_hanning', alpha=1)
        window_visu(64, 'poisson_hanning')

    .. seealso:: :func:`window_poisson`, :func:`window_hann`
    """
    w1 = window_hann(N)
    w2 = window_poisson(N, alpha=alpha)
    return w1*w2


def window_cauchy(N, alpha=3):
    r"""Cauchy tapering window
    
    :param int N: window length
    :param float alpha: parameter of the poisson window
    
    .. math:: w(n) = \frac{1}{1+\left(\frac{\alpha*n}{N/2}\right)**2}
    
    .. plot::
        :width: 80%
        :include-source:

        from spectrum import window_visu
        window_visu(64, 'cauchy', alpha=3)
        window_visu(64, 'cauchy', alpha=4)
        window_visu(64, 'cauchy', alpha=5)
        
        
    .. seealso:: :func:`window_poisson`, :func:`window_hann`
    """
    n = linspace(-N/2., (N)/2., N)
    w = 1./(1.+ (alpha*n/(N/2.))**2)
    return w


