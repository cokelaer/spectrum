"""
.. topic:: The :mod:`datasets` module provides data sets to test the
    **Spectrum** functionalities.

    .. autosummary::

        data_cosine
        marple_data
        TimeSeries

    .. codeauthor:: Thomas Cokelaer 2011

    :Reference: [Marple]_
"""
import numpy as np
from numpy import arange, cos, pi, sin
from numpy.random import randn

# example marple app 3.A
A_cholesky = np.array(
    [[2 + 0.0j, 0.5 - 0.5j, -0.2 + 0.1j], [0.5 + 0.5j, 1, 0.3 - 0.2j], [-0.2 - 0.1j, 0.3 + 0.2j, 0.5]], dtype=complex
)
a_cholesky = np.array([2 + 0.0j, 0.5 - 0.5j, 1.0, -0.2 + 0.1j, 0.3 - 0.2j, 0.5], dtype=complex)
B_cholesky = np.array([1 + 3j, 2 - 1j, 0.5 + 0.8j], dtype=complex)
# should return
sol_cholesky = np.array([0.95945946 + 5.25675676j, 4.41891892 - 7.04054054j, -5.13513514 + 6.35135135j])


#: 64-complex data length from Marple reference [Marple]_
marple_data = [
    1.349839091 + 2.011167288j,
    -2.117270231 + 0.817693591j,
    -1.786421657 - 1.291698933j,
    1.162236333 - 1.482598066j,
    1.641072035 + 0.372950256j,
    0.072213709 + 1.828492761j,
    -1.564284801 + 0.824533045j,
    -1.080565453 - 1.869776845j,
    0.927129090 - 1.743406534j,
    1.891979456 + 0.972347319j,
    -0.105391249 + 1.602209687j,
    -1.618367076 + 0.637513280j,
    -0.945704579 - 1.079569221j,
    1.135566235 - 1.692269921j,
    1.855816245 + 0.986030221j,
    -1.032083511 + 1.414613724j,
    -1.571600199 + 0.089229003j,
    -0.243143231 - 1.444692016j,
    0.838980973 - 0.985756695j,
    1.516003132 + 0.928058863j,
    0.257979959 + 1.170676708j,
    -2.057927608 + 0.343388647j,
    -0.578682184 - 1.441192508j,
    1.584011555 - 1.011150956j,
    0.614114344 + 1.508176208j,
    -0.710567117 + 1.130144477j,
    -1.100205779 - 0.584209621j,
    0.150702029 - 1.217450142j,
    0.748856127 - 0.804411888j,
    0.795235813 + 1.114466429j,
    -0.071512341 + 1.017092347j,
    -1.732939839 - 0.283070654j,
    0.404945314 - 0.781708360j,
    1.293794155 - 0.352723092j,
    -0.119905084 + 0.905150294j,
    -0.522588372 + 0.437393665j,
    -0.974838495 - 0.670074046j,
    0.275279552 - 0.509659231j,
    0.854210198 - 0.008278057j,
    0.289598197 + 0.506233990j,
    -0.283553183 + 0.250371397j,
    -0.359602571 - 0.135261074j,
    0.102775671 - 0.466086507j,
    -0.009722650 + 0.030377999j,
    0.185930878 + 0.808869600j,
    -0.243692726 - 0.200126961j,
    -0.270986766 - 0.460243553j,
    0.399368525 + 0.249096692j,
    -0.250714004 - 0.362990230j,
    0.419116348 - 0.389185309j,
    -0.050458215 + 0.702862442j,
    -0.395043731 + 0.140808776j,
    0.746575892 - 0.126762003j,
    -0.559076190 + 0.523169816j,
    -0.344389260 - 0.913451135j,
    0.733228028 - 0.006237417j,
    -0.480273813 + 0.509469569j,
    0.033316225 + 0.087501869j,
    -0.321229130 - 0.254548967j,
    -0.063007891 - 0.499800682j,
    1.239739418 - 0.013479125j,
    0.083303742 + 0.673984587j,
    -0.762731433 + 0.408971250j,
    -0.895898521 - 0.364855707j,
]


def data_cosine(N=1024, A=0.1, sampling=1024.0, freq=200):
    r"""Return a noisy cosine at a given frequency.

    :param N:           the final data size
    :param A:           the strength of the noise
    :param float sampling: sampling frequency of the input :attr:`data`.
    :param float freq:  the frequency :math:`f_0` of the cosine.

    .. math:: x[t] = cos(2\pi t * f_0) + A w[t]

    where w[t] is a white noise of variance 1.

    .. doctest::

        >>> from spectrum import data_cosine
        >>> a = data_cosine(N=1024, sampling=1024, A=0.5, freq=100)

    """
    t = arange(0, float(N) / sampling, 1.0 / sampling)
    x = cos(2.0 * pi * t * freq) + A * randn(t.size)
    return x


def data_two_freqs(N=200):
    """A simple test example with two close frequencies"""
    nn = arange(N)
    xx = cos(0.257 * pi * nn) + sin(0.2 * pi * nn) + 0.01 * randn(nn.size)
    return xx


def spectrum_data(filename):
    """Simple utilities to retrieve data sets from"""
    import inspect
    from pathlib import Path

    import spectrum

    spectrum_path = Path(inspect.getfile(spectrum)).parent

    filename = spectrum_path / "data" / filename
    if filename.exists():
        return str(filename)
    else:
        raise Exception("unknown file %s" % filename)


#: filename of a WAV data file 150,000 data points
dolphin_filename = spectrum_data("DOLPHINS.wav")


class TimeSeries(object):
    """A simple Base Class for various data sets.

    .. doctest::

        >>> from spectrum import TimeSeries
        >>> data = [1, 2, 3, 4, 3, 2, 1, 0 ]
        >>> ts = TimeSeries(data, sampling=1)
        >>> ts.plot()
        >>> ts.dt
        1.0

    """

    def __init__(self, data, sampling=1):
        """

        :param array data: input data (list or numpy.array)
        :param sampling: the sampling frequency of the data (default 1Hz)

        """
        self.data = data
        self.N = len(data)
        self.sampling = sampling
        self.dt = 1.0 / sampling

    def plot(self, **kargs):
        """Plot the data set, using the sampling information to set the x-axis
        correctly."""
        from pylab import grid, linspace, plot, xlabel, ylabel

        time = linspace(1 * self.dt, self.N * self.dt, self.N)
        plot(time, self.data, **kargs)
        xlabel("Time")
        ylabel("Amplitude")
        grid(True)
