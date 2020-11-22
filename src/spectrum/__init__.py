from __future__ import absolute_import

try:
    import pkg_resources
    version = pkg_resources.require("spectrum")[0].version
    __version__ = version
except:
    __version__ = "?"



import logging
def spectrum_set_level(level):
    assert level in ['DEBUG', 'INFO', 'CRITICAL', 'ERROR', 'WARNING']
    logging.getLogger().setLevel(level)


#: default number of samples used to compute FFT
default_NFFT = 4096

from .arma import *
from .burg import *
from .cholesky import *
from .correlation import *
from .correlog import *
from .covar import *
from .criteria import *
from .datasets import *
from .eigen import *
from .eigenfre import *
from .io import *
from .levinson import *
from .linear_prediction import *
from .linalg import *
#from lms import *
from .lpc import *
from .minvar import *
from .modcovar import *
from .mtm import *
from .periodogram import *
from .psd import *
from .spectrogram import *
from .tools import *
from .toeplitz import *
from .transfer import *
from .window import *
from .waveform import *
from .yulewalker import *

