from __future__ import absolute_import



import logging
def spectrum_set_level(level):
    assert level in ['DEBUG', 'INFO', 'CRITICAL', 'ERROR', 'WARNING']
    logging.getLogger().setLevel(level)


#: default number os samples used to compute FFT
default_NFFT = 4096

from .mtm import *
from .transfer import *
from .arma import *
from .burg import *
from .cholesky import *
from .correlation import *
from .correlog import *
from .covar import *
from .criteria import *
from .eigen import *
from .eigenfre import *
from .levinson import *
from .linear_prediction import *
from .linalg import *
#from lms import *
from .lpc import *
from .minvar import *
from .modcovar import *
from .periodogram import *
from .psd import *
from .datasets import *
from .tools import *
from .toeplitz import *
from .window import *
from .waveform import *
from .yulewalker import *


