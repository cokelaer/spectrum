from importlib import metadata


def get_package_version(package_name):
    try:
        version = metadata.version(package_name)
        return version
    except metadata.PackageNotFoundError:  # pragma no cover
        return f"{package_name} not found"


version = get_package_version("versionix")

import logging


def spectrum_set_level(level):
    assert level in ["DEBUG", "INFO", "CRITICAL", "ERROR", "WARNING"]
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
from .linalg import *
from .linear_prediction import *

# from lms import *
from .lpc import *
from .minvar import *
from .modcovar import *
from .mtm import *
from .periodogram import *
from .psd import *
from .spectrogram import *
from .toeplitz import *
from .tools import *
from .transfer import *
from .waveform import *
from .window import *
from .yulewalker import *
