import spectrum
from spectrum import *


def test_spectrum():
    assert spectrum.__version__
    assert spectrum.default_NFFT
