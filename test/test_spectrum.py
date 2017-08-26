import spectrum
from spectrum import *


def test_spectrum():
    assert spectrum.default_NFFT
    spectrum_set_level("DEBUG")
