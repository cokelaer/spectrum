from spectrum.io import readwav
from spectrum.datasets import dolphin_filename

def test_readwav():
    data, fs = readwav(dolphin_filename)
    assert fs == 22050
