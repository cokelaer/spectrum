
from spectrum import Spectrogram, dolphin_filename, readwav


def test_spectrogram():

    data, samplerate = readwav(dolphin_filename)
    p = Spectrogram(data, ws=128, W=4096, sampling=samplerate)
    p.periodogram()
    p.plot()

