
__all__ = ['readwav']


def readwav(filename):
    """Read a WAV file and returns the data and sample rate

    ::

        from spectrum.io import readwav
        readwav()

    """
    from scipy.io.wavfile import read as readwav
    samplerate, signal = readwav(filename)
    return signal, samplerate
