import matplotlib.pyplot as plt
import numpy as np
import pytest

from spectrum import WelchPeriodogram


@pytest.mark.parametrize("sampling", [1, 8000, 44100])
@pytest.mark.parametrize("nfft", [None, 128, 512])
@pytest.mark.parametrize("pad_to", [None, 1023])
def test_welch_spectrum_frequency_metadata(sampling, nfft, pad_to):
    data = np.sin(np.arange(1024) * 0.2)
    (psd, frequencies), estimate = WelchPeriodogram(
        data, NFFT=nfft, sampling=sampling, pad_to=pad_to
    )
    assert estimate.sampling == sampling
    assert estimate.NFFT == (pad_to or nfft or 256)
    np.testing.assert_array_equal(estimate.psd, psd)
    np.testing.assert_allclose(estimate.frequencies(), frequencies, atol=1e-12)
    plt.close("all")
