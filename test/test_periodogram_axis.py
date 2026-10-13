import numpy as np
import pytest

from spectrum import speriodogram


@pytest.mark.parametrize("axis", [0, 1, -1, -2])
@pytest.mark.parametrize("window", ["hamming", "hann", "rectangle"])
@pytest.mark.parametrize("nfft", [None, 11, 32])
@pytest.mark.parametrize("detrend", [False, True])
@pytest.mark.parametrize("complex_data", [False, True])
def test_matrix_periodogram_matches_separate_channels(
    axis, window, nfft, detrend, complex_data
):
    data = np.random.default_rng(5).normal(size=(13, 3))
    if complex_data:
        data = data + 1j * np.random.default_rng(8).normal(size=data.shape)
    columns = [
        speriodogram(data[:, i], NFFT=nfft, window=window, detrend=detrend)
        for i in range(data.shape[1])
    ]
    expected = np.column_stack(columns)
    if axis in (1, -1):
        data = data.T
        expected = expected.T
    before = data.copy()
    actual = speriodogram(data, NFFT=nfft, window=window, axis=axis, detrend=detrend)
    assert actual.shape == expected.shape
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(data, before)
