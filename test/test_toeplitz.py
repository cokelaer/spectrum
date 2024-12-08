import numpy
import pytest

from spectrum import toeplitz
from spectrum.toeplitz import HERMTOEP, TOEPLITZ


def test_toeplitz():
    T0 = 3 + 0j
    TC = numpy.array([-2.0 + 0.5j, 0.7 - 1j])
    TR = numpy.array([-0.2 - 0.4j, 0.3 - 0.6j])
    Z = numpy.array([1.0 + 3.0j, 2.0 - 1.0j, 0.5 + 0.8j])
    X = TOEPLITZ(T0, TC, TR, Z)
    expected_X = numpy.array([0.23518549 + 1.24372203j, 1.03025807 + 0.53575117j, 0.47334662 + 0.24031779j])
    assert X == pytest.approx(expected_X)


def test_hermtoep():
    T0 = 3.0
    T = numpy.array([-2.0 + 0.5j, 0.7 - 1j])
    Z = numpy.array([1.0 + 3.0j, 2.0 - 1.0j, 0.5 + 0.8j])
    X = HERMTOEP(T0, T, Z)
    expected_X = numpy.array([2.29697452 + 1.69904459j, 3.31847134 + 1.29697452j, 1.49283439 + 0.94745223j])
    assert X == pytest.approx(expected_X)
