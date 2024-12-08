import numpy
from numpy.testing import assert_array_almost_equal
from pylab import linspace, log10, plot, savefig

from spectrum import (
    data_cosine,
    data_two_freqs,
    ev,
    marple_data,
    music,
    pev,
    pmusic,
    spectrum_set_level,
)
from spectrum.eigenfre import eigen


def test_pmusic():
    p = pmusic(marple_data, 15, NSIG=11)
    p()
    p = pmusic(data_cosine(), 15, NSIG=11, verbose=True)
    p()
    print(p)

    # test verbosity of the _get_signal_space function
    spectrum_set_level("DEBUG")
    pmusic(data_two_freqs(), 15, threshold=1)()
    pmusic(data_two_freqs(), 15, NSIG=11, verbose=True)()
    pmusic(data_two_freqs(), 15, criteria="mdl")()

    #
    pmusic(data_two_freqs(), 15, NSIG=0)()


def test_pev():
    p = pev(marple_data, 15, NSIG=11)
    p()
    p = pev(data_cosine(), 15, NSIG=11, verbose=True)
    p()
    print(p)


def test_eigen_constructor():
    try:
        eigen(marple_data, 8, 4, method="dummy")
        assert False
    except:
        assert True

    # threshold and NSIG cannot be used together
    p = pmusic(marple_data, 15, NSIG=11, threshold=2)

    # NSIG must be less than IP and > 0
    try:
        p = pmusic(marple_data, 15, NSIG=110, threshold=2)
        assert False
    except:
        assert True
    p = pmusic(marple_data, 15, NSIG=-10)

    # NSIG and threshold cannot be provided together
    try:
        eigen(marple_data, 8, NSIG=4, threshold=2)
        assert False
    except:
        assert True

    # NSIG must be positive
    try:
        eigen(marple_data, 8, NSIG=-10)
        assert False
    except:
        assert True

    # NSIG must be less than P (8)
    try:
        eigen(marple_data, 8, NSIG=40)
        assert False
    except:
        assert True


def test_eigenfre_music():

    psd, s = music(marple_data, 15, NSIG=11)


def test_eigenfre_ev():
    psd, s = ev(marple_data, 15, NSIG=11)
    assert_array_almost_equal(
        s,
        numpy.array(
            [
                4.45510959e01,
                1.01451096e01,
                8.37309134e00,
                7.17637043e00,
                6.62545637e00,
                5.83043837e00,
                4.16284271e00,
                3.69224764e00,
                3.64345761e00,
                3.07519938e00,
                2.05618798e00,
                1.53143913e00,
                8.21242005e-01,
                1.10463229e-01,
                1.02225490e-02,
            ]
        ),
    )


def test_eigen_parameters():
    psd, s = ev(data_cosine(), 15)
    psd, s = ev(data_cosine(), 15, NSIG=11)
    psd, s = ev(data_cosine(), 15, threshold=2)


def create_figure():
    psd = test_eigenfre_music()
    f = linspace(-0.5, 0.5, len(psd))
    plot(f, 10 * log10(psd / max(psd)), "--", label="MUSIC 15")
    savefig("psd_eigenfre_music.png")

    psd = test_eigenfre_ev()
    f = linspace(-0.5, 0.5, len(psd))
    plot(f, 10 * log10(psd / max(psd)), "--", label="EV 15")
    savefig("psd_eigenfre_ev.png")


if __name__ == "__main__":
    create_figure()
