from spectrum import CORRELATION, marple_data, xcorr
from nose.tools import assert_almost_equal
from numpy.testing import assert_array_almost_equal
from numpy import array
from spectrum.correlation import * # just for nosetests to scan the docstrings

def test_CORRELATION():
    R15 = CORRELATION(marple_data, maxlags=15, norm='biased')
    R15 = CORRELATION(marple_data, maxlags=15, norm='unbiased')
    R15 = CORRELATION(marple_data, maxlags=15, norm='coeff')
    R15 = CORRELATION(marple_data, maxlags=15, norm=None)



def test_CORRELATION_biased():
    R15 = CORRELATION(marple_data, maxlags=15, norm='biased')

    assert_almost_equal(R15[0], 1.7804598944893049+0j)
    assert_almost_equal(R15[1], 0.32076613+1.50586147j)
    assert_almost_equal(R15[2], -1.29947785+0.74815755j, places=6)


    R30 = CORRELATION(marple_data, maxlags=30, norm='biased')
    assert_almost_equal(R30[0], 1.7804598944893049+0j, places=7)
    assert_almost_equal(R30[1], R15[1])
    assert_almost_equal(R30[2], R15[2])

def test_CORRELATION_unbiased():
    R15 = CORRELATION(marple_data, maxlags=15, norm='unbiased')
    assert_almost_equal(R15[0], 1.7804598944893049+0j)
    assert_almost_equal(R15[1],  0.32585765+1.52976403j)
    assert_almost_equal(R15[2],  -1.34139649+0.77229166j)


    R30 = CORRELATION(marple_data, maxlags=30, norm='unbiased')
    assert_almost_equal(R30[0], 1.7804598944893049+0j, places=7)
    assert_almost_equal(R30[1], R15[1])
    assert_almost_equal(R30[2], R15[2])


def test_xcorr():

    x = array([1,2,3,4,5])
    corr,l = xcorr(x, x, maxlags=4, norm='coeff')
    corr_coeff = array([ 0.09090909,  0.25454545,  0.47272727,  0.72727273,  1.        ,    0.72727273,  0.47272727,  0.25454545,  0.09090909])
    assert_array_almost_equal(corr, corr_coeff)

    corr,l = xcorr(x, x, maxlags=4, norm='biased')
    corr_biased = array([  1. ,   2.8,   5.2,   8. ,  11. ,   8. ,   5.2,   2.8,   1. ])
    assert_array_almost_equal(corr, corr_biased)

    corr,l = xcorr(x, x, maxlags=4, norm='unbiased')
    corr_unbiased = array([  5.        ,   7.        ,   8.66666667,  10. ,   11.        ,  10.        ,   8.66666667,   7. ,   5. ])
    assert_array_almost_equal(corr, corr_unbiased)

    corr,l = xcorr(x, x, maxlags=4, norm=None)
    corr_none = array([ 5, 14, 26, 40, 55, 40, 26, 14,  5])
    assert_array_almost_equal(corr, corr_none)

    # check default behaviour of maxlags
    corr1,l = xcorr(x, x)
    corr2,l = xcorr(x, x, maxlags=4)
    assert_array_almost_equal(corr1, corr2)

    # check default behaviour of maxlags
    corr1,l = xcorr(x)
    corr2,l = xcorr(x, x)
    assert_array_almost_equal(corr1, corr2)


def test_xcorr_versus_CORRELATION_real_data():
    from spectrum.tools import twosided as func
    x = array([1,2,3,4,5])

    for norm in ['biased', 'unbiased', 'coeff',None]:
        corr1,l = xcorr(x, x, maxlags=4, norm=norm)
        corr2 = CORRELATION(x, x, maxlags=4, norm=norm)
        assert_array_almost_equal(corr1, func(corr2))

def _test_xcorr_versus_CORRELATION_imag_data():
    from spectrum.tools import twosided as func
    x = array([1,2,3,4,5+1.j])

    for norm in ['biased', 'unbiased', 'coeff',None]:
        corr1,l = xcorr(x, x, maxlags=4, norm=norm)
        corr2 = CORRELATION(x, x, maxlags=4, norm=norm)
        assert_array_almost_equal(corr1, func(corr2))
