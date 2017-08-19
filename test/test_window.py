from spectrum import *
from numpy.testing import assert_almost_equal
import pytest




def test_class_Window():
    w = Window(65, name='hann')
    w.enbw
    w.mean_square
    w.frequencies
    w.plot_time_freq()
    w.compute_response(NFFT=32)

    try:
        w = Window(64, name='wrong')
        assert False
    except:
        assert True
    w.info()
    print(w)

    # recompute response
    w = Window(65, name='hann')
    w.response
    w.plot_frequencies(maxdB=120, mindB=-100)

    try:
        w = Window(65, name="dummy")
        assert False
    except:
        assert True


#unittest of create_window
def test_create_window_error():
    try:
        create_window('dummy')
        assert False
    except:
        assert True


#test that create_window(N, name) works for all valid names
@pytest.mark.parametrize('test_window_name,length', 
    [(name, size) for name, size in zip(window_names.keys(), [1,51,52])])
def test_create_window(test_window_name, length):
    create_window(length, name=test_window_name)



#test that create_window(N, name) is indeed equivalent to the direct call window_name
@pytest.mark.parametrize("name,param",
    [('blackman', {'alpha':2}),
      ('kaiser',   {'beta':8.6}),
      ('gaussian', {'alpha':2.5}),
      ('chebwin',  {'attenuation': 50}),
      ('flattop',  {'mode':'symmetric'}),
      ('tukey',    {'r': 0.5}),
      ('poisson',      {'alpha': 2}),
      ('poisson_hanning', {'alpha': 2}),
      ('cauchy',   {'alpha': 3})])


def test_check_window_switch(name, param):
    f = eval('window_'+name)
    w1 = f(64, **param)
    w2 = create_window(64, name, **param)
    for x,y in zip(w1,w2):
        assert x==y

def test_create_window_others():
    try:
        create_window(11, "hamming", dummy=1)
        assert False
    except ValueError:
        assert True

    try:
        create_window(11, "kaiser", beta=0.5)
        create_window(11, "kaiser", dummy=1)
        assert False
    except ValueError:
        assert True



def test_bartlett():
    """unit and functional test window_bartlett"""
    vec7 = array([ 0.,  0.33333333,  0.66666667,  1.,  0.66666667,  0.33333333,  0.])
    vec8 = array([ 0.,  0.28571429,  0.57142857,  0.85714286,  0.85714286,  0.57142857,  0.28571429,  0.])
    for x, y in zip(window_bartlett(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_bartlett(7), vec7):
        assert_almost_equal(x, y)


def test_kaiser():
    """unit and functional test window_kaiser"""
    vec8 = array([ 0.00133251,  0.09113651,  0.45964377,  0.92046158,  0.92046158,  0.45964377,  0.09113651,  0.00133251])
    vec7 = array([ 0.00133251,  0.13040195,  0.63041193,  1.        ,  0.63041193,     0.13040195,  0.00133251])
    for x, y in zip(window_kaiser(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_kaiser(7), vec7):
        assert_almost_equal(x, y)

    window_kaiser(8, method='other')
    window_kaiser(8, method='numpy')

    assert window_kaiser(1, method='numpy') == np.array([1])


def test_blackman():
    """unit and functional test window_bartlett"""
    vec8 = array([ -1.38777878e-17,   9.04534244e-02,   4.59182958e-01,  9.20363618e-01,   9.20363618e-01,   4.59182958e-01, 9.04534244e-02,  -1.38777878e-17])
    vec7 = array([ -1.38777878e-17,   1.30000000e-01,   6.30000000e-01, 1.00000000e+00,   6.30000000e-01,   1.30000000e-01, -1.38777878e-17])
    for x, y in zip(window_blackman(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_blackman(7), vec7):
        assert_almost_equal(x, y)
    for x, y in zip(window_blackman(7, alpha=0.16), vec7):
        assert_almost_equal(x, y)

    assert window_blackman(1, alpha=0.16) == np.array([1])


def test_hann():
    """unit and functional test window_bartlett"""
    vec7 = array([ 0.  ,  0.25,  0.75,  1.  ,  0.75,  0.25,  0.  ])
    vec8 = array([ 0.        ,  0.1882551 ,  0.61126047,  0.95048443,  0.95048443, 0.61126047,  0.1882551 ,  0.        ])
    for x, y in zip(window_hann(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_hann(7), vec7):
        assert_almost_equal(x, y)


def test_hammming():
    """unit and functional test window_hamming"""
    vec8 = array([ 0.08      ,  0.25319469,  0.64235963,  0.95444568,  0.95444568,   0.64235963,  0.25319469,  0.08      ])
    vec7 = array([ 0.08,  0.31,  0.77,  1.  ,  0.77,  0.31,  0.08])
    for x, y in zip(window_hamming(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_hamming(7), vec7):
        assert_almost_equal(x, y)


def test_chebwin():
    """unit and functional test chebwin"""
    vec7 = array([ 0.1116911 ,  0.41962999,  0.81377359,  1.        ,  0.81377359,   0.41962999,  0.1116911 ])
    vec8 = array([ 0.09455132,  0.34937508,  0.71822375,  1.        ,  1.        ,   0.71822375,  0.34937508,  0.09455132])
    for x, y in zip(window_checbwin(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_chebwin(7), vec7):
        assert_almost_equal(x, y)


def test_gaussian():
    """unit and functional test gaussian"""
    vec8 = array([ 0.09139376,  0.29502266,  0.64438872,  0.9523448 ,  0.9523448 ,  0.64438872,  0.29502266,  0.09139376])
    vec7 = array([ 0.1006689 ,  0.36044779,  0.77483743,  1.        ,  0.77483743,  0.36044779,  0.1006689 ])
    for x, y in zip(window_gaussian(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_gaussian(7), vec7):
        assert_almost_equal(x, y)


def test_cauchy():
    window_cauchy(64)


def test_cosine():
    window_cosine(64)
    assert window_cosine(1) == np.array([1.])


def test_riemann():
    window_riemann(64)


def test_lanczos():
    window_lanczos(64)
    assert window_lanczos(1) == np.array([1.])


def test_poisson():
    window_poisson(64)
    #assert window_poisson(1) == np.array([1.])


def test_poisson_hanning():
    window_poisson_hanning(64)


def test_bartlett_hann():
    vec7 = array([ 0.  ,  0.27,  0.73,  1.  ,  0.73,  0.27,  0.  ])
    vec8 = array([ 0., 0.2116453, 0.60170081, 0.92808246, 0.92808246, 0.60170081, 0.2116453, 0.])
    for x, y in zip(window_bartlett_hann(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_bartlett_hann(7), vec7):
        assert_almost_equal(x, y)

    assert window_bartlett_hann(1) == np.array([1.])

def test_window_visu():
    window_visu(64, 'hamming')

def test_enbw():
    N = 64
    w = create_window(N, 'rectangle')
    assert enbw(w) == 1.

def test_window_parzen():
    vec7  = array([ 0.0058309 ,  0.1574344 ,  0.65014577,  1.        ,  0.65014577, 0.1574344 ,  0.0058309 ])
    vec8 = array([ 0.00390625,  0.10546875,  0.47265625,  0.91796875,  0.91796875, 0.47265625,  0.10546875,  0.00390625])
    for x, y in zip(window_parzen(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_parzen(7), vec7):
        assert_almost_equal(x, y)


def test_bohman():
    vec8 = array([  3.89804309e-17,   7.07247468e-02,   4.37484012e-01, 9.10368513e-01,   9.10368513e-01,   4.37484012e-01, 7.07247468e-02,   3.89804309e-17])
    vec7 = array([  3.89804309e-17,   1.08997781e-01,   6.08997781e-01, 1.00000000e+00,   6.08997781e-01,   1.08997781e-01, 3.89804309e-17])
    for x, y in zip(window_bohman(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_bohman(7), vec7):
        assert_almost_equal(x, y)


def test_chebwin():
    vec8 = array([ 0.09455132,  0.34937508,  0.71822375,  1.,  1.,0.71822375,  0.34937508,  0.09455132])
    vec7 = array([ 0.1116911 ,  0.41962999,  0.81377359,  1.,  0.81377359, 0.41962999,  0.1116911 ])
    for x, y in zip(window_chebwin(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_chebwin(7), vec7):
        assert_almost_equal(x, y)


def test_nuttall():
    window_nuttall(64)
    assert window_nuttall(1) == np.array([1.])


def test_blackman_nuttall():
    vec7 = array([  3.62800000e-04,   6.13345000e-02,   5.29229800e-01,   1.00000000e+00,   5.29229800e-01,   6.13345000e-02,  3.62800000e-04])
    vec8 = array([  3.62800000e-04,   3.77757690e-02,   3.42727620e-01, 8.91851861e-01,   8.91851861e-01,   3.42727620e-01, 3.77757690e-02,   3.62800000e-04])
    for x, y in zip(window_blackman_nuttall(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_blackman_nuttall(7), vec7):
        assert_almost_equal(x, y)


def test_blackman_harris():
    vec7 = array([  6.00000000e-05,   5.56450000e-02,   5.20575000e-01, 1.00000000e+00,   5.20575000e-01,   5.56450000e-02, 6.00000000e-05])
    vec8 = array([  6.00000000e-05,   3.33917235e-02,   3.32833504e-01, 8.89369772e-01,   8.89369772e-01,   3.32833504e-01, 3.33917235e-02,   6.00000000e-05])
    for x, y in zip(window_blackman_harris(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_blackman_harris(7), vec7):
        assert_almost_equal(x, y)


def test_flattop():
    vec8 = array([  9.05133399e-04,  -2.64123651e-02,  -5.55579501e-02, 4.43549557e-01,   1.00000000e+00,   4.43549557e-01, -5.55579501e-02,  -2.64123651e-02])
    for x, y in zip(window_flattop(8, 'periodic', precision='octave'), vec8):
        assert_almost_equal(x, y)

    assert window_flattop(1) == np.array([1.])
    window_flattop(1, "periodic") 


def test_tukey():
    vec8 = array([ 0.        ,  0.61126047,  1.        ,  1.        ,  1.        ,     1.        ,  0.61126047,  0.        ])
    vec7 = array([ 0.  ,  0.75,  1.  ,  1.  ,  1.  ,  0.75,  0.  ])

    for x, y in zip(window_tukey(8), vec8):
        assert_almost_equal(x, y)
    for x, y in zip(window_tukey(7), vec7):
        assert_almost_equal(x, y)

    window_tukey(64, r=0)
    window_tukey(64, r=1)
    assert window_tukey(1, r=1) == np.array([1.])


def test_window_rietz():
    window_riesz(64)
