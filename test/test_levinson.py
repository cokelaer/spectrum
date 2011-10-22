from spectrum import *
from spectrum.levinson import *
from nose.tools import assert_almost_equal
import numpy
from numpy.testing import assert_array_almost_equal


#LEVINSON
T0 = 3
T = numpy.array([-2+0.5j, .7-1j])


expected_P = 1.3221
expected_A = numpy.array([.86316+0.03158j, .34737+0.21053j])




def test_levinson():

    data = [T0]
    data.extend(T)
    A, P,k = LEVINSON(data)

    #use only 4 digits to compare float P
    assert_almost_equal(P, expected_P, 4)

    #compare element by element of vector A
    for a1, a2 in zip(A, expected_A):
        assert_almost_equal(a1, a2, 3)

def test_levinson_real():
    r = [5.0000, -1.5450, -3.9547, 3.9331, 1.4681, -4.7500]
    A, E, K = LEVINSON(r)
    assert_array_almost_equal(A, numpy.array([  6.14739427e-01,   9.89813712e-01,   4.20968656e-04,
         3.44472001e-03,  -7.70967347e-03]))
    assert_almost_equal(0.1791451516, E)


    A, E,k = LEVINSON([1, 0.5,0.3],1)
    for a1, a2 in zip(A, [-0.5]):
        assert_almost_equal(a1, a2)
    assert_almost_equal(E, 0.75)


    A, E, K = LEVINSON([1,0.5,0.1,0.05])
    assert_array_almost_equal(A, (array([-0.625,  0.275, -0.125])))
    assert_array_almost_equal(E, 0.708749)
    assert_array_almost_equal(K, array([-0.5  ,  0.2  , -0.125]))

