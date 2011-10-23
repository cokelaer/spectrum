from spectrum import lms
from lms import LMS



def test_lms():
    U = 0.98
    x = [1,2,3,4,5]
    ar = [1,0.5,0.2]
    LMS(x, U, ar)
