from spectrum.errors import *

def test_is_positive():
    assert True == is_positive_integer(1)

    try:
        is_positive_integer(-1)
        assert False
    except:
        assert True
    try:
        is_positive_integer(1.)
        assert False
    except:
        assert True

def test_errors():

    a = SpectrumError()
    print a

    a = SpectrumChoiceError("dummy", ['valid'])
    print a

    a = SpectrumPSDError()
    print a

    a = SpectrumModifiedError()
    print a

    a = SpectrumARMAError()
    print a

    a = SpectrumMAError()
    print a

    a = SpectrumARError()
    print a
    a = SpectrumOrder()
    print a
    a = SpectrumNPSD()
    print a
    

    try:
        raise SpectrumModifiedError
        assert False
    except:
        assert True
    try:
        raise SpectrumError
        assert False
    except:
        assert True
