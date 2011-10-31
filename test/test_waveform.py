from spectrum import *
from spectrum import waveform
from waveform import *
from numpy import linspace


def test_morlet():
    x = morlet(0,1,100)
    try:
        x = morlet(0,1,-100)
        assert False
    except:
        assert True
    

def test_chirp():
    x = chirp(linspace(0,1,1000))
    x = chirp(linspace(0,1,1000), form='linear')
    x = chirp(linspace(0,1,1000), form='quadratic')
    x = chirp(linspace(0,1,1000), form='logarithmic')
    try:
        x = chirp(linspace(0,1,1000), form='dummy')
        assert False
    except:
        assert True


def test_mexican():
    mexican(0,1,10)
    try:
        mexican(0,1,0)
        assert False
    except:
        assert True

def test_meyeraux():
    meyeraux(10)
