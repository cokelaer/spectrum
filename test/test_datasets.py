from spectrum import datasets
from spectrum.datasets import *

def test_cos():
    d = data_cosine()

def test_marple_data():
    d = marple_data
    assert len(d) == 64
