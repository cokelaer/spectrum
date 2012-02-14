#from spectrum import datasets
from spectrum.datasets import *

def test_marple_data():
    d = marple_data
    assert len(d) == 64


def test_timeseries():
    data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)
    ts = TimeSeries(data, sampling=1)
    #assert ts.N == 1024
    #assert ts.sampling == 1024
    ts.plot()


def test_data_cosine():
    data = data_cosine(N=1024, A=0.1, sampling=1024, freq=200)
    
