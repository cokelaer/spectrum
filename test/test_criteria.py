from spectrum.criteria import *
from spectrum import *



def test_aic():
    order = arange(1, 25)
    rho = [aryule(marple_data, i, norm='biased')[1] for i in order]
    AIC(len(marple_data), rho, order)

    AIC(len(marple_data), rho, order)
    KIC(len(marple_data), rho, order)
    AKICc(len(marple_data), rho, order)
    FPE(len(marple_data), rho, order)
    MDL(len(marple_data), rho, order)
    CAT(len(marple_data), rho, order)


    crit = Criteria("AIC", 20)
    try:
        crit = Criteria("dummy", 20)
        assert False
    except:
        assert True

    try:
        crit = Criteria("AIC", 20)
        crit.N = 0
        assert False
    except:
        assert True
