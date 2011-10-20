Criteria for Parametric methods
===================================



yule walker 
----------------


.. plot::
    :width: 80%
    :include-source:

    from spectrum import *
    from pylab import *

    order = arange(1, 25)
    rho = [aryule(marple_data, i, norm='biased')[1] for i in order]
    plot(order, AIC(len(marple_data), rho, order), label='AIC')
