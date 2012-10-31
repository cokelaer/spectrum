Criteria for Parametric methods
===================================

In order to estimate the order of a parametric model, one chose a PSD method such as the :func:`~spectrum.yulewalker.aryule` function. This function (when given an order) returns a list of AR parameters. The order selected may not be optimal (too low or too high). One tricky question is then to find a criteria to select this order in an optiaml way. Criteria are available and the following example illustrate their usage. 


Example 1
----------
Let us consider a data set (the Marple data already used earlier). We use the aryule function to estimate the AR parameter. This function also returns a parameter called `rho`. This parameter together with the length of the data and the selected order can be used by criteria functions such as the :func:`~spectrum.criteria.AIC` function to figure out the optimal order.

.. plot::
    :width: 80%
    :include-source:

    from spectrum import *
    from pylab import *

    order = arange(1, 25)
    rho = [aryule(marple_data, i, norm='biased')[1] for i in order]
    plot(order, AIC(len(marple_data), rho, order), label='AIC')


The optimal order corresponds to the minimal of the plotted function.

Example 2
-----------
We can look at another example that was look at earlier with a AR(4):

.. plot::
    :width: 80%
    :include-source:

    import scipy.signal
    from spectrum import *
    from pylab import *

    # Define AR filter coefficients and some data accordingly
    a = [1, -2.2137, 2.9403, -2.1697, 0.9606];
    x = scipy.signal.lfilter([1], a, randn(1,256))

    # study different order
    order = arange(1, 25)
    rho = [aryule(x[0], i, norm='biased')[1] for i in order]
    plot(order, AIC(len(x[0]), rho, order), label='AIC')


Here, is appears that an order of 4 (at least) should be used, which correspond indeed to the original choice.
