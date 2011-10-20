variance outputs
================

.. plot::
    :width: 80%
    :include-source:

    from pylab import *
    import scipy.signal
    from spectrum import *

    # Define AR filter coefficients
    a = [1, -2.2137, 2.9403, -2.1697, 0.9606];

    true_variance = linspace(0.1, 1, 20)
    estimated_variance = []
    for tv in true_variance:
        x = scipy.signal.lfilter([1], a, tv**0.5 * randn(1,256))
        AR, v, k = arburg(x[0], 4)
        estimated_variance.append(v)

    plot(true_variance, estimated_variance, 'o')
    xlabel('true variance')
    ylabel('estimated variance')


