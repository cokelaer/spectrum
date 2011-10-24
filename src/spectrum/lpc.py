from numpy.fft import fft, ifft
from tools import nextpow2
from numpy import real
from spectrum.levinson import LEVINSON

__all__ = ['lpc']

def lpc(x, N=None):
    """Linear Predictor Coefficients.

    :param x:
    :param int N: default is length(X) - 1
    
    :Details:

    Finds the coefficients :math:`A=(1, a(2), \dots a(N+1))`, of an Nth order 
    forward linear predictor that predicts the current value value of the 
    real-valued time series x based on past samples:
    
    .. math:: \hat{x}(n) = -a(2)*x(n-1) - a(3)*x(n-2) - ... - a(N+1)*x(n-N)

    such that the sum of the squares of the errors

    .. math:: err(n) = X(n) - Xp(n)

    is minimized. This function  uses the Levinson-Durbin recursion to 
    solve the normal equations that arise from the least-squares formulation.  

    .. seealso:: :func:`levinson`, :func:`aryule`, :func:`prony`, :func:`stmcb`

    .. todo:: matrix case, references
    
    :Example:

    ::
    
        from scipy.signal import lfilter
        noise = randn(50000,1);  % Normalized white Gaussian noise
        x = filter([1], [1 1/2 1/3 1/4], noise)
        x = x[45904:50000]
        x.reshape(4096, 1)
        x = x[0]

    Compute the predictor coefficients, estimated signal, prediction error, and autocorrelation sequence of the prediction error:
   

    1.00000 + 0.00000i   0.51711 - 0.00000i   0.33908 - 0.00000i   0.24410 - 0.00000i

    ::
 
        a = lpc(x, 3)
        est_x = lfilter([0 -a(2:end)],1,x);    % Estimated signal
        e = x - est_x;                        % Prediction error
        [acs,lags] = xcorr(e,'coeff');   % ACS of prediction error

    """

    m = len(x)
    if N == None:
        N = m - 1 #default value if N is not provided
    elif N > m-1:
        #disp('Warning: zero-padding short input sequence')
        x.resize(N+1)
        #todo: check this zero-padding. 

    X = fft(x, 2**nextpow2(2.*len(x)-1))
    R = real(ifft(abs(X)**2))
    R = R/(m-1.) #Biased autocorrelation estimate
    a, e, ref = LEVINSON(R, N)
    return a, e
