Concepts
===========

Frequency range
-------------------

Frequency range for the PSD estimate can one of 'onesided', 'twosided',
'centered'. The default is 'onesided' for real-valued signals and 'twosided' for
complex signals.

More details follows this convention:

* 'onesided': returns the one-sided PSD estimate of a real-valued input signal,
  x. If NFFT is even, PSD has length NFFT/2 + 1 over the interval [0,pi]. If NFFT is odd, the length of
  PSD is (NFFT+1)/2 and the interval is [0, pi)
* 'twosided' returns the two-sided PSD estimate for either real or complex
  values. PSD has the length NFFT and is computed over the interval [0,2pi).
* 'centered' returns the centered two-sided PSD estimate for either real or complex values. PSD has length NFFT and is computed over (-pi,pi] for even length NFFT and (-pi,pi] for odd length NFFT.


Padding
----------
NFFT is used to compute the fft with NFFT points. If the data is shorter, the
data is padded with zeros so that the frequency signal has a length equal to
NFFT. The PSD's length is not equal to NFFT (see above). 


PSD length for real data case
------------------------------


N = 8; FFT gives X0, X1,X2,X3, XN/2, X3,X2,X1 --> PSD length=2+6=5=NFFT/2+1
N = 9; FFT gives X0, X1,X2,X3, X4,X4, X3,X2,X1 --> PSD length=1+3+1=5=(NFFT+1)/2

if one sets the psd, the NFFT is set manually to 2+ (N - 2) * 2
if one sets the psd, we assume that the original data is even.



