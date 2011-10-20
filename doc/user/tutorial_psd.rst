Spectrum object
=================

The normal usage to estimate a PSD is to use of the PSD estimate starting with the letter p such as parma, pminvar, pburg, (exception: use Periodgram instead of pperiodogram). All the PSD estimates are based upon the Spectrum class, which can be used as illustrated in this short example.

Let us create a Spectrum instance::

    from spectrum import *
    p = Spectrum(data_cosine(), sampling=1024)

This instance stored the original data size and sampling::

    p.N
    p.sampling


However, it contains no information about the PSD estimation method. If you type::

    p.psd

it should return a warning message telling that psd is not yet computed. 


You can either compute the PSD independantly, and then populate the psd attribute::

    psd = speriodogram(p.data)

Or set the PSD method::

    p.method = minvar

and call the function with the proper optional arguments::

    p(15, NFFT=4096)

Now, psd should contain the PSD vector. 

Again, the best way in such case is to use::

    p =  pminvar(data, 15)
    p()

that takes care of setting the relevant attributes. 

