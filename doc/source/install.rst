Installation and source 
=========================

**Spectrum** is available on `PYPi <http://pypi.python.org/pypi>`_, so you should be able to type::

    easy_install -U spectrum 


Since spectrum depends on other python packages such as Numpy, Matplotlib and Scipy they will be installed automatically (if not already installed). 

You can also install them yourself by typing::

    easy_install -U numpy
    easy_install -U matplotlib
    easy_install -U scipy


.. note:: the -U option updates the libraries that are already installed.


**Spectrum** source code is available on Github https://github.com/cokelaer/spectrum


Conda installation
========================

Spectrum is not available on CONDA. However, you can install dependencies as
follows::

    conda install numpy matplotlib scipy 

I put a todo issue on the github account if you want to do the full portage on
conda-forge for instance (https://github.com/cokelaer/spectrum/issues/28)


