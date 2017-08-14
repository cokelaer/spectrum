Installation and source 
=========================

**Spectrum** is available on `PYPi <http://pypi.python.org/pypi>`_, so you should be able to type::

    pip install spectrum 


Since spectrum depends on other python packages such as Numpy, Matplotlib and Scipy they will be installed automatically (if not already installed). 

You can also install the dependencies yourself by typing::

    pip install numpy matplotlib scipy

**Spectrum** source code is available on Github https://github.com/cokelaer/spectrum




Conda installation
========================

Spectrum is not available on CONDA. However, you can install dependencies as
follows::

    conda install numpy matplotlib scipy 

I put a todo issue on the github account if you want to do the full portage on
conda-forge for instance (https://github.com/cokelaer/spectrum/issues/28)

For Linux and MAC users, if you prefer to use conda (avoiding 
compilation of dependencies), please use::

    conda config --add channels conda-forge 
    conda install spectrum



