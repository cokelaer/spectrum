Installation
=================

Using pip 
--------------

**Spectrum** is available on `PYPi <http://pypi.python.org/pypi>`_, so you should be able to type::

    pip install spectrum


Since **spectrum** depends on other python packages such as Numpy, Matplotlib and Scipy they will be installed automatically (if not already installed).

You can also install them yourself by typing::

    pip numpy matplotlib scipy

**Spectrum** source code is available on Github https://github.com/cokelaer/spectrum


Conda installation
---------------------

**Spectrum** is now available on CONDA. For Linux and MAC users, if you prefer to use conda, please use::

    conda config --add channels conda-forge
    conda install spectrum


From source and notes for developers
-----------------------------------------

Developpers who want to get the source code can clone the repository::

    git clone git@github.com:cokelaer/spectrum.git
    cd spectrum
    python setup.py install


Then, you can test the library using **pytest** or compile the documentation
with Sphinx. To do so, install sphinx and other dependencies::

    pip install --file requirements-dev.txt



