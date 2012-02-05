import os, sys, glob
pj = os.path.join

from setuptools import setup, find_packages

name='spectrum'
os.sys.path.append((pj(os.path.curdir, 'src', name)))
try:
    import spectrum
    version = spectrum.__version__
except:
    version = 0.9

data_files = [ (pj('share', 'data'), glob.glob(os.path.join('share','data', '*.dat')))]

setup(
    name=name,
    version=version,
    description="Spectrum Analysis Tools",
    long_description="""This package provides functions/classes to estimate Power Spectral Densities using methods based on Fourier transform, Parametric methods or eigenvalues analysis. The Fourier methods are based upon correlogram, periodogram and Welch estimates. Standard tapering windows (Hann, Hamming, Blackman) and more exotic ones are available (DPSS, Taylor, ...). The parametric methods are based on Yule-Walker, BURG, MA and ARMA, covariance and modified covariance methods. Non-parametric methods based on eigen analysis (e.g., MUSIC) and minimum variance analysis are also implemented""",
    author="Thomas Cokelaer",
    author_email="cokelaer@gmail.com",
    url='http://www.assembla.com/spaces/PySpectrum/wiki',
    license='LGPL',

    packages = find_packages('src'),
    package_dir={ '' : 'src' },

    # Dependencies
    install_requires = ['matplotlib', 'numpy', 'scipy'],
    data_files = data_files,
    platforms=["Linux"],
    classifiers=["Development Status :: 1 - Planning",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Financial and Insurance Industry",
        "Intended Audience :: Information Technology",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Telecommunications Industry",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering",
        ]   
)

