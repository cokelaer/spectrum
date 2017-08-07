import os, sys, glob
pj = os.path.join

from setuptools import setup, find_packages
from distutils.core import Extension




_MAJOR               = 0
_MINOR               = 6
_MICRO               = 6
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)


data_files = [ (pj('share', 'data'), glob.glob(os.path.join('share','data', '*.dat')))]

with open('README.rst') as f:
    readme = f.read()


setup(
    name="spectrum",
    version=version,
    description="Spectrum Analysis Tools",
    long_description=readme,
    author="Thomas Cokelaer",
    author_email="cokelaer@gmail.com",
    url='http://github.com/cokelaer/spectrum',
    license='new BSD',

    ext_modules=[Extension('spectrum.mydpss', ['src/cpp/mydpss.c', ],
        export_symbols=['multitap'])],

    packages = find_packages('src'),
    package_dir={ '' : 'src' },

    # Dependencies
    install_requires=['numpy', 'scipy', 'easydev'],
    extras_require={
        'plot': ['matplotlib']
    },
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
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        ]   
)

