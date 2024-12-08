import glob
import os
import sys

pj = os.path.join

from distutils.core import Extension

from setuptools import find_packages, setup

_MAJOR = 0
_MINOR = 9
_MICRO = 0
version = "%d.%d.%d" % (_MAJOR, _MINOR, _MICRO)
release = "%d.%d" % (_MAJOR, _MINOR)


with open("README.rst") as f:
    readme = f.read()


setup(
    name="spectrum",
    version=version,
    description="Spectrum Analysis Tools",
    long_description=readme,
    author="Thomas Cokelaer",
    author_email="cokelaer@gmail.com",
    url="http://github.com/cokelaer/spectrum",
    license="new BSD",
    ext_modules=[
        Extension(
            "spectrum.mydpss",
            [
                "src/cpp/mydpss.c",
            ],
            export_symbols=["multitap"],
        )
    ],
    packages=find_packages("src"),
    package_dir={"": "src"},
    # Dependencies
    install_requires=open("requirements.txt").read(),
    extras_require={
        "plot": ["matplotlib"],
        "testing": [
            "pytest",
            "pytest-cov",
            "pytest-xdist",
            "pytest-mock",
            "pytest-timeout",
            "pytest-runner",
            "coveralls",
        ],
        "doc": ["sphinx", "sphinx_rtd_theme"],
    },
    package_data={
        "spectrum.data": ["*"],
    },
    platforms=["Linux"],
    classifiers=[
        "Development Status :: 1 - Planning",
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
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
    ],
)
