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
print data_files

setup(
    name=name,
    version=version,
    description="spectrum analysis tools",
    author="Thomas Cokelaer",
    author_email="cokelaer@gmail.com",
    url='',
    license='LGPL',


    packages = find_packages('src'),
    package_dir={ '' : 'src' },

    # Dependencies
    install_requires = ['matplotlib', 'numpy', 'scipy'],
    data_files = data_files,

)

