import os
import sys
from distutils.util import convert_path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Check Python version
if sys.version_info < (3, 6):
    raise Exception('openpnm requires Python 3.6 or greater to run')

sys.path.append(os.getcwd())

main_ = {}
ver_path = convert_path('openpnm/__init__.py')
with open(ver_path) as f:
    for line in f:
        if line.startswith('__version__'):
            exec(line, main_)

setup(
    name='OpenPNM',
    description = 'A framework for conducting pore network modeling simulations of multiphase transport in porous materials',
    version=main_['__version__'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    packages=[
        'openpnm',
        'openpnm.core',
        'openpnm.network',
        'openpnm.geometry',
        'openpnm.phases',
        'openpnm.physics',
        'openpnm.utils',
        'openpnm.models',
        'openpnm.algorithms',
        'openpnm.topotools'
    ],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-image',
        'pandas',
        'sympy',
        'numba',
        'networkx',
        'h5py',
        'porespy',
        'transforms3d',
        'flatdict',
        'scikit-umfpack'],

    author='OpenPNM Team',
    author_email='jgostick@uwaterloo.ca',
    download_url='https://github.com/pmeal/OpenPNM/',
    url='http://openpnm.org'
)
