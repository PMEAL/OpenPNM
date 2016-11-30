#!/usr/bin/env python3

import os
import sys
from distutils.util import convert_path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Check Python version
if sys.version_info < (3, 4):
    raise Exception('OpenPNM requires Python 3.4 or greater to run')

sys.path.append(os.getcwd())

main_ = {}
ver_path = convert_path('OpenPNM/__init__.py')
with open(ver_path) as f:
    for line in f:
        if line.startswith('__version__'):
            exec(line, main_)

setup(
    name='OpenPNM',
    description = 'A framework for conducting pore network modeling simulations of multiphase transport in porous materials.',
    version=main_['__version__'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    packages=[
        'OpenPNM',
        'OpenPNM.Base',
        'OpenPNM.Network',
        'OpenPNM.Network.models',
        'OpenPNM.Geometry',
        'OpenPNM.Geometry.models',
        'OpenPNM.Phases',
        'OpenPNM.Phases.models',
        'OpenPNM.Physics',
        'OpenPNM.Physics.models',
        'OpenPNM.Utilities',
        'OpenPNM.Algorithms',
        'OpenPNM.Postprocessing'
    ],
    install_requires=[
        'numpy',
        'scipy>=0.14.0',
        'matplotlib',
        'scikit-image',
        'transforms3d',
        'dill',
        'pandas',
        'pyyaml'
    ],
    author='OpenPNM Team',
    author_email='jeff.gostick@mcgill.ca',
    download_url='https://github.com/pmeal/OpenPNM/',
    url='https://github.com/pmeal/OpenPNM'
)
