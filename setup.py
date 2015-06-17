#!/usr/bin/env python3

import os
import sys
import OpenPNM

sys.path.append(os.getcwd())

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name = 'OpenPNM',
    description = 'A framework for conducting pore network modeling simulations of multiphase transport in porous materials.',
    version = OpenPNM.__version__,
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Physics'
    ],
    packages = [
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
    install_requires = [
        'numpy',
        'scipy>=0.14.0',
        'matplotlib'
    ],
    author = 'OpenPNM Team',
    author_email = 'jeff.gostick@mcgill.ca',
    download_url = 'https://github.com/pmeal/OpenPNM/',
    url = 'https://github.com/pmeal/OpenPNM'
)
