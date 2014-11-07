#!/usr/bin/env python

import os
import sys
sys.path.append(os.getcwd())

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='OpenPNM',
    packages=['OpenPNM',
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
              'OpenPNM.Postprocessing'],
    version='1.1-beta',
    description="A framework for conducting pore network modeling simulations of multiphase transport in porous materials.",
    author='OpenPNM Team',
    author_email='jeff.gostick@mcgill.ca',
    download_url='https://github.com/pmeal/OpenPNM/',
    url='https://github.com/pmeal/OpenPNM',
)
