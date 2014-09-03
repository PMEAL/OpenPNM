#!/usr/bin/env python

import os
import sys
sys.path.append(os.getcwd())
import OpenPNM

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='OpenPNM',
    packages=['OpenPNM',
              'OpenPNM.Base',
              'OpenPNM.Network',
              'OpenPNM.Geometry',
              'OpenPNM.Geometry.models',
              'OpenPNM.Phases',
              'OpenPNM.Phases.models',
              'OpenPNM.Physics',
              'OpenPNM.Physics.models',
              'OpenPNM.Utilities',
              'OpenPNM.Algorithms',
              'OpenPNM.Postprocessing',
              'OpenPNM.Postprocessing.Export'],
    version=OpenPNM.__version__,
    description="A framework for conducting pore network modeling simulations of multiphase transport in porous materials.",
    author='OpenPNM Team',
    author_email='jeff.gostick@mcgill.ca',
    download_url='https://github.com/pmeal/OpenPNM/tarball/v1.0alpha4',
    url='https://github.com/pmeal/OpenPNM',

)
