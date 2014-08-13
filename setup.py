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
    version=OpenPNM.__version__,
    description="A framework for conducting pore network modeling simulations of multiphase transport in porous materials.",
    author='OpenPNM Team',
    author_email='jeff.gostick@mcgill.ca',
    url='http://openpnm.org/',
    license='MIT',

)
