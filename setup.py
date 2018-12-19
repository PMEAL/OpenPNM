import os
import sys
from distutils.util import convert_path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

sys.path.append(os.getcwd())

main_ = {}
ver_path = convert_path('openpnm/__init__.py')
with open(ver_path) as f:
    for line in f:
        if line.startswith('__version__'):
            exec(line, main_)

setup(
    name='openpnm',
    description = 'A framework for conducting pore network modeling simulations of multiphase transport in porous materials',
    version=main_['__version__'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
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
        'openpnm.io',
        'openpnm.models',
        'openpnm.models.misc',
        'openpnm.models.geometry',
        'openpnm.models.phases',
        'openpnm.models.physics',
        'openpnm.algorithms',
        'openpnm.topotools',
        'openpnm.materials',
    ],
    install_requires=[
        'numpy>=1.15',
        'scipy>=1.1',
        'scikit-image>=0.14',
        'networkx>=2',
        'h5py>=2.8',
        'sympy',
        'matplotlib',
        'pandas',
        'numba',
        'porespy',
        'transforms3d',
        'flatdict',
        'gitpython',
        'jsonschema',
        'unyt'
        ],
    author='OpenPNM Team',
    author_email='jgostick@uwaterloo.ca',
    download_url='https://github.com/pmeal/OpenPNM/',
    url='http://openpnm.org'
)
