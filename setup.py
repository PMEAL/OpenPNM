import os
import sys
from distutils.util import convert_path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

sys.path.append(os.getcwd())

about = {}
ver_path = convert_path('openpnm/__version__.py')
with open(ver_path) as f:
    for line in f:
        if line.startswith('__version__'):
            exec(line, about)

# Read the contents of README file
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='openpnm',
    description='A framework for conducting pore network modeling simulations '
    + 'of multiphase transport in porous materials',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=about['__version__'],
    zip_safe=False,
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
        'openpnm.phases.mixtures',
        'openpnm.phases.mixtures.species.ions',
        'openpnm.phases.mixtures.species.liquids',
        'openpnm.phases.mixtures.species.gases',
        'openpnm.physics',
        'openpnm.utils',
        'openpnm.io',
        'openpnm.models',
        'openpnm.models.topology',
        'openpnm.models.misc',
        'openpnm.models.geometry',
        'openpnm.models.phases',
        'openpnm.models.physics',
        'openpnm.algorithms',
        'openpnm.algorithms.metrics',
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
        'transforms3d',
        'flatdict',
        'gitpython',
        'jsonschema',
        'unyt',
        'terminaltables',
        'docrep',
    ],
    author='OpenPNM Team',
    author_email='jgostick@uwaterloo.ca',
    download_url='https://github.com/pmeal/OpenPNM/',
    url='http://openpnm.org',
    project_urls={
        'Documentation': 'https://openpnm.readthedocs.io/en/master/',
        'Source': 'https://github.com/PMEAL/OpenPNM',
        'Tracker': 'https://github.com/PMEAL/OpenPNM/issues',
    },
)
