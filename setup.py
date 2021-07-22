import os
import sys
import codecs
import os.path
from distutils.util import convert_path
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

sys.path.append(os.getcwd())
ver_path = convert_path('openpnm/__version__.py')


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


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
    version=get_version(ver_path),
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
        'chemicals',
        'docrep>=0.3',
        'flatdict',
        'gitpython',
        'h5py',
        'ipython',
        'jsonschema',
        'json-tricks',
        'matplotlib',
        'networkx',
        'numba',
        'numpy',
        'pandas',
        'pypardiso',
        'scikit-image',
        'scipy',
        'sympy',
        'terminaltables',
        'tqdm',
        'transforms3d',
    ],
    author='OpenPNM Team',
    author_email='jgostick@uwaterloo.ca',
    download_url='https://github.com/PMEAL/OpenPNM/',
    url='http://openpnm.org',
    project_urls={
        'Documentation': 'https://pmeal.github.io/OpenPNM',
        'Source': 'https://github.com/PMEAL/OpenPNM',
        'Tracker': 'https://github.com/PMEAL/OpenPNM/issues',
    },
)
