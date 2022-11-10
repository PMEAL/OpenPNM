import os
import sys
import codecs
import os.path
from distutils.util import convert_path
from setuptools import setup, find_packages

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
            ver = line.split(delim)[1].split(".")
            if "dev0" in ver:
                ver.remove("dev0")
            return ".".join(ver)
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
    packages=find_packages("."),
    install_requires=[
        'chemicals',
        'docrep',
        'h5py',
        'jsonschema',
        'matplotlib',
        'networkx',
        'numba',
        'numpy',
        'pandas',
        'pyamg',
        'pypardiso',
        'rich',
        'scikit-image',
        'scipy',
        'sympy',
        'thermo',
        'tqdm',
        'transforms3d',
    ],
    author='OpenPNM Team',
    author_email='jgostick@uwaterloo.ca',
    download_url='https://github.com/PMEAL/OpenPNM/',
    url='http://openpnm.org',
    project_urls={
        'Documentation': 'https://openpnm.org',
        'Source': 'https://github.com/PMEAL/OpenPNM',
        'Tracker': 'https://github.com/PMEAL/OpenPNM/issues',
    },
)
