[![](https://github.com/PMEAL/OpenPNM/actions/workflows/nightly.yml/badge.svg)](https://github.com/PMEAL/OpenPNM/actions/workflows/nightly.yml)
[![](https://codecov.io/gh/PMEAL/OpenPNM/branch/dev/graph/badge.svg)](https://codecov.io/gh/PMEAL/OpenPNM)
[![](https://img.shields.io/badge/Documentation-Read-blue.svg)](https://pmeal.github.io/OpenPNM/)
[![](https://badge.fury.io/py/openpnm.svg)](https://pypi.python.org/pypi/openpnm)
[![](https://anaconda.org/conda-forge/openpnm/badges/installer/conda.svg)](https://anaconda.org/conda-forge/openpnm)

-----

> VERSION 3.0 of OpenPNM is now out. All the [examples on the website](https://openpnm.org/_examples/index.html) are now using the features and idioms of V3. For a description of the main changes please see our [recent blog post](http://pmeal.com/posts/2022-10-10-notebook-post/). 

# Overview of OpenPNM

OpenPNM is a comprehensive framework for performing pore network simulations of porous materials.

## More Information

For more details about the package can be found in the [online documentation](https://openpnm.org)

## Installation and Requirements

### [pip](https://pypi.org/project/openpnm/)
OpenPNM can be installed using `pip` by running the following command in a terminal:

```shell
pip install openpnm
```

### [conda-forge](https://anaconda.org/conda-forge/openpnm)
OpenPNM can also be installed from the [conda-forge](https://anaconda.org/conda-forge/openpnm) repository using:

```
conda install -c conda-forge openpnm
```

### For developers
For developers who intend to change the source code or contribute to OpenPNM, the source code can be downloaded from [Github](https://github.com/PMEAL/OpenPNM/) and installed by running:

```
pip install -e 'path/to/downloaded/files'
```

The advantage to installing from the source code is that you can edit the files and have access to your changes each time you import OpenPNM.

OpenPNM requires the Scipy Stack (Numpy, Scipy, Matplotlib, etc), which is most conveniently obtained by installing the [Anaconda Distribution](https://www.anaconda.com/download/).

## Asking Questions and Getting Help

Github now has a [Discussions](https://github.com/PMEAL/OpenPNM/discussions) function, which works similarly to [stack overflow](https://www.stackoverflow.com).  Please post your question in the [Q&A category](https://github.com/PMEAL/OpenPNM/discussions?discussions_q=category%3AQ%26A) so devs or users can provide answers, vote on accepted answers, improve on each other's answers, and generally discuss things. Most importantly, all answers are searchable so eventually, once enough questions have been posted and answered, you can find what you're looking for with a simple search.

## Contact

OpenPNM is developed by the Porous Materials Engineering and Analysis Lab [(PMEAL)](http://pmeal.com), in the [Department of Chemical Engineering](https://uwaterloo.ca/chemical-engineering/) at the [University of Waterloo](https://uwaterloo.ca/) in Waterloo, Ontario, Canada.

The lead developer for this project is Prof. Jeff Gostick (jgostick@gmail.com).

## Acknowledgements

OpenPNM is grateful to [CANARIE](https://canarie.ca) for their generous funding over the past few years.  We would also like to acknowledge the support of [NSERC of Canada](https://www.nserc-crsng.gc.ca/) for funding many of the student that have contributed to OpenPNM since it's inception in 2011.

## Citation

If you use OpenPNM in a publication, please cite the following paper:

> _Gostick et al._ "**OpenPNM: a pore network modeling package.**" Computing in Science & Engineering 18, no. 4 (2016): 60-74.
> [doi:10.1109/MCSE.2016.49](https://ieeexplore.ieee.org/document/7478437)

Also, we ask that you "star" :star: this repository so we can track the number of users who are interested in this project, which is helpful for securing future grant funding.
