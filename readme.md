[![](https://badge.fury.io/py/openpnm.svg)](https://pypi.python.org/pypi/openpnm)  [![](https://travis-ci.org/PMEAL/OpenPNM.svg?branch=master)](https://travis-ci.org/PMEAL/OpenPNM)
[![](https://codecov.io/gh/PMEAL/OpenPNM/branch/master/graph/badge.svg)](https://codecov.io/gh/PMEAL/OpenPNM)
[![](https://img.shields.io/badge/ReadTheDocs-GO-blue.svg)](http://openpnm.readthedocs.io/en/master/)
[![](https://ci.appveyor.com/api/projects/status/7lmjmjbq09p8l5dn/branch/master?svg=true&passingText=windows%20-%20OK)](https://ci.appveyor.com/project/jgostick/openpnm/branch/master)

# Overview of OpenPNM

*OpenPNM* is an open source project aiming to provide porous media researchers with a comprehensive framework for performing pore network simulations on a wide range of materials.

## Installation and Requirements

OpenPNM can be installed from the Python Package index using:

```
pip install openpnm
```

Or the source code can be downloaded from [Github](https://github.com/pmeal/OpenPNM/) and installed by running:

```
pip install -e 'path/to/downloaded/files'
```

The advantage to installing from the source code is that you can edit the files and have access to your changes each time you import *OpenPNM*.

OpenPNM requires the *Scipy Stack* (Numpy, Scipy, Matplotlib, etc), which is most conveniently obtained by installing the [Anaconda Distribution](https://conda.io/docs/user-guide/install/download.html).

## Example Usage

The following code block illustrates how to use OpenPNM to perform a mercury intrusion porosimetry simulation:

``` python

import openpnm as op
pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
Hg = op.phases.Mercury(network=pn)
phys = op.physics.Standard(network=pn, phase=Hg, geometry=geo)
mip = op.algorithms.Porosimetry(network=pn)
mip.setup(phase=Hg)
mip.set_inlets(pores=pn.pores(['left', 'right', 'top', 'bottom']))
mip.run()

```

The network can be visualized in [Paraview](http://www.paraview.org) giving the following:

![](https://i.imgur.com/mSDrIBOm.png)

The drainage curve can be visualized with `MIP.plot_intrusion_curve()` giving something like this:

![](https://i.imgur.com/1C2uXt9m.png)

A collection of examples is available as a separate Github repository: [OpenPNM-Examples](https://www.github.com/PMEAL/OpenPNM-Examples)

## Release Management and Versioning

OpenPNM uses [Semantic Versioning](http://semver.org) (i.e. X.Y.Z) to label releases.  All major and minor versions (X.Y.z) are available on [PyPI](https://pypi.python.org/pypi), but bugfixe releases (x.y.Z) are not generally pushed unless the bug is important.

OpenPNM uses the [Github Flow](https://guides.github.com/introduction/flow/) system of Git branching, except instead of merging PRs into *master*, they are merged into a branch called *dev*. Any code added to *dev* is done via Pull Requests (PRs).  When new PRs are merged into the *dev* branch, they are *not* given a new version number. Once enough new features have been added, the *dev* branch is merged into the *master* branch, and the minor release number (x.Y.z) will be incremented. An exception to this rule are bugfixes which may be found on *master*.  In these cases a PR can be merged into *master* and the version number wil be incremented (x.y.Z) to indicate the fix.

OpenPNM depends on several other packages widely known as the [Scipy Stack](https://www.scipy.org/stackspec.html).  It is our policy to always support the latest version of all these packages and their dependencies.

The main developer for this project is Prof. Jeff Gostick (jgostick@gmail.com).

## Licence and Citation

OpenPNM is free to use and is offered under the permissive [MIT License](http://opensource.org/licenses/MIT)

If you do use OpenPNM in an academic work, the developers ask that you cite the following paper, which outlines the design principles and general uses of OpenPNM:

    Gostick et al. OpenPNM: A pore network modeling package. Computing in Science & Engineering. 18(4), p60-74.

A link to this article can be found [here](http://doi.org/10.1109/MCSE.2016.49).
