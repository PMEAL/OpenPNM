[![](https://github.com/PMEAL/OpenPNM/workflows/Ubuntu/badge.svg)](https://github.com/PMEAL/OpenPNM/actions)
[![](https://github.com/PMEAL/OpenPNM/workflows/macOS/badge.svg)](https://github.com/PMEAL/OpenPNM/actions)
[![](https://github.com/PMEAL/OpenPNM/workflows/Windows/badge.svg)](https://github.com/PMEAL/OpenPNM/actions)
[![](https://github.com/PMEAL/OpenPNM/workflows/Examples/badge.svg)](https://github.com/PMEAL/OpenPNM/actions)<br/>
[![](https://codecov.io/gh/PMEAL/OpenPNM/branch/master/graph/badge.svg)](https://codecov.io/gh/PMEAL/OpenPNM)
[![](https://img.shields.io/badge/ReadTheDocs-GO-blue.svg)](http://openpnm.readthedocs.io/en/master/)
[![](https://badge.fury.io/py/openpnm.svg)](https://pypi.python.org/pypi/openpnm)
[![](https://anaconda.org/conda-forge/openpnm/badges/installer/conda.svg)](https://anaconda.org/conda-forge/openpnm)


-----

**Please cite as:**
> _Gostick et al._ "**OpenPNM: a pore network modeling package.**" Computing in Science & Engineering 18, no. 4 (2016): 60-74.
> [doi:10.1109/MCSE.2016.49](https://ieeexplore.ieee.org/document/7478437)

# Overview of OpenPNM

*OpenPNM* is a comprehensive framework for performing pore network simulations of porous materials.

For more details about the package can be found in the [on-line documentation](https://openpnm.readthedocs.io/en/master/index.html)

## Installation and Requirements

OpenPNM can be installed from the [Python Package Index](https://pypi.org/project/openpnm/) using:

```
pip install openpnm
```

or from [Anaconda Cloud](https://anaconda.org/conda-forge/openpnm) using:

```
conda install -c conda-forge openpnm
```

or the source code can be downloaded from [Github](https://github.com/pmeal/OpenPNM/) and installed by running:

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

The network can be visualized in [`ParaView`](http://www.paraview.org) giving the following:

<p align="center">
  <img src="https://user-images.githubusercontent.com/14086031/77927983-dc3dd000-7275-11ea-8700-c96c2d51aa1f.png" width="60%"></img>
</p>

The drainage curve can be visualized with `mip.plot_intrusion_curve()` giving something like this:

<p align="center">
  <img src="https://user-images.githubusercontent.com/14086031/77930201-96363b80-7278-11ea-95fd-4a55fb1d6148.png" width="60%"></img>
</p>

A collection of examples is available in the *examples* folder of this repository: [Examples](https://www.github.com/PMEAL/OpenPNM/tree/master/examples)

## Contact

OpenPNM is developed by the Porous Materials Engineering and Analysis Lab [(PMEAL)](http://pmeal.com), in the [Department of Chemical Engineering](https://uwaterloo.ca/chemical-engineering/) at the [University of Waterloo](https://uwaterloo.ca/) in Waterloo, Ontario, Canada.

The lead developer for this project is Prof. Jeff Gostick (jgostick@gmail.com).
