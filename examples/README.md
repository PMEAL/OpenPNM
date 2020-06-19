## OpenPNM: Gallery of Examples

This page links to examples in the github repo at [github.com/PMEAL/OpenPNM/examples/notebooks](https://www.github.com/PMEAL/OpenPNM/examples/notebooks)


[//]: # (This line and the one below are not rendered in the final file, so basically act as comments)
[//]: # (https://github.com/PMEAL/OpenPNM/blob/master/examples/XXXX.ipynb)


### Tutorials

| Link | Description |
|:---|:---|
| [Tutorial 1 - Basic](/examples/notebooks/intro_to_openpnm_basic.ipynb) | An overview of OpenPNM in terms of basic manual calculations |
| [Tutorial 2 - Intermediate](/examples/notebooks/tutorials/intro_to_openpnm_intermediate.ipynb) | A repeat of tutorial 1, but using the features of OpenPNM correctly |
| [Tutorial 3 - Advanced](/examples/notebooks/tutorials/intro_to_openpnm_advanced.ipynb) | A deeper dive into OpenPNM including defining custom pore-scale models and phases |
| [Concise but Detailed Overview](/examples/notebooks/tutorials/concise_overview_of_openpnm.ipynb)  | This notebook goes over most of the key aspects OpenPNM with minimal description, favoring links to other relevant notebooks and resources instead for further reading. |
| [Storage of Network Data and Topology](/examples/notebooks/tutorials/data_and_topology_storage.ipynb) | Provides a explanation of how data is stored in OpenPNM, including the storage of the topological information |
| [Querying the Topology to Find Neighbors](/examples/notebooks/tutorials/finding_neighbor_pores_and_throats.ipynb) | Illustrates how to find neighboring pores and throats using topological and spatial information. |
| [Using and Creating Labels](/examples/notebooks/tutorials/using_and_creating_labels.ipynb) |  |
| [Defining Multiple Subdomains](/examples/notebooks/tutorials/defining_multiple_subdomains.ipynb) |  |
| [The Workspace and Projects](/examples/notebooks/tutorials/overview_of_workspace_and_projects.ipynb) |  |

### Network

#### Generation

| Link | Description |
|:---|:---|
| [Cubic Lattice]() | `Done` |
| [Dual Cubic Lattice]() | `Done` |
| [Cubic Template]() | `Jeff` |
| [Delaunay and Voronoi Tessellation]() | `Jeff` |

#### Manipulation

| Link | Description |
|:---|:---|
| [Adding Boundary Pores](/examples/notebooks/networks/manipulation/) | `Done` |
| [Adding Pores and Throats](/examples/notebooks/networks/manipulation/) |  |
| [Merging Networks](/examples/notebooks/networks/manipulation/) |  |
| [Stitching Networks](/examples/notebooks/networks/manipulation/) |  |
| [Joining Pore Network and Continuum Domains](/examples/notebooks/networks/manipulation/) |  |

#### Extraction

| Link | Description |
|:---|:---|
| [Benthiemer ICL Benchmark](/examples/notebooks/networks/extraction/) | `Done` |
| [Doddington ICL Benchmark](/examples/notebooks/networks/extraction/) | `Done` |
| [Berea ICL Benchmark](/examples/notebooks/networks/extraction/) | `Done` |
| [Working with Extracted Networks](/examples/notebooks/networks/extraction/) |`Niloo` |

### Geometry Calculations

| Link | Description |
|:---|:---|
| [Basic Stick and Ball](/examples/notebooks/geometry/) |  |
| [Defining Continuum Regions](/examples/notebooks/geometry/) |  |
| [Overview of Shape Factors](/examples/notebooks/geometry/) | `Zohaib and/or Amin` |
| [Adjusting Pore Size Distributions](/examples/notebooks/geometry/) | `Jeff` |

### Predefined Materials

| Link | Description |
|:---|:---|
| [Fibrous Media with Voronoi Tessellations](/examples/notebooks/materials/) | `Done` |

### Thermophysical Properties

| Link | Description |
|:---|:---|
| [Creating a Custom Phase](/examples/notebooks/phases/) | `Done` |
| [Working with Mixtures](/examples/notebooks/phases/) | `Jeff` |

### Simulations

#### Percolation

| Link | Description |
|:---|:---|
| [Ordinary Percolation](/examples/notebooks/algorithms/percolation/) |  |
| [Invasion Percolation](/examples/notebooks/algorithms/percolation/) |  |
| [Mixed Invasion Percolation](/examples/notebooks/algorithms/percolation/) |  |
| [Meniscus Model Comparison](/examples/notebooks/algorithms/percolation/) |  |

#### Single Phase Transport

| Link | Description |
|:---|:---|
| [Basic Fickian Diffusion](/examples/notebooks/algorithms/single_phase/) | `Jeff` |
| [Permeability Tensor](/examples/notebooks/algorithms/single_phase/) | `Jeff` |
| [Deep Dive into Conductance Models](/examples/notebooks/algorithms/single_phase/) | `Amin and Zohaib` |
| [Diffusion with Concentration Dependent Diffusivity](/examples/notebooks/algorithms/single_phase/) | `Amin` |
| [1D Heat Transfer](/examples/notebooks/algorithms/single_phase/) | `Remove?` |

#### Multiphase Transport

| Link | Description |
|:---|:---|
| [Relative Diffusivity](/examples/notebooks/algorithms/multiphase/) | `Niloo` |
| [Relative Permeability in 2D](/examples/notebooks/algorithms/multiphase/) | `Niloo` |

#### Reactive Transport

| Link | Description |
|:---|:---|
| [Diffusion with Source and Sink Terms](/examples/notebooks/algorithms/reactive/) | `Mike` |
| [Heat Transfer with Source Terms](/examples/notebooks/algorithms/reactive/) | `Done` |

#### Transient Transport

| Link | Description |
|:---|:---|
| [Transient Fickian Diffusion](/examples/notebooks/algorithms/transient/) | `Mike` |
| [Transient Advection-Diffusion](/examples/notebooks/algorithms/transient/) | `Stephen` |

#### Multiphysics

| Link | Description |
|:---|:---|
| [Advection-Diffusion](/examples/notebooks/algorithms/multiphysics/) | `Amin` |
| [Nernst-Planck-Poisson](/examples/notebooks/algorithms/multiphysics/) | `Mehrez` |

#### Solvers and Settings

| Link | Description |
|:---|:---|
| [Overview of Basic Solver Settings](/examples/notebooks/algorithms/general/) | `Mehrez` |
| [Overview of Reactive Solver Settings](/examples/notebooks/algorithms/general/) | `Mehrez` |
| [Overview of Transient Solver Settings](/examples/notebooks/algorithms/general/) | `Mehrez` |
| [Overview of Available Matrix Solvers](/examples/notebooks/algorithms/general/) | `Amin` |

### Import, Export and Visualization

| Link | Description |
|:---|:---|
| [Quick Plotting Networks](/examples/notebooks/io/) | `Amin` |
| [Rendering in Paraview](/examples/notebooks/io/) | `Jeff` |
| [Loading a Statoil File and Calculating K](/examples/notebooks/io/) | `Done` |
