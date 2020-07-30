## OpenPNM: Gallery of Examples

This page links to examples in the github repo at [github.com/PMEAL/OpenPNM/examples/notebooks](https://www.github.com/PMEAL/OpenPNM/examples/notebooks)


[//]: # (This line and the one below are not rendered in the final file, so basically act as comments)
[//]: # (It's possible to insert images into the cells using <img src="../docs/static/logo.png" width="100" align="left"> )


### [Tutorials](/examples/notebooks/tutorials)

| Link | Description |
|:-----|:------------|
| [Tutorial 1 - Basic](/examples/notebooks/tutorials/intro_to_openpnm_basic.ipynb) | An overview of OpenPNM in terms of basic manual calculations |
| [Tutorial 2 - Intermediate](/examples/notebooks/tutorials/intro_to_openpnm_intermediate.ipynb) | A repeat of Tutorial 1, but using the features of OpenPNM correctly |
| [Concise but Detailed Overview](/examples/notebooks/tutorials/concise_overview_of_openpnm.ipynb)  | This notebook goes over most of the key aspects OpenPNM with minimal description, favoring links to other relevant notebooks and resources instead for further reading |
| [Storage of Network Data and Topology](/examples/notebooks/tutorials/data_and_topology_storage.ipynb) | Provides a explanation of how numerical data is stored in OpenPNM, including the storage of the topological information using sparse adjacency matrices |
| [Querying the Topology to Find Neighbors](/examples/notebooks/tutorials/finding_neighbor_pores_and_throats.ipynb) | Illustrates how to find neighboring pores and throats using topological and spatial information |
| [Using and Creating Labels](/examples/notebooks/tutorials/using_and_creating_labels.ipynb) | Labels are used to mark pores and throats for easier selection elsewhere in the simulation, such as applying boundary conditions.  This tutorial illustrates how to use the labels that are automatically included on generated networks and how to add new user-defined labels. |
| [Defining Multiple Subdomains](/examples/notebooks/tutorials/defining_multiple_subdomains.ipynb) | One of the "power" features of OpenPNM is the ability to define multiple subdomains, enabling the simulation of layered or hierarchical materials with different pore sizes.  This tutorial provides an explanation of this feature. |
| [The Workspace and Projects](/examples/notebooks/tutorials/overview_of_workspace_and_projects.ipynb) | OpenPNM simulations called Projects and several Projects can be open within a single Workspace.  This tutorial illustrates the relationship between these to constructs and their features. |

### [Network](/examples/notebooks/networks)

#### [Generation](/examples/notebooks/networks/generation)

| Link | Description |
|:-----|:------------|
| [Cubic Lattice](/examples/notebooks/networks/generation/cubic_networks.ipynb) | The Cubic lattice is the classic pore network topology, and provides an excellent starting point for any investigation.  The ability to alter the coodination number up to 26, combined with deleting random pores and throats creates an even more realistic topology. |
| [Dual Cubic Lattice](/examples/notebooks/networks/generation/dual_cubic_lattices.ipynb) | The simultaneous simulation transport in the void and solid, and their interactions, can be done with two interpenetrating cubic networks as illustrated in this tutorial |
| [Cubic Template](/examples/notebooks/networks/generation/cubic_templates.ipynb) | Cubic lattices do not have to be contrained to cubic shaped domains.  This tutorial explains how arbitrary shaped domains (e.g. a spherical particle) can be created. |
| [Delaunay and Voronoi Tessellation](/examples/notebooks/networks/generation/random_networks_based_on_delaunay_and_voronoi_tessellations.ipynb) | If truly random topology is necessary then this can be accomplished by performaing a Delaunay tessellation on randomly (or not so randomly) distributed points |

#### [Manipulation](/examples/notebooks/networks/manipulation)

| Link | Description |
|:-----|:------------|
| [Adding Boundary Pores](/examples/notebooks/networks/manipulation/adding_boundary_pores.ipynb) | Boundary pores are useful when conducting transport simulations, but they are not added to generated OpenPNM network automatically.  This tutorial explains their importance as well as how to add them to the network. |
| [Manually Adding and Removing Pores and Throats](/examples/notebooks/networks/manipulation/manually_adding_pores_and_throats.ipynb) | It's possible to alter a pore network by manually adding/removing pores and/or throats.  This tutorial desribes the available tools in OpenPNM and explains what occurs behind the scenes. |
| [Stitching and Merging Networks](/examples/notebooks/networks/manipulation/stitching_and_merging_networks_together.ipynb) | Given two or more networks, it is necessary to somehow join them into a single domain to perform simulations.  This tutorial explains the difference between merging and stitching, and gives some examples on how to do both. |
| [Joining Pore Network and Continuum Domains](/examples/notebooks/networks/manipulation/coupling_continuum_regions_with_pore_networks.ipynb) | This is a variation on stitching networks together, where one has much smaller spacing than the other |

#### [Extraction](/examples/notebooks/networks/extraction)

| Link | Description |
|:-----|:------------|
| [Pore-scale Imaging and Modeling](/examples/notebooks/networks/extraction/Pore-scale Imaging and Modeling-MGambier.ipynb) | `Done` |
| [Benthiemer ICL Benchmark](/examples/notebooks/networks/extraction/benthiemer_ICL_benchmark.ipynb) | `Done` |
| [Doddington ICL Benchmark](/examples/notebooks/networks/extraction/doddington_ICL_benchmark.ipynb) | `Done` |
| [Berea ICL Benchmark](/examples/notebooks/networks/extraction/berea_ICL_benchmark.ipynb) | `Done` |
| [Predicting Permeability of Berea](/examples/notebooks/networks/extraction/predicting_effective_permeability_of_berea.ipynb) | `Done` |
| [Working with Extracted Networks](/examples/notebooks/networks/extraction/working_with_extracted_networks.ipynb) |`Done` |

### [Geometry](/examples/notebooks/geometry)

| Link | Description |
|:-----|:------------|
| [Basic Stick and Ball](/examples/notebooks/geometry/stick_and_ball.ipynb) |  |
| [Defining Continuum Regions](/examples/notebooks/geometry/defining_continuum_regions.ipynb) |  |
| [Overview of Shape Factors](/examples/notebooks/geometry/overview_of_shape_factors.ipynb) | `Zohaib and/or Amin` |
| [Adjusting Pore Size Distributions](/examples/notebooks/geometry/adjusting_pore_size_distributions.ipynb) |  |

### [Materials](/examples/notebooks/materials)

| Link | Description |
|:-----|:------------|
| [Fibrous Media with Voronoi Tessellations](/examples/notebooks/materials/fibrous_media_based_on_voronoi_tessellation.ipynb) | `Done` |

### [Thermophysical Properties](/examples/notebooks/phases)

| Link | Description |
|:-----|:------------|
| [Creating a Custom Phase](/examples/notebooks/phases/creating_a_custom_phase.ipynb) | `Done` |
| [Working with Mixtures](/examples/notebooks/phases/working_with_mixtures.ipynb) | `Jeff` |

### [Simulations](/examples/notebooks/algorithms)

#### [Percolation](/examples/notebooks/algorithms/percolation)

| Link | Description |
|:-----|:------------|
| [Ordinary Percolation](/examples/notebooks/algorithms/percolation/A_ordinary_percolation.ipynb) |  |
| [Invasion Percolation](/examples/notebooks/algorithms/percolation/B_invasion_percolation.ipynb) |  |
| [Mixed Invasion Percolation](/examples/notebooks/algorithms/percolation/C_mixed_invasion_percolation.ipynb) |  |
| [Meniscus Model Comparison](/examples/notebooks/algorithms/percolation/D_meniscus_model_comparison.ipynb) |  |
| [Mercury Intrusion Porosimetry](/examples/notebooks/algorithms/percolation/capillary_pressure_curves.ipynb) |  |

#### [Single Phase Transport](/examples/notebooks/algorithms/single_phase)

| Link | Description |
|:-----|:------------|
| [Basic Fickian Diffusion, Tortuosity, and Formation Factor](/examples/notebooks/algorithms/single_phase/fickian_diffusion_and_tortuosity.ipynb) | `Jeff` |
| [Permeability Tensor](/examples/notebooks/algorithms/single_phase/stokes_flow_and_permeability_tensor.ipynb) | `Jeff` |
| [Deep Dive into Conductance Models](/examples/notebooks/algorithms/single_phase/deep_dive_into_conductance_models.ipynb) | `Amin and Zohaib` |
| [Diffusion with Concentration Dependent Diffusivity](/examples/notebooks/algorithms/single_phase/diffusion_with_concentration_dependent_diffusivity.ipynb) | `Amin` |
| [1D Heat Transfer](/examples/notebooks/algorithms/single_phase/one_dimensional_heat_transfer.ipynb) | `Remove?` |

#### [Multiphase Transport](/examples/notebooks/algorithms/multiphase)

| Link | Description |
|:-----|:------------|
| [Relative Diffusivity](/examples/notebooks/algorithms/multiphase/relative_diffusivity.ipynb) | `Done` |
| [Relative Permeability in 2D](/examples/notebooks/algorithms/multiphase/relative_permeability_2D.ipynb) | `Done` |

#### [Reactive Transport](/examples/notebooks/algorithms/reactive)

| Link | Description |
|:-----|:------------|
| [Diffusion with Source and Sink Terms](/examples/notebooks/algorithms/reactive/diffusion_with_source_and_sink_terms.ipynb) | OpenPNM is capable of simulating chemical reactions in pores by adding source and sink terms. This example shows how to add source and sink terms to a steady state fickian diffusion simulation. |
| [Heat Transfer with Source Terms](/examples/notebooks/algorithms/reactive/one_dimensional_continuum_heat_transfer_with_source_term.ipynb) | `Done` |

#### [Transient Transport](/examples/notebooks/algorithms/transient)

| Link | Description |
|:-----|:------------|
| [Transient Fickian Diffusion](/examples/notebooks/algorithms/transient/transient_fickian_diffusion.ipynb) | The package OpenPNM allows for the simulation of many transport phenomena in porous media such as Stokes flow, Fickian diffusion, advection-diffusion, transport of charged species, etc. Transient and steady-state simulations are both supported. An example of a transient Fickian diffusion simulation through a Cubic pore network is shown here. |
| [Transient Fickian Diffusion with Reaction](/examples/notebooks/algorithms/transient/transient_fickian_diffusion_with_reaction.ipynb) | OpenPNM supports adding reaction terms to both steady state and transient simulations. OpenPNM already includes many different source term models that can be added to simulate a reaction. In this example, we show how to add a powerlaw source term model to a transient fickian diffusion simulation. |
| [Transient Advection-Diffusion](/examples/notebooks/algorithms/transient/transient_advection_diffusion.ipynb) | `Stephen` |

#### [Multiphysics](/examples/notebooks/algorithms/multiphysics)

| Link | Description |
|:-----|:------------|
| [Advection-Diffusion](/examples/notebooks/algorithms/multiphysics/advection_diffusion.ipynb) | `Amin` |
| [Nernst-Planck-Poisson](/examples/notebooks/algorithms/multiphysics/nernst_planck_poisson.ipynb) | `Mehrez` |

#### [Solvers and Settings](/examples/notebooks/algorithms/general)

| Link | Description |
|:-----|:------------|
| [Understanding Basic Solver Settings](/examples/notebooks/algorithms/general/understanding_basic_transport_settings.ipynb) | `Mehrez` |
| [Understanding Reactive Solver Settings](/examples/notebooks/algorithms/general/understanding_reactive_transport_settings.ipynb) | `Mehrez` |
| [Understanding Transient Solver Settings](/examples/notebooks/algorithms/general/understanding_transient_transport_settings.ipynb) | `Mehrez` |
| [Comparison of Available Matrix Solvers](/examples/notebooks/algorithms/general/available_matrix_solvers.ipynb) | `Amin` |
| [Overview of Inheritance](/examples/notebooks/algorithms/general/overview_of_inheritance.ipynb) | `Amin` |

### [Import, Export and Visualization](/examples/notebooks/io)

| Link | Description |
|:-----|:------------|
| [Quick Plotting Networks](/examples/notebooks/io/quick_plotting_networks.ipynb) |  |
| [Rendering in Paraview](/examples/notebooks/io/rendering_in_paraview.ipynb) |  |
| [Loading a Statoil File and Calculating K](/examples/notebooks/io/loading_statoil_finding_permeability.ipynb) |  |
