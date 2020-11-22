# OpenPNM Tutorials

The following table lists the available tutorials in the approximate order they should be tackled, along with a brief description of each.

| Link | Description |
|:-----|:------------|
| [Concise Overview](concise_overview_of_openpnm.ipynb) | A whirlwind tour of OpenPNM many features and function. This run-through attempts to show all the core features of OpenPNM is a shortest possible way, but includes links to the relevant examples throughout. |
| [Learning OpenPNM - Introductory](intro_to_openpnm_basic.ipynb) | Learn how to compute permeability of a network *without* using most of OpenPNMs features. This provides a baseline understanding of the type calculations are done in a pore network model, without the added confusion learning OpenPNM. |
| [Learning OpenPNM - Intermediate](intro_to_openpnm_intermediate.ipynb) | A repeat of the basic introductory tutorial, but utilizing the features of OpenPNM.  The hope is that after completing the basic tutorial, it will be clear why OpenPNM works the way it does. |
| [Learning OpenPNM - Advanced](intro_to_openpnm_advanced.ipynb) | See the advanced usage of OpenPNM such as customizing the code and defining pore-scale models. |
| [Data Storage](data_and_topology_storage.ipynb) | Dig into the details of how data is stored by OpenPNM, including both properties like pore diameter, as well as network topology. |
| [Topological Lookups](finding_neighbor_pores_and_throats.ipynb) | It is quite often necessary to find pores and throats that are neighbors to a given set of pores, for instance to assign boundary conditions, or to trim pores from the network.  This tutorial oulines the various means of finding and selecting pores and throats. |
| [Labels](using_and_creating_labels.ipynb) | Labels are a highly useful feature of OpenPNM, and this tutorial outlines their creation and use. Applying labels to pores allows for easy look-up later. |
| [Pore Scale Models](pore_scale_models.ipynb) | Pore-scale models are the main way to customize OpenPNM, for instance by defining a new capillary pressure relationship or a specifying the properties of a special fluid. |
| [Workspace and Projects](overview_of_workspace_and_projects.ipynb) | The workspace and projects are the way to manage several simulations at the same time. |
| [Data Exchange](data_exchange_between_objects.ipynb) | To simplify the interaction with data, OpenPNM provides a few short-cuts for reading data, as explained in this tutorial. |
| [Multiple Subdomains](defining_multiple_subdomains.ipynb) | The ability to model multiscale materials such as layered structures of hierarchical media was built-in to the design of OpenPNM.  This accomplished by assigning the pores of each scale to a their own subdomain, as illustrated in this tutorial. |