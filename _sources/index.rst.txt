.. _front_page:

###############################
Pore-network modeling made easy
###############################

.. sidebar:: Highlight

    The image below was extracted using the SNOW algorithm from PoreSpy, then
    imported into OpenPNM for use, then exported to Paraview for visualization.

    .. image:: /../docs/static/images/extracted_berea.png
        :width: 600px
        :align: center

*OpenPNM* is an open source project to provide porous media researchers with a ready-made framework for performing a wide range of pore network simulations.

Highly customizable
###################

*OpenPNM* is highly customizable, meaning that you can customize it almost at every level. Usually, it's sufficient to use of our built-in classes, and customize it by adding your own models. However, if need be, you can write your own *Network*, *Geometry*, *Phase*, *Physics*, and *Algorithm* from scratch easily by using our generic classes.

Rich IO ecosystem
#################
If you already have a network model in a typical file format, but want to use OpenPNM to perform simulations, that's not a problem! OpenPNM supports a wide variety of file formats :mod:`~openpnm.io` to import from and export to such as :class:`~openpnm.io.Statoil`, :class:`~openpnm.io.NetworkX`, etc.

.. toctree::
    :maxdepth: 1

    installation.rst
    quick_start.rst
    userguide.rst
    modules/index.rst
    citation.rst
