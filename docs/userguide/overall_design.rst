.. _overall_design:

###############################################################################
Overall Design
###############################################################################

OpenPNM separates different types of properties between different objects.  There are 5 types: **Network**, **Geometry**, **Phase**, **Physics**, and **Algorithms**.  Each of these are described in more detail below, but their names hopefully indicate what sort of data or roles are assigned to each.

===============================================================================
Base
===============================================================================

The OpenPNM **Base** class, from which all other main OpenPNM objects descend, is a subclass of the Python Dictionary or ``dict``.

In addition to the methods included on every ``dict`` (e.g. ``pop``, ``keys``, etc), the OpenPNM **Base** class has quite a few additional methods for working specifically with OpenPNM data.
