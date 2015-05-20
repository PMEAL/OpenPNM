r"""
###############################################################################
:mod:`OpenPNM.Phases` -- Phase Property Estimation Methods
###############################################################################

Contents
--------
GenericPhase: The basic class which defines how a Phase is instantiated.  It
also has a few specific method for querying the health of mixtures or physics
objects.

Subclasses: OpenPNM includes a few pre-written subclasses that describe the
most commonly used materials, like Air, Water and Mercury.  Creating a custom
Phase subclass simply requires placing a file in the Phases directory and it
will be automatically loaded.

Classes
-------

.. autoclass:: GenericPhase
   :members:

.. autoclass:: Air
   :members:

.. autoclass:: Water
   :members:

.. autoclass:: Mercury
   :members:

"""
from .__GenericPhase__ import GenericPhase
from .__Air__ import Air
from .__Water__ import Water
from .__Mercury__ import Mercury
from .__TestPhase__ import TestPhase
from . import models
