r"""
###############################################################################
:mod:`OpenPNM.Physics` -- Pore Scale Physics Models
###############################################################################

Contents
--------
GenericPhysics: This base class is essentially an init statement that ensures
the Physics object registers itself with Phase and Network objects correctly.

Subclasses: OpenPNM includes one subclass called Standard which invokes the
typical pore scale physics models such as the Hagen-Poiseiulle model for
hydraulic conductance through a tube.  Customized Physics classes can be
created by adding a file to the Physics directory, which will be imported
automatically.

Classes
-------
.. autoclass:: GenericPhysics
   :members:

.. autoclass:: Standard
   :members:

"""

from .__GenericPhysics__ import GenericPhysics
from .__Standard__ import Standard
from .__TestPhysics__ import TestPhysics
from . import models
