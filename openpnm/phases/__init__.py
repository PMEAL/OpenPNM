r"""

**openpnm.phases**

----

This module contains the GenericPhase class, plus several subclasses which
are preconfigured to have the properties of specific fluids

----

**The GenericPhase Class**

The ``GenericPhase`` class is a direct child of the ``Base`` class, so contains
the usual methods such as find pore indices based on labels.  It does, however,
also inherit from ``ModelsMixin`` so has methods for add, removing and
regenerating models.

----

**Library of Preconfigured Phase Classes**

OpenPNM include a few Phase subclasses that contain a suite of pre-configured
models that predict the thermophysical properties of certain common phases.

+-------------+---------------------------------------------------------------+
| Class       | Comments                                                      |
+=============+===============================================================+
| Water       | Most models include the impact of salinity                    |
+-------------+---------------------------------------------------------------+
| Air         | A mixture of O2 and N2, but no humidity                       |
+-------------+---------------------------------------------------------------+
| Mercury     | Useful for porosimetry simulations, assumed theta is 140      |
+-------------+---------------------------------------------------------------+

----

**Customizing a GenericPhase Instance**

The ``GenericPhase`` class has no pore-scale models attached, so is a blank
slate for creating custom Phases.  The following code snippet illustrates how
to do this to create an *oil* phase with a temperature dependent viscosity
model based on a simple 2nd order polynomial:

>>> import openpnm as op
>>> import numpy as np
>>> pn = op.network.Cubic([5, 5, 5])
>>> oil = op.phases.GenericPhase(network=pn)

Now add the pore-scale model for viscosity from the ``models`` module:

>>> oil.add_model(propname='pore.viscosity',
...               model=op.models.misc.polynomial,
...               a=[50000, 1, -.1],
...               prop='pore.temperature')

Upon adding the model, values are immediately calculated at the Phase's
current temperature:

>>> np.round(oil['pore.viscosity'][0])
41418.0

If the temperature is changed, the model can be regenerated to update the
viscosity values:

>>> oil['pore.temperature'] = 355
>>> np.round(oil['pore.viscosity'][0])
41418.0

Note that the oil viscosity has NOT changed!  To propigate the new temperature
to all the other calculated pore-scale properties, call the
``regenerate_models`` function:

>>> oil.regenerate_models()
>>> np.round(oil['pore.viscosity'][0])
37752.0

"""
from .GenericPhase import GenericPhase
from .Air import Air
from .Water import Water
from .Mercury import Mercury
from .MultiPhase import MultiPhase
