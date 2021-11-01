r"""

**openpnm.algorithms**

----

The ``algorithms`` module contains classes for conducting transport simulations
on pore networks.

"""

from .GenericAlgorithm import GenericAlgorithm, SettingsGenericAlgorithm
from .GenericTransport import GenericTransport, SettingsGenericTransport

from .ReactiveTransport import ReactiveTransport, SettingsReactiveTransport
from .TransientReactiveTransport import TransientReactiveTransport

from .StokesFlow import StokesFlow
from .NonNewtonianStokesFlow import NonNewtonianStokesFlow

from .FickianDiffusion import FickianDiffusion
from .TransientFickianDiffusion import TransientFickianDiffusion

from .AdvectionDiffusion import AdvectionDiffusion
from .TransientAdvectionDiffusion import TransientAdvectionDiffusion

from .FourierConduction import FourierConduction
from .OhmicConduction import OhmicConduction

from .OrdinaryPercolation import OrdinaryPercolation
from .InvasionPercolation import InvasionPercolation
from .MixedInvasionPercolation import MixedInvasionPercolation
from .MixedInvasionPercolationCoop import MixedInvasionPercolationCoop

from .Porosimetry import Porosimetry

from .IonicConduction import IonicConduction
from .TransientIonicConduction import TransientIonicConduction

from .NernstPlanck import NernstPlanck
from .TransientNernstPlanck import TransientNernstPlanck

from .NernstPlanckMultiphysicsSolver import NernstPlanckMultiphysicsSolver
from .TransientNernstPlanckMultiphysicsSolver import (
    TransientNernstPlanckMultiphysicsSolver
)

from . import metrics
