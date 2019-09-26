from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging
import numpy as np
logger = logging.getLogger(__name__)


class FickianDiffusion(ReactiveTransport):
    r"""
    A class to simulate binary diffusion.

    Parameters
    ----------
    network : OpenPNM Network object
        The network on which this algorithm operates

    project : OpenPNM Project object
        Either a network or a project must be specified

    name : string, optional
        A unique name to give the object for easier identification.  If not
        given, one is generated.

    Notes
    -----
    Fickian diffusion in porous materials occurs in the void space, but
    becuase the diffusion is defined to pores it is impacted by the porosity
    and tortuosity of the network.  Thus the total diffusive flux through the
    network is reduced.  This class can be used to simualte diffusion-reaction
    in domains with arbitrarily complex boundary conditions, or it can be used
    to calculate the effective diffusivity of the network by applying
    controlled boundary conditions on opposing faces, calculate the diffusion
    rate, and inverting Fick's first law:

    .. math::

        D_{eff} = N_{A}*L/(A*\Delta C_{A})

    This class includes a method for calculating Deff automatically assuming
    appropriate boundary conditions were applied (``calc_eff_diffusivity``).
    The length and area of the domain should be supplied, but if they are
    not an attempt is made to calculate them.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration',
                   'conductance': 'throat.diffusive_conductance',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': ''},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase on which the algorithm is to be run.

        quantity : string
            (default is ``'pore.mole_fraction'``)  The name of the physical
            quantity to be calculated.

        conductance : string
            (default is ``'throat.diffusive_conductance'``) The name of the
            pore-scale transport conductance values.  These are typically
            calculated by a model attached to a *Physics* object associated
            with the given *Phase*.

        Notes
        -----
        Any additional arguments are added to the ``settings`` dictionary of
        the object.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        super().setup(**kwargs)

    def set_continuity_BC(self, ps1, ps2, K12=1.0, mode="merge"):
        r"""
        Apply continuity boundary conditon between two phases.

        This boundary condition enforces c[ps1] = c[ps2] * K12

        Parameters
        ----------
        ps1 : array_like
            The pore indices in phase 1 where the condition should be applied

        ps2 : array_like
            The pore indices in phase 2 where the condition should be applied

        K12 : scalar or array_like
            Partition coefficient; relates the concentrations at two phases.
            If a scalar is supplied it is assigne to all locations, and if a
            vector is supplied, it must be the same size as the indices given
            in ``ps1`` and ``ps2``.

        mode : string, optional
            Controls how the boundary conditions are applied.  Options are:

            - ``'merge'``: (Default) Adds supplied boundary conditions to
            already existing conditions

            - ``'overwrite'``: Deletes all boundary condition on object then
            adds the given ones
        """
        # Hijack the parse_mode function to verify bctype argument
        ps1 = self._parse_indices(ps1)
        ps2 = self._parse_indices(ps2)
        K12 = np.array(K12)
        if ps1.size != ps2.size:
            if ps1.size != 1 and ps2.size != 1:
                raise Exception("Inconsistent array length: ps1 and ps2")
        if ps1.size < ps2.size:
            ps1, ps2, K12 = ps2, ps1, 1/K12
        if ps2.size == 1:
            ps2 = ps2.repeat(ps1.size)

        # Store partition coefficient (K12) values
        if ('pore.bc_continuity' not in self.keys()) or (mode == 'overwrite'):
            self['pore.bc_continuity'] = np.empty((self.Np, 2)) * np.nan
        self['pore.bc_continuity'][ps1, 0] = ps2
        self['pore.bc_continuity'][ps1, 1] = K12

    def calc_effective_diffusivity(self, inlets=None, outlets=None,
                                   domain_area=None, domain_length=None):
        r"""
        This calculates the effective diffusivity in this linear transport
        algorithm.

        Parameters
        ----------
        inlets : array_like
            The pores where the inlet composition boundary conditions were
            applied.  If not given an attempt is made to infer them from the
            algorithm.

        outlets : array_like
            The pores where the outlet composition boundary conditions were
            applied.  If not given an attempt is made to infer them from the
            algorithm.

        domain_area : scalar, optional
            The area of the inlet (and outlet) boundary faces.  If not given
            then an attempt is made to estimate it, but it is usually
            underestimated.

        domain_length : scalar, optional
            The length of the domain between the inlet and outlet boundary
            faces.  If not given then an attempt is made to estimate it, but it
            is usually underestimated.

        Notes
        -----
        The area and length of the domain are found using the bounding box
        around the inlet and outlet pores which do not necessarily lie on the
        edge of the domain, resulting in underestimation of sizes.
        """
        return self._calc_eff_prop(inlets=inlets, outlets=outlets,
                                   domain_area=domain_area,
                                   domain_length=domain_length)
