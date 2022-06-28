import logging
import numpy as np
from openpnm.utils import Docorator


__all__ = [
    'BCsMixin',
]


docstr = Docorator()
logger = logging.getLogger(__name__)


class BCsMixin:
    """
    Mixin class to add boundary condition functionality to algorithms.
    """

    def set_value_BC(self, pores, values, mode='merge'):
        r"""
        Applies constant value boundary conditons to the specified pores.

        These are sometimes referred to as Dirichlet conditions.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        values : float or array_like
            The value to apply in each pore. If a scalar is supplied
            it is assigne to all locations, and if a vector is applied is
            must be the same size as the indices given in ``pores``.
        mode : str, optional
            Controls how the boundary conditions are applied. The default
            value is 'merge'. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'merge'      (default) Adds supplied boundary conditions to the
                         existing conditions, including replacing any
                         existing conditions they may be present. This is
                         equivalent to calling ``'remove'`` on the given
                         locations followed by ``'add'``.
            'overwrite'  Deletes all boundary conditions of the given type
                         then adds the specified new ones.  This is equivalent
                         to called ``'clear'`` followed by ``'add'``.
            'add'        Adds the supplied boundary conditions to the
                         existing conditions but does *not* overwrite any
                         conditions they are already present.
            'remove'     Removes boundary conditions from the specified
                         locations.
            'clear'      Removes all boundary conditions from the object.
            ============ =====================================================

            If BCs of the another type already exist in the given locations,
            those values are kept.

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quantity'] = 'pore.pressure'``.

        """
        self._set_BC(pores=pores, bctype='value', bcvalues=values, mode=mode)

    def set_rate_BC(self, pores, rates=None, total_rate=None, mode='merge',
                    **kwargs):
        r"""
        Apply constant rate boundary conditons to the specified locations.

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        rates : float or array_like, optional
            The rates to apply in each pore. If a scalar is supplied that
            rate is assigned to all locations, and if a vector is supplied
            it must be the same size as the indices given in ``pores``.
        total_rate : float, optional
            The total rate supplied to all pores. The rate supplied by
            this argument is divided evenly among all pores. A scalar must
            be supplied! Total_rate cannot be specified if rate is
            specified.
        mode : str, optional
            Controls how the boundary conditions are applied. The default
            value is 'merge'. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'merge'      (default) Adds supplied boundary conditions to the
                         existing conditions, including replacing any
                         existing conditions they may be present. This is
                         equivalent to calling ``'remove'`` on the given
                         locations followed by ``'add'``.
            'overwrite'  Deletes all boundary conditions of the given type
                         then adds the specified new ones.  This is equivalent
                         to called ``'clear'`` followed by ``'add'``.
            'add'        Adds the supplied boundary conditions to the
                         existing conditions but does *not* overwrite any
                         conditions they are already present.
            'remove'     Removes boundary conditions from the specified
                         locations.
            'clear'      Removes all boundary conditions from the object.
            ============ =====================================================

        Notes
        -----
        The definition of ``quantity`` is specified in the algorithm's
        ``settings``, e.g. ``alg.settings['quantity'] = 'pore.pressure'``.

        """
        # handle total_rate feature
        if total_rate is not None:
            if not np.isscalar(total_rate):
                raise Exception('total_rate argument accepts scalar only!')
            if rates is not None:
                raise Exception('Cannot specify both arguments: rate and '
                                + 'total_rate')
            pores = self._parse_indices(pores)
            rates = total_rate/pores.size
        self._set_BC(pores=pores, bctype='rate', bcvalues=rates, mode=mode)

    @docstr.get_sections(base='GenericTransport._set_BC',
                         sections=['Parameters', 'Notes'])
    def _set_BC(self, pores, bctype, bcvalues=None, mode='merge'):
        r"""
        This private method is called by public facing BC methods, to
        apply boundary conditions to specified pores

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied
        bctype : str
            Specifies the type or the name of boundary condition to apply.
            The types can be one one of the following:

            ============ =====================================================
            bctype       meaning
            ============ =====================================================
            'value'      Specify the value of the quantity in each pore
            'rate'       Specify the flow rate into each pore
            ============ =====================================================

        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.
            If a single value is given, it's assumed to apply to all
            locations unless the 'total_rate' bc.type is supplied whereby
            a single value corresponds to a total rate to be divded evenly
            among all pores. Otherwise, different values can be applied to
            all pores in the form of an array of the same length as
            ``pores``.
        mode : str, optional
            Controls how the boundary conditions are applied. The default
            value is 'merge'. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'merge'      (default) Adds supplied boundary conditions to the
                         existing conditions, including overwriting any
                         existing conditions they may be present. This is
                         equivalent to calling ``'remove'`` on the given
                         locations followed by ``'add'``.
            'overwrite'  Deletes all boundary conditions of the given type
                         then adds the specified new ones.  This is equivalent
                         to called ``'clear'`` followed by ``'add'``.
            'add'        Adds the supplied boundary conditions to the
                         existing conditions but does *not* overwrite any
                         conditions they are already present.
            'remove'     Removes boundary conditions from the specified
                         locations.
            'clear'      Removes all boundary conditions from the object.
            ============ =====================================================

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm. Use ``remove_BCs`` to
        clear existing BCs before applying new ones.

        """
        # Hijack the parse_mode function to verify bctype argument
        bc_types = list(self['pore.bc'].keys())
        bctype = self._parse_mode(
            bctype,
            allowed=bc_types,
            single=True
        )
        othertypes = np.setdiff1d(bc_types, bctype).tolist()
        mode = self._parse_mode(
            mode,
            allowed=['merge', 'overwrite', 'add', 'remove', 'clear'],
            single=True
        )
        pores = self._parse_indices(pores)

        values = np.array(bcvalues)
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of values must match the number of locations')

        # Find locations where other bc types are defined
        other_bcs = np.zeros(self.Np, dtype=bool)
        for item in othertypes:
            other_bcs += np.isfinite(self[f"pore.bc.{item}"])
        other_inds = pores[other_bcs[pores]]

        # Find locations which are unique to the current bc type
        current_inds = np.setdiff1d(pores, other_inds)

        # Catch pores with existing BCs
        if mode == 'add':
            temp = np.setdiff1d(pores, current_inds)
            self[f"pore.bc.{bctype}"][temp] = np.nan
        elif mode == 'remove':
            self[f"pore.bc.{bctype}"][pores] = np.nan
        elif mode == 'clear':
            self[f"pore.bc.{bctype}"] = np.nan
        elif mode == 'merge':
            self[f"pore.bc.{bctype}"][current_inds] = values
        elif mode == 'overwrite':   # Remove existing BCs and write new ones
            self[f"pore.bc.{bctype}"] = np.nan
            self[f"pore.bc.{bctype}"][current_inds] = values

    def remove_BC(self, pores=None, bctype='all'):
        r"""
        Removes boundary conditions from the specified pores.

        Parameters
        ----------
        pores : array_like, optional
            The pores from which boundary conditions are to be removed. If
            no pores are specified, then BCs are removed from all pores. No
            error is thrown if the provided pores do not have any BCs assigned.
        bctype : str, or List[str]
            Specifies which type of boundary condition to remove. The
            default value is 'all'. Options are:

            ============ =====================================================
            bctype       meaning
            ============ =====================================================
            'all'        Removes all boundary conditions
            'value'      Removes only value conditions
            'rate'       Removes only rate conditions
            ============ =====================================================

        """
        if isinstance(bctype, str):
            bctype = [bctype]
        if 'all' in bctype:
            bctype = ['value', 'rate']
        if pores is None:
            pores = self.Ps
        if ('pore.bc.value' in self.keys()) and ('value' in bctype):
            self['pore.bc.value'][pores] = np.nan
        if ('pore.bc.rate' in self.keys()) and ('rate' in bctype):
            self['pore.bc.rate'][pores] = np.nan
