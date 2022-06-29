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

    def set_value_BC(self, pores, values, mode='merge', force=False):
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
        self.set_BC(pores=pores, bctype='value', bcvalues=values,
                    mode=mode, force=force)

    def set_rate_BC(self, pores, rates=None, total_rate=None, mode='merge',
                    force=False):
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
        self.set_BC(pores=pores, bctype='rate', bcvalues=rates, mode=mode,
                    force=force)

    @docstr.get_sections(base='GenericTransport._set_BC',
                         sections=['Parameters', 'Notes'])
    def set_BC(self, pores, bctype, bcvalues=None, mode='replace', force=False):
        r"""
        The main method for setting and adjusting boundary conditions.

        This method is called by other more convenient wrapper functions like
        ``set_value_BC``.

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied
        bctype : str
            Specifies the type or the name of boundary condition to apply. This
            can be anything, but normal options are 'rate' and 'value'.
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
            'replace'    (default) Adds supplied boundary conditions to the
                         existing conditions, including overwriting any
                         existing conditions they may be present. This is
                         equivalent to calling ``'remove'`` on the given
                         locations followed by ``'add'``. If ``force=True``
                         this also overwrites any BCs of other types.
            'add'        Adds the supplied boundary conditions to the
                         existing conditions but does *not* overwrite any
                         conditions they are already present. If ``force=True``
                         this will overwrite any locations where other BC
                         types are present.
            'remove'     Removes boundary conditions from the specified
                         locations. if ``force=True`` this also removes
                         any BCs of the other types.
            'clear'      Removes all boundary conditions from the object of
                         the of the specified type. If ``force=True`` this
                         clears all BCs of the other types as well.
            ============ =====================================================

        force : bool, optional
            If ``True`` then the ``'mode'`` is applied to all other bctypes as
            well.

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm. Use ``remove_BCs`` to
        clear existing BCs before applying new ones.

        """
        if not isinstance(bctype, str):
            raise Exception('bctype must be a single string')
        bc_types = list(self['pore.bc'].keys())
        other_types = np.setdiff1d(bc_types, bctype).tolist()
        mode = self._parse_mode(
            mode,
            allowed=['replace', 'add', 'remove', 'clear'],
            single=True
        )
        pores = self._parse_indices(pores)

        values = np.array(bcvalues)
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of values must match the number of locations')

        # Find locations where other bc types are defined
        other_bcs = np.zeros(self.Np, dtype=bool)
        for item in other_types:
            other_bcs += np.isfinite(self[f"pore.bc.{item}"])
        other_inds = pores[other_bcs[pores]]

        # Find locations which are unique to the current bc type
        current_inds = np.setdiff1d(pores, other_inds)

        # Catch pores with existing BCs
        if mode == 'add':
            mask = np.ones_like(pores, dtype=bool)
            for item in other_types:
                mask[np.isfinite(self[f'pore.bc.{item}'][pores])] = False
            self[f"pore.bc.{bctype}"][pores[mask]] = values
            if force:
                for item in other_types:
                    self[f"pore.bc.{bctype}"][pores] = values
        elif mode == 'remove':
            self[f"pore.bc.{bctype}"][pores] = np.nan
        elif mode == 'clear':
            self[f"pore.bc.{bctype}"] = np.nan
        elif mode == 'merge':
            self[f"pore.bc.{bctype}"][current_inds] = values
        elif mode == 'overwrite':   # Remove existing BCs and write new ones
            self[f"pore.bc.{bctype}"] = np.nan
            self[f"pore.bc.{bctype}"][current_inds] = values


if __name__ == "__main__":
    import openpnm as op
    from openpnm.models import collections
    pn = op.network.Cubic(shape=[3, 3, 1], spacing=1e-4)
    pn.add_model_collection(collections.geometry.cones_and_cylinders())
    pn.regenerate_models()
    air = op.phase.Air(network=pn, name="air")
    air.add_model_collection(collections.physics.standard())
    air.regenerate_models()
    fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
    # check mode add, with and without force
    fd['pore.bc.rate'][1] = 1.0
    fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add')
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 2
    assert fd['pore.bc.value'][mask].sum() == 6.0
    fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add', force=True)
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 3
    assert fd['pore.bc.value'][mask].sum() == 9.0














