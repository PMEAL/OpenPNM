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

    def set_value_BC(self, pores=[], values=[], mode='overwrite', force=False):
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
            value is 'merge'. For definition of various modes, see the
            docstring for ``set_BC``.
        force : bool, optional
            If ``True`` then the ``'mode'`` is applied to all other bctypes as
            well. The default is ``False``.

        """
        self.set_BC(pores=pores, bctype='value', bcvalues=values,
                    mode=mode, force=force)

    def set_rate_BC(self, pores=[], rates=None, total_rate=None, mode='overwrite',
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
            value is 'merge'. For definition of various modes, see the
            docstring for ``set_BC``.
        force : bool, optional
            If ``True`` then the ``'mode'`` is applied to all other bctypes as
            well. The default is ``False``.

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
    def set_BC(self, pores, bctype, bcvalues=None, mode='replace',
               force=False):
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
            'overwrite'  (default) Adds supplied boundary conditions to the
                         existing conditions, including overwriting any
                         existing conditions that may be present. This is
                         equivalent to calling ``'remove'`` on the given
                         locations followed by ``'add'``. If ``force=True``
                         this also removes any BCs of other types in the
                         given locations.
            'add'        Adds the supplied boundary conditions to the
                         existing conditions but does *not* overwrite any
                         conditions that are already present. If ``force=True``
                         this will remove values from  any locations where
                         other BC types are present.
            'remove'     Removes boundary conditions from the specified
                         locations. if ``force=True`` this also removes
                         any BCs of the other types from the specified
                         locations.
            'clear'      Removes all boundary conditions from the object of
                         the of the specified type from all locations. If
                         ``force=True`` this clears all BCs of the other types
                         as well.
            ============ =====================================================

        force : bool, optional
            If ``True`` then the ``'mode'`` is applied to all other bctypes as
            well. The default is ``False``.

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm.

        """
        if not isinstance(bctype, str):
            raise Exception('bctype must be a single string')
        bc_types = list(self['pore.bc'].keys())
        other_types = np.setdiff1d(bc_types, bctype).tolist()
        mode = self._parse_mode(
            mode,
            allowed=['overwrite', 'add', 'remove', 'clear'],
            single=True
        )
        pores = self._parse_indices(pores)

        values = np.array(bcvalues)
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of values must match the number of locations')

        mask = np.ones_like(pores, dtype=bool)  # Indices of pores to keep
        if mode == 'add':
            # Remove indices that are already present for given bc type
            mask[np.isfinite(self[f'pore.bc.{bctype}'][pores])] = False
            if force:  # Set locs on other bcs to nan
                for item in other_types:
                    self[f"pore.bc.{item}"][pores[mask]] = np.nan
            # Now remove indices that are present for other BCs
            for item in other_types:
                mask[np.isfinite(self[f'pore.bc.{item}'][pores])] = False
            self[f"pore.bc.{bctype}"][pores[mask]] = values
        elif mode == 'overwrite':
            if force:  # Set locs on other bcs to nan
                for item in other_types:
                    self[f"pore.bc.{item}"][pores[mask]] = np.nan
            # Now remove indices that are present for other BCs
            for item in other_types:
                mask[np.isfinite(self[f'pore.bc.{item}'][pores])] = False
            self[f"pore.bc.{bctype}"][pores[mask]] = values
        elif mode == 'remove':
            if force:  # Set locs on other bcs to nan
                for item in other_types:
                    self[f"pore.bc.{item}"][pores[mask]] = np.nan
            # Now remove indices that are present for other BCs
            for item in other_types:
                mask[np.isfinite(self[f'pore.bc.{item}'][pores])] = False
            self[f"pore.bc.{bctype}"][pores[mask]] = np.nan
        elif mode == 'clear':
            self[f"pore.bc.{bctype}"] = np.nan
            if force:  # Set locs on other bcs to nan
                for item in other_types:
                    self[f"pore.bc.{item}"] = np.nan


if __name__ == "__main__":
    import openpnm as op
    from openpnm.models import collections
    pn = op.network.Cubic(shape=[3, 3, 1], spacing=1e-4)
    pn.add_model_collection(collections.geometry.cones_and_cylinders())
    pn.regenerate_models()
    air = op.phase.Air(network=pn, name="air")
    air.add_model_collection(collections.physics.standard())
    air.regenerate_models()
    # check mode add, with and without force
    fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
    fd['pore.bc.rate'][1] = 1.0
    fd['pore.bc.value'][0] = 1.0
    fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add')
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 2
    assert fd['pore.bc.value'][mask].sum() == 4.0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 1
    fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='add', force=True)
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 3
    assert fd['pore.bc.value'][mask].sum() == 7.0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 0
    # check mode overwrite, with and without force
    fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
    fd['pore.bc.rate'][1] = 1.0
    fd['pore.bc.value'][0] = 1.0
    fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='overwrite')
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 2
    assert fd['pore.bc.value'][mask].sum() == 6.0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 1
    fd.set_value_BC(pores=[0, 1, 2], values=3.0, mode='overwrite', force=True)
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 3
    assert fd['pore.bc.value'][mask].sum() == 9.0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 0
    # check mode remove, with and without force
    fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
    fd['pore.bc.rate'][1] = 1.0
    fd['pore.bc.value'][0] = 1.0
    fd.set_value_BC(pores=[0, 1], mode='remove')
    mask = np.isfinite(fd['pore.bc.value'])
    assert mask.sum() == 0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 1
    fd.set_value_BC(pores=[0, 1], mode='remove', force=True)
    assert np.isfinite(fd['pore.bc.value']).sum() == 0
    assert np.isfinite(fd['pore.bc.value']).sum() == 0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 0
    # check mode clear, with and without force
    fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
    fd['pore.bc.rate'][1] = 1.0
    fd['pore.bc.value'][0] = 1.0
    fd.set_value_BC(mode='clear')
    assert np.isfinite(fd['pore.bc.value']).sum() == 0
    assert np.isfinite(fd['pore.bc.rate']).sum() == 1
    fd.set_value_BC(mode='clear', force=True)
    assert np.isfinite(fd['pore.bc.rate']).sum() == 0












