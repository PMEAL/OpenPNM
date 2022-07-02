import logging
import numpy as np
from openpnm.core import ParserMixin, LabelMixin, Base2
from openpnm.utils import Docorator
import functools


__all__ = ['GenericAlgorithm']


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='GenericAlgorithmSettings', sections=docstr.all_sections)
@docstr.dedent
class GenericAlgorithmSettings:
    r"""

    Parameters
    ----------
    %(BaseSettings.parameters)s

    """


@docstr.get_sections(base='GenericAlgorithm', sections=['Parameters'])
@docstr.dedent
class GenericAlgorithm(ParserMixin, LabelMixin, Base2):
    r"""
    Generic class to define the foundation of Algorithms

    Parameters
    ----------
    %(Base.parameters)s

    """

    def __init__(self, network, name='alg_#', **kwargs):
        super().__init__(network=network, name=name, **kwargs)
        self.settings._update(GenericAlgorithmSettings())
        self['pore.all'] = np.ones(network.Np, dtype=bool)
        self['throat.all'] = np.ones(network.Nt, dtype=bool)

    # @functools.cached_property
    @property
    def iterative_props(self):
        r"""
        Finds and returns properties that need to be iterated while
        running the algorithm.
        """
        import networkx as nx
        phase = self.project[self.settings.phase]
        # Generate global dependency graph
        dg = nx.compose_all([x.models.dependency_graph(deep=True)
                             for x in [phase]])
        variable_props = self.settings["variable_props"].copy()
        variable_props.add(self.settings["quantity"])
        base = list(variable_props)
        # Find all props downstream that depend on base props
        dg = nx.DiGraph(nx.edge_dfs(dg, source=base))
        if len(dg.nodes) == 0:
            return []
        iterative_props = list(nx.dag.lexicographical_topological_sort(dg))
        # "variable_props" should be in the returned list but not "quantity"
        if self.settings.quantity in iterative_props:
            iterative_props.remove(self.settings["quantity"])
        return iterative_props

    def _update_iterative_props(self, iterative_props=None):
        """
        Regenerates phase, geometries, and physics objects using the
        current value of ``quantity``.

        Notes
        -----
        The algorithm directly writes the value of 'quantity' into the
        phase, which is against one of the OpenPNM rules of objects not
        being able to write into each other.

        """
        if iterative_props is None:
            iterative_props = self.iterative_props
        if not iterative_props:
            return
        # Fetch objects associated with the algorithm
        phase = self.project[self.settings.phase]
        # Update 'quantity' on phase with the most recent value
        quantity = self.settings['quantity']
        phase[quantity] = self.x
        # Regenerate all associated objects
        phase.regenerate_models(propnames=iterative_props)

    def set_BC(self, pores, bctype, bcvalues=[], mode='overwrite'):
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
            can be anything, but normal options are 'rate' and 'value'. If a
            list of strings is provided, then each mode in the list is
            handled in order, so if ``mode='clear' and ``bctype=['value',
            'rate']`` then these two bc arrays will be set to all ``nans`.
        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.
            If a single value is given, it's assumed to apply to all
            locations unless the 'total_rate' bc.type is supplied whereby
            a single value corresponds to a total rate to be divded evenly
            among all pores. Otherwise, different values can be applied to
            all pores in the form of an array of the same length as
            ``pores``.
        mode : str or list of str, optional
            Controls how the boundary conditions are applied. The default
            value is 'merge'. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'overwrite'  (default) Adds supplied boundary conditions to the
                         existing conditions, including overwriting any
                         existing conditions that may be present. This is
                         equivalent to calling ``'remove'`` on the given
                         locations followed by ``'add'``.
            'add'        Adds the supplied boundary conditions to the
                         existing conditions but does *not* overwrite any
                         conditions that are already present.
            'remove'     Removes boundary conditions from the specified
                         locations.
            'clear'      Removes all boundary conditions from the object of
                         the of the specified type from all locations.
            ============ =====================================================

            If a list of strings is provided, then each mode in the list is
            handled in order, so that ``['remove', 'add']`` will give the same
            results add ``'overwrite'``.  Another option would be ``['clear',
            'add']``, which would remove all existing bcs and add the supplied
            ones.

        Notes
        -----
        It is not possible to have multiple boundary conditions for a
        specified location in one algorithm.

        """
        # If a list of modes was given, handle them each in order
        if not isinstance(mode, str):
            for item in mode:
                self.set_BC(pores=pores, bctype=bctype,
                            bcvalues=bcvalues, mode=item)
            return
        # If a list of bctypes was given, handle them each in order
        if not isinstance(bctype, str):
            for item in bctype:
                self.set_BC(pores=pores, bctype=item,
                            bcvalues=bcvalues, mode=mode)
            return

        # Begin method
        bc_types = list(self['pore.bc'].keys())
        other_types = np.setdiff1d(bc_types, bctype).tolist()

        mode = self._parse_mode(
            mode,
            allowed=['overwrite', 'add', 'remove', 'clear'],
            single=True)

        pores = self._parse_indices(pores)

        no_bc = np.nan if self[f'pore.bc.{bctype}'].dtype in (float, int) else False

        values = np.array(bcvalues)
        if values.size == 1:  # Expand to array if scalar given
            values = np.ones_like(pores, dtype=values.dtype)*values

        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of values must match the number of locations')

        mask = np.ones_like(pores, dtype=bool)  # Indices of pores to keep
        if mode == 'add':
            # Remove indices that are already present for given bc type
            mask[isfinite(self[f'pore.bc.{bctype}'][pores])] = False
            # Now remove indices that are present for other BCs
            for item in other_types:
                mask[isfinite(self[f'pore.bc.{item}'][pores])] = False
            if mask.sum() > 0:
                self[f"pore.bc.{bctype}"][pores[mask]] = values[mask]
            else:
                logger.warning('No valid pore locations were specified')
        elif mode == 'overwrite':
            # Now remove indices that are present for other BCs
            for item in other_types:
                mask[isfinite(self[f'pore.bc.{item}'][pores])] = False
            if mask.sum() > 0:
                self[f"pore.bc.{bctype}"][pores[mask]] = values[mask]
            else:
                logger.warning('No valid pore locations were specified')
        elif mode == 'remove':
            # Now remove indices that are present for other BCs
            for item in other_types:
                mask[isfinite(self[f'pore.bc.{item}'][pores])] = False
            if mask.sum():
                self[f"pore.bc.{bctype}"][pores[mask]] = no_bc
            else:
                logger.warning('No valid pore locations were specified')
        elif mode == 'clear':
            self[f"pore.bc.{bctype}"] = no_bc


def isfinite(arr, inf=None):
    if arr.dtype in (bool, ):
        results = arr == True
    elif arr.dtype in (float, ):
        results = ~np.isnan(arr)
    else:
        results = arr != inf
    return results
