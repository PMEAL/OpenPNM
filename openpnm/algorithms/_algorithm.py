import logging
import numpy as np
from openpnm.core import ParserMixin, LabelMixin, Base2
from openpnm.utils import Docorator


__all__ = ['Algorithm']


logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='AlgorithmSettings', sections=docstr.all_sections)
@docstr.dedent
class AlgorithmSettings:
    r"""

    Parameters
    ----------
    %(BaseSettings.parameters)s

    """


@docstr.get_sections(base='Algorithm', sections=['Parameters'])
@docstr.dedent
class Algorithm(ParserMixin, LabelMixin, Base2):
    r"""
    Generic class to define the foundation of Algorithms

    Parameters
    ----------
    %(Base.parameters)s

    """

    def __init__(self, network, name='alg_?', **kwargs):
        super().__init__(network=network, name=name, **kwargs)
        self.settings._update(AlgorithmSettings())
        self['pore.all'] = np.ones([network.Np, ], dtype=bool)
        self['throat.all'] = np.ones([network.Nt, ], dtype=bool)

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

    def clear_BCs(self, bctype=[]):
        r"""
        Clear all BCs of the given type(s)

        Parameters
        ----------
        bctype : str or list(str)
            The specific type of boundary condition to clear. If not provided
            all types will be cleared.
        """
        self.set_BC(pores=None, bctype=bctype, mode='remove')

    def set_BC(self, pores=None, bctype=[], bcvalues=[], mode='add'):
        r"""
        The main method for setting and adjusting boundary conditions.

        This method is called by other more convenient wrapper functions like
        ``set_value_BC``.

        Parameters
        ----------
        pores : array_like
            The pores where the boundary conditions should be applied. If
            ``None`` is given then *all* pores are assumed.  This is useful
            when ``mode='remove'``.
        bctype : str
            Specifies the type or the name of boundary condition to apply. This
            can be anything, but normal options are 'rate' and 'value'. If
            an empty list is provided, then all bc types will be assumed. This
            is useful for clearing all bcs if ``mode='remove'`` and ``pores=
            None``.
        bcvalues : int or array_like
            The boundary value to apply, such as concentration or rate.
            If a single value is given, it's assumed to apply to all
            locations. Different values can be applied to all pores in the
            form of an array of the same length as ``pores``. Note that using
            ``mode='add'`` and ``values=np.nan`` is equivalent to removing
            bcs from the given ``pores``.
        mode : str or list of str, optional
            Controls how the boundary conditions are applied. Options are:

            ============ =====================================================
            mode         meaning
            ============ =====================================================
            'add'        (default) Adds the supplied boundary conditions to
                         the given locations. Raises an exception if values
                         of any type already exist in the given locations.
            'overwrite'  Adds supplied boundary conditions to the given
                         locations, including overwriting conditions of the
                         given type or any other type that may be present in
                         the given locations.
            'remove'     Removes boundary conditions of the specified type
                         from the specified locations. If ``bctype`` is not
                         specified then *all* types are removed. If no
                         locations are given then values are remvoed from
                         *all* locations.
            ============ =====================================================

            If a list of strings is provided, then each mode in the list is
            handled in order, so that ``['remove', 'add']`` will give the same
            results add ``'overwrite'``.

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
        if len(bctype) == 0:
            bctype = self['pore.bc'].keys()
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
            allowed=['overwrite', 'add', 'remove'],
            single=True)

        # Infer the value that indicates "no bc" based on array dtype
        no_bc = get_no_bc(self[f'pore.bc.{bctype}'])

        if pores is None:
            pores = self.Ps
        pores = self._parse_indices(pores)

        # Deal with size of the given bcvalues
        values = np.array(bcvalues)
        if values.size == 1:
            values = np.ones_like(pores, dtype=values.dtype)*values
        # Ensure values and pores are the same size
        if values.size > 1 and values.size != pores.size:
            raise Exception('The number of values must match the number of locations')

        # Finally adjust the BCs according to mode
        if mode == 'add':
            mask = np.ones_like(pores, dtype=bool)  # Indices of pores to keep
            for item in self['pore.bc'].keys():  # Remove pores that are taken
                mask[isfinite(self[f'pore.bc.{item}'][pores])] = False
            if not np.all(mask):  # Raise exception if some conflicts found
                msg = "Some of the given locations already have BCs, " \
                    + "either use mode='remove' first or " \
                    + "use mode='overwrite' instead"
                raise Exception(msg)
            self[f"pore.bc.{bctype}"][pores[mask]] = values[mask]
        elif mode == 'overwrite':
            # Put given values in specified BC, sort out conflicts below
            self[f"pore.bc.{bctype}"][pores] = values
            # Collect indices that are present for other BCs for removal
            mask = np.ones_like(pores, dtype=bool)
            for item in other_types:
                self[f"pore.bc.{item}"][pores] = get_no_bc(self[f"pore.bc.{item}"])
                # Make a note of any BCs values of other types
                mask[isfinite(self[f'pore.bc.{item}'][pores])] = False
            if not np.all(mask):  # Warn that other values were overwritten
                msg = 'Some of the given locations already have BCs of ' \
                    + 'another type, these will be overwritten'
                logger.warning(msg)
        elif mode == 'remove':
            self[f"pore.bc.{bctype}"][pores] = no_bc


def get_no_bc(arr):
    no_bc = np.nan if arr.dtype in (float, int) else False
    return no_bc


def isfinite(arr, inf=None):
    if arr.dtype in (bool, ):
        results = arr == True
    elif arr.dtype in (float, ):
        results = ~np.isnan(arr)
    else:
        results = arr != inf
    return results
