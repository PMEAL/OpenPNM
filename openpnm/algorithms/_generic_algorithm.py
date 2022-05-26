import logging
import numpy as np
from openpnm.core import ParserMixin, LabelMixin, Base2
from openpnm.utils import Docorator, SettingsAttr
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
    network : OpenPNM Network
        The network object to which this algorithm will apply
    name : str, optional
        Name of the algorithm. If not provided one is generated.

    """

    def __init__(self, network, settings=None, **kwargs):
        self._settings = SettingsAttr(GenericAlgorithmSettings)
        self.settings._update(settings)
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'alg_01'
        super().__init__(network=network, settings=self.settings, **kwargs)
        self['pore.all'] = np.ones(network.Np, dtype=bool)
        self['throat.all'] = np.ones(network.Nt, dtype=bool)

    @functools.cached_property
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

    def _update_iterative_props(self):
        """
        Regenerates phase, geometries, and physics objects using the
        current value of ``quantity``.

        Notes
        -----
        The algorithm directly writes the value of 'quantity' into the
        phase, which is against one of the OpenPNM rules of objects not
        being able to write into each other.

        """
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
