import openpnm.models as mods
from openpnm.geometry import GenericGeometry
import openpnm.models as mods


class Continuum(GenericGeometry):
    r"""


    """
    def __init__(self, porosity, **kwargs):
        # Ensure throats were not specified, since these are found automatically
        if 'throats' in kwargs.keys():
            raise Exception('Throats are automatically determined')
        # Begin process of finding neighboring throats
        if 'network' in kwargs.keys():
            network = kwargs.pop('network')
        else:
            project = kwargs.pop('project')
            network = project.network
        Ps = kwargs.pop('pores')
        Ps = self._parse_indices(Ps)
        Ts = network.find_neighbor_throats(pores=Ps, mode='xnor')
        # Finally initialize the class
        super().__init__(network=network, pores=Ps, throats=Ts, **kwargs)
        # Now add pore-scale models
        self['pore.porosity'] = porosity
        self['throat.length'] = 1e-12
        f = mods.geometry.pore_size.largest_sphere
        self.add_model(propname='throat.total_length',
                       model=mods.geometry.throat_length.ctc)
        # Now need to find the dimensions of each pore in each direction
