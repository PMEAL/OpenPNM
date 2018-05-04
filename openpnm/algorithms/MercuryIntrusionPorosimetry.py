import openpnm as op
from openpnm.algorithms import GenericAlgorithm, Drainage
from openpnm.core import logging
logger = logging.getLogger()


class MercuryIntrusionPorosimetry(GenericAlgorithm):

    def __init__(self, network, phase, **kwargs):
        super().__init__(network=network, **kwargs)

        Hg = op.phases.Mercury(network=network)
        Hg['throat.surface_tension'] = 0.480
        Hg['throat.contact_angle'] = 140
        Hg.add_model(propname='throat.capillary_pressure',
                     model=op.physics.models.capillary_pressure.washburn,
                     surface_tension='throat.surface_tension',
                     contact_angle='throat.contact_angle',
                     regen_mode='normal')
        mip = Drainage(network=network)
        mip.setup(invading_phase=Hg)
        if 'pore.surface' not in network.keys():
            op.topotools.find_surface_pores(network=network)
        mip.set_inlets(network.pores('surface'))
        mip.run()
        self.plot_data = mip.plot_drainage_curve
        self.get_data = mip.get_drainage_data
        self.update(mip)
