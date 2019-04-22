from openpnm.algorithms import GenericAlgorithm, InvasionPercolation
from openpnm.algorithms import DirectionalRelativePermeability
from openpnm.utils import logging
# from openpnm import models
import numpy as np
import openpnm
# import matplotlib.pyplot as plt
# from collections import namedtuple
# import csv
# import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


default_settings = {
                    'inv_inlets': dict(),
                    'flow_inlets': dict(),
                    'flow_outlets': dict(),
                    'mode': 'strict',
                    'sat': dict(),
                    'relperm_wp': dict(),
                    'relperm_nwp': dict(),
                    'perm_wp': dict(),
                    'perm_nwp': dict(),
                    'wp': [],
                    'nwp': [],
                    'BP_1': dict(),
                    'BP_2': dict()}


class RelativePermeability(GenericAlgorithm, DirectionalRelativePermeability):
    r"""
    A subclass of Generic Algorithm to calculate relative permeabilities of
    fluids in a drainage process. The main roles of this subclass are to
    implement Invasion Percolation if no invasion sequence is given as an
    argument and to implement a method for calculating the relative
    permeabilities of the fluids.
    """
    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        self.settings.update(settings)

    def setup(self, input_vect=['xx', 'xy', 'xz', 'yx', 'yy', 'yz',
                                'zx', 'zy', 'zz']):
        input_vect.sort()
        self.settings['input_vect']=input_vect
        inlet_dict={'x': 'left', 'y': 'front', 'z': 'top'}
        outlet_dict={'x': 'right', 'y': 'back', 'z': 'bottom'}
        for i in range(len(input_vect)):
            inv=input_vect[i][0]
            flow=input_vect[i][1]
            self.settings['BP_1'].update({flow: (inlet_dict[flow])})
            self.settings['BP_2'].update({flow: (outlet_dict[flow])})
            self.settings['inv_inlets'].update({inv: (inlet_dict[inv])})
            self.settings['flow_inlets'].update({flow: (inlet_dict[flow])})
            self.settings['flow_outlets'].update({flow: (outlet_dict[flow])})


    def run(self,Snw_num=None):
        net= self.project.network
        oil = openpnm.phases.GenericPhase(network=net, name='oil')
        water = openpnm.phases.GenericPhase(network=net, name='water')
        oil['pore.viscosity']=0.547
        oil['throat.surface_tension'] = 0.072
        oil['pore.surface_tension']=0.072
        oil['pore.contact_angle']=110
        water['throat.contact_angle'] = 70
        water['pore.contact_angle'] = 70
        water['throat.surface_tension'] = 0.0483
        water['pore.surface_tension'] = 0.0483
        water['pore.viscosity']=0.4554
        mod = openpnm.models.physics.hydraulic_conductance.hagen_poiseuille
        oil.add_model(propname='throat.hydraulic_conductance',
                          model=mod)
        oil.add_model(propname='throat.entry_pressure',
                          model=openpnm.models.physics.capillary_pressure.washburn)
        water.add_model(propname='throat.hydraulic_conductance',
                          model=mod)
        water.add_model(propname='throat.entry_pressure',
                          model=openpnm.models.physics.capillary_pressure.washburn)
        self.settings['nwp']= oil
        self.settings['wp']= water
        # define inlet/outlets
        Iinlets_init=dict()
        for dim in self.settings['inv_inlets']:
            Iinlets_init.update({dim: net.pores(self.settings['inv_inlets'][dim])}) 
        Iinlets=dict()
        for key in Iinlets_init.keys():
            Iinlets.update({key: ([Iinlets_init[key][x] for x in \
                                   range(0, len(Iinlets_init[key]), 2)])})
        Foutlets_init=dict()
        for dim in self.settings['flow_outlets']:
            Foutlets_init.update({dim: net.pores(self.settings['flow_outlets'][dim])})
        Foutlets=dict()
        for key in Foutlets_init.keys():
            Foutlets.update({key: ([Foutlets_init[key][x] for x in \
                                    range(0, len(Foutlets_init[key]), 2)])})
        Finlets_init=dict()
        for dim in self.settings['flow_inlets']:
            Finlets_init.update({dim: net.pores(self.settings['flow_inlets'][dim])})
        Finlets=dict()
        for key in Finlets_init.keys():
            Finlets.update({key: ([Finlets_init[key][x] for x in \
                                   range(0, len(Finlets_init[key]), 2)])})
        inv=InvasionPercolation(network=net)
        inv.setup(phase=self.settings['nwp'],
                  entry_pressure='throat.entry_pressure',
                  pore_volume='pore.volume',
                  throat_volume='throat.volume')
        K_dir=set(self.settings['flow_inlets'].keys())
        for dim in K_dir:
            B_pores=[net.pores(self.settings['BP_1'][dim]),net.pores(self.settings['BP_2'][dim])]
            in_outlet_pores=[Finlets_init[dim],Foutlets_init[dim]]
            [Kw,Knw]=super(DirectionalRelativePermeability,self).abs_perm_calc(B_pores,in_outlet_pores)
            self.settings['perm_wp'].update({dim:Kw})
            self.settings['perm_nwp'].update({dim: Knw})
        for i in range(len(self.settings['input_vect'])):
            invasion=self.settings['input_vect'][i][0]
            flow=self.settings['input_vect'][i][1]
            relperm_wp=[]
            relperm_nwp=[]
            inv.set_inlets(pores=Iinlets[invasion])
            inv.run()
            if Snw_num is None:
                Snw_num=10
            max_seq = np.max([np.max(self.settings['pore.invasion_sequence']),
                          np.max(self.settings['throat.invasion_sequence'])])
            start=max_seq//Snw_num
            stop=max_seq
            step=max_seq//Snw_num
            Snwparr = []
            B_pores=[net.pores(self.settings['BP_1'][flow]),net.pores(self.settings['BP_2'][flow])]
            in_outlet_pores=[Finlets_init[flow],Foutlets_init[flow]]
            for i in range(start, stop, step):
                sat=self._sat_occ_update(i)
                Snwparr.append(sat)
                [Kewp,Kenwp]=super(DirectionalRelativePermeability,self).rel_perm_calc(B_pores,
                                                                                        in_outlet_pores)
                relperm_wp.append(Kewp/self.settings['perm_wp'][flow])
                relperm_nwp.append(Kenwp/self.settings['perm_nwp'][flow])
            self.settings['relperm_wp'].update({self.settings['input_vect'][i]:\
                                                   relperm_wp})
            self.settings['relperm_nwp'].update({self.settings['input_vect'][i]:\
                                                    relperm_nwp})
            self.settings['sat'].update({self.settings['input_vect'][i]: Snwparr})