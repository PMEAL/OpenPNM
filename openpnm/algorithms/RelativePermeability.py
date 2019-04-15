from openpnm.algorithms import GenericAlgorithm, StokesFlow, InvasionPercolation
from openpnm.utils import logging
from openpnm import models
import numpy as np
import openpnm
# import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


default_settings = {
                    'inv_inlets':dict(),
                    'flow_inlets':dict(),
                    'flow_outlets':dict(),
                    'mode': 'strict',
                    'sat': dict(),
                    'relperm_wp': dict(),
                    'relperm_nwp':dict(),
                    'perm_wp':dict(),
                    'perm_nwp':dict(),
                    'wp': [],
                    'nwp': [],
                    'BP_1': dict(),
                    'BP_2': dict(),
                    }
class RelativePermeability(GenericAlgorithm):
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

    def setup(self, inv_inlets={'x': 'left', 'y': 'front', 'z': 'top'},
                    flow_inlets={'x': 'left', 'y': 'front', 'z': 'top'},
                    flow_outlets={'x': 'right', 'y': 'back', 'z': 'bottom'},
                    ):
        self.settings['BP_1']=flow_inlets
        self.settings['BP_2']=flow_outlets
        self.settings['inv_inlets'] = inv_inlets
        self.settings['flow_inlets'] = flow_inlets
        self.settings['flow_outlets'] = flow_outlets
        #self._setup_ip_algs()
        #self.run()
               
#    def run(self):
#        self._setup_ip_algs()
#        return {'ky': self.settings['relperm_wp']['y'], 'Saty': self.settings['sat']['y']}
    #def _setup_ip_algs(self):
    def setup_ip_algs(self):
        network = self.project.network
        print('printing now',self.settings['inv_inlets'])
        oil = openpnm.phases.GenericPhase(network=network, name='oil')
        water = openpnm.phases.GenericPhase(network=network, name='water')
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
        Iinlets_init = {'x':network.pores(self.settings['inv_inlets']['x']), 
                   'y':network.pores(self.settings['inv_inlets']['y']),
                   'z':network.pores(self.settings['inv_inlets']['z'])}
        Iinlets=dict()
        for key in Iinlets_init.keys():
            Iinlets.update({key:([Iinlets_init[key][x] for x in range(0, len(Iinlets_init[key]), 2)])})
            
        Foutlets_init = {'x':network.pores(self.settings['flow_outlets']['x']), 
                   'y':network.pores(self.settings['flow_outlets']['y']),
                   'z':network.pores(self.settings['flow_outlets']['z'])}
        Foutlets=dict()
        for key in Foutlets_init.keys():
            Foutlets.update({key:([Foutlets_init[key][x] for x in range(0, len(Foutlets_init[key]), 2)])})
            
        Finlets_init = {'x':network.pores(self.settings['flow_inlets']['x']), 
                   'y':network.pores(self.settings['flow_inlets']['y']),
                   'z': network.pores(self.settings['flow_inlets']['z'])}
        Finlets=dict()
        for key in Finlets_init.keys():
            Finlets.update({key:([Finlets_init[key][x] for x in range(0, len(Finlets_init[key]), 2)])})
        network=self.project.network
        inv=InvasionPercolation(network=network)
        inv.setup(phase=self.settings['nwp'],
                  entry_pressure='throat.entry_pressure',
                  pore_volume='pore.volume',
                  throat_volume='throat.volume')
        Sarr=np.linspace(0,1,num=20)

        for dim in self.settings['inv_inlets']:
            relperm_wp=[]
            relperm_nwp=[]
            Stokes_alg_wp = StokesFlow(network=network, phase=self.settings['wp'])
            Stokes_alg_wp.set_value_BC(pores=network.pores(self.settings['BP_1'][dim]),values=1)
            Stokes_alg_wp.set_value_BC(pores=network.pores(self.settings['BP_2'][dim]),values=0)
            Stokes_alg_wp.run()
            val=Stokes_alg_wp.calc_effective_permeability(inlets=Finlets_init[dim], outlets=Foutlets_init[dim])
            self.settings['perm_wp'].update({dim:val})
            self.project.purge_object(obj=Stokes_alg_wp)
            Stokes_alg_nwp = StokesFlow(network=network, phase=self.settings['nwp'])
            Stokes_alg_nwp.set_value_BC(pores=network.pores(self.settings['BP_1'][dim]),values=1)
            Stokes_alg_nwp.set_value_BC(pores=network.pores(self.settings['BP_2'][dim]),values=0)
            Stokes_alg_nwp.run()
            val=Stokes_alg_nwp.calc_effective_permeability(inlets=Finlets_init[dim], outlets=Foutlets_init[dim])
            self.settings['perm_nwp'].update({dim:val})
            self.project.purge_object(obj=Stokes_alg_nwp)
            #print(self.settings['perm_wp'])
            #print(self.settings['perm_wp'])
            inv.set_inlets(pores=Iinlets[dim])
            inv.run()
            Snwparr =  []
            for Snw in Sarr:
                res1=inv.results(Snwp=Snw)
                #print(res1['pore.occupancy'])
                occ_ts=res1['throat.occupancy']
                if np.any(occ_ts):
                    #print(Snw)
                    Snwparr.append(Snw)
                    self.settings['nwp']['pore.occupancy'] = res1['pore.occupancy']
                    self.settings['wp']['pore.occupancy'] = 1-res1['pore.occupancy']
                    self.settings['nwp']['throat.occupancy'] = res1['throat.occupancy']
                    self.settings['wp']['throat.occupancy'] = 1-res1['throat.occupancy']
                    #phase regenerate models
                    mode=self.settings['mode']
                    self.settings['wp'].add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)
                    self.settings['nwp'].add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)
                    Stokes_alg_mp_wp = StokesFlow(network=network,phase=self.settings['wp'])
                    Stokes_alg_mp_wp.setup(conductance='throat.conduit_hydraulic_conductance')
                    Stokes_alg_mp_wp.set_value_BC(pores=network.pores(self.settings['BP_1'][dim]),values=1)
                    Stokes_alg_mp_wp.set_value_BC(pores=network.pores(self.settings['BP_2'][dim]),values=0)
                    Stokes_alg_mp_nwp = StokesFlow(network=network,phase=self.settings['nwp'])
                    Stokes_alg_mp_nwp.set_value_BC(pores=network.pores(self.settings['BP_1'][dim]),values=1)
                    Stokes_alg_mp_nwp.set_value_BC(pores=network.pores(self.settings['BP_2'][dim]),values=0)
                    Stokes_alg_mp_nwp.setup(conductance='throat.conduit_hydraulic_conductance')
                    Stokes_alg_mp_wp.run()
                    Stokes_alg_mp_nwp.run()
                    Keff_mp_wp = Stokes_alg_mp_wp.calc_effective_permeability(inlets=Finlets_init[dim], outlets=Foutlets_init[dim])
                    Keff_mp_nwp = Stokes_alg_mp_nwp.calc_effective_permeability(inlets=Finlets_init[dim], outlets=Foutlets_init[dim])
                    relperm_wp.append(Keff_mp_wp/self.settings['perm_wp'][dim])
                    relperm_nwp.append(Keff_mp_nwp/self.settings['perm_nwp'][dim])
                    self.project.purge_object(obj=Stokes_alg_mp_wp)
                    self.project.purge_object(obj=Stokes_alg_mp_nwp)
                    #print(Keff_mp_wp)
                    #print(Keff_mp_nwp)
                self.settings['relperm_wp'].update({dim:relperm_wp})
                self.settings['relperm_nwp'].update({dim:relperm_nwp})
                self.settings['sat'].update({dim:Snwparr})
            print('relperm_wp',self.settings['relperm_wp'][dim])
            print('relperm_nwp',self.settings['relperm_nwp'][dim])
            #print('sat',self.settings['sat'][dim])
            
            
                
    

