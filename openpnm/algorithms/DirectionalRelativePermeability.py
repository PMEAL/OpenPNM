from openpnm.algorithms import GenericAlgorithm, StokesFlow
from openpnm.utils import logging
from openpnm import models
import numpy as np
import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


default_settings = {'sat': [],
                    'relperm_wp': [],
                    'relperm_nwp': [],
                    'perm_wp': [],
                    'perm_nwp': [],
                    'wp': [],
                    'nwp': [],
                    'pore.invasion_sequence': [],
                    'throat.invasion_sequence': [],
                    'flow_inlet': [],
                    'flow_outlet': [],
                    'pore_volume': [],
                    'throat_volume': [],
                    'results':{'sat':[], 'krw':[], 'krnw':[]}
                    }


class DirectionalRelativePermeability(GenericAlgorithm):
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

    def setup(self, invading_phase=None, defending_phase=None,
              pore_invasion_sequence=None,
              throat_invasion_sequence=None):
        self.settings['nwp']=invading_phase
        self.settings['wp']= defending_phase
        self.settings['pore.invasion_sequence']=pore_invasion_sequence
        self.settings['throat.invasion_sequence']=throat_invasion_sequence
        self.settings['throat_volume']='pore.volume'
        self.settings['pore_volume']='throat.volume'

    def set_inlets(self, pores):
        self.settings['flow_inlet'] = pores

    def set_outlets(self, pores):
        self.settings['flow_outlet'] = pores

    def _regenerate_models(self):
        self.settings['wp'].add_model(model=models.physics.multiphase.conduit_conductance,
                                      propname='throat.conduit_hydraulic_conductance',
                                      throat_conductance='throat.hydraulic_conductance')
        self.settings['nwp'].add_model(model=models.physics.multiphase.conduit_conductance,
                                       propname='throat.conduit_hydraulic_conductance',
                                       throat_conductance='throat.hydraulic_conductance')

    def abs_perm_calc(self,B_pores,in_outlet_pores):
        network=self.project.network
        St_wp = StokesFlow(network=network, phase=self.settings['wp'])
        St_wp.set_value_BC(pores=B_pores[0], values=1)
        St_wp.set_value_BC(pores=B_pores[1], values=0)
        St_wp.run()
        val=St_wp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                              outlets=in_outlet_pores[1])
        Kwp=val
        self.project.purge_object(obj=St_wp)
        St_nwp = StokesFlow(network=network, phase=self.settings['nwp'])
        St_nwp.set_value_BC(pores=B_pores[0], values=1)
        St_nwp.set_value_BC(pores=B_pores[1], values=0)
        St_nwp.run()
        val=St_nwp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                               outlets=in_outlet_pores[1])
        Knwp=val
        self.project.purge_object(obj=St_nwp)
        return [Kwp,Knwp]
    
    def rel_perm_calc(self,B_pores,in_outlet_pores):
        network=self.project.network
        self._regenerate_models()
        St_mp_wp = StokesFlow(network=network, phase=self.settings['wp'])
        St_mp_wp.setup(conductance='throat.conduit_hydraulic_conductance')
        St_mp_wp.set_value_BC(pores=B_pores[0], values=1)
        St_mp_wp.set_value_BC(pores=B_pores[1], values=0)
        St_mp_nwp = StokesFlow(network=network, phase=self.settings['nwp'])
        St_mp_nwp.set_value_BC(pores=B_pores[0], values=1)
        St_mp_nwp.set_value_BC(pores=B_pores[1], values=0)
        St_mp_nwp.setup(conductance='throat.conduit_hydraulic_conductance')
        St_mp_wp.run()
        St_mp_nwp.run()
        Kewp=St_mp_wp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                                        outlets=in_outlet_pores[1])
        Kenwp=St_mp_nwp.calc_effective_permeability(inlets=in_outlet_pores[0],                         
                                                        outlets=in_outlet_pores[1])
        self.project.purge_object(obj=St_mp_wp)
        self.project.purge_object(obj=St_mp_nwp)
        return [Kewp,Kenwp]
    
    def _sat_occ_update(self,i):
        network=self.project.network
        pore_mask=self.settings['pore.invasion_sequence']<i
        throat_mask=self.settings['throat.invasion_sequence']<i
        sat_p=np.sum(network['pore.volume'][pore_mask])
        sat_t=np.sum(network['throat.volume'][throat_mask])
        sat1=sat_p+sat_t
        bulk=(np.sum(network['pore.volume']) + np.sum(network['throat.volume']))
        sat=sat1/bulk
        self.settings['nwp']['pore.occupancy'] = pore_mask
        self.settings['wp']['pore.occupancy'] = 1-pore_mask
        self.settings['nwp']['throat.occupancy'] = throat_mask
        self.settings['wp']['throat.occupancy'] = 1-throat_mask
        return sat

    def run(self, Snw_num=None, IP_pores=None):
        [Kwp,Knwp]=self.abs_perm_calc(B_pores=[self.settings['flow_inlet'],self.settings['flow_outlet']],
                            in_outlet_pores=[self.settings['flow_inlet'],self.settings['flow_outlet']])
        self.settings['perm_wp']=Kwp
        self.settings['perm_nwp']=Knwp
        self.settings['IP_pores']=IP_pores
        if Snw_num is None:
            Snw_num=10
        max_seq = np.max([np.max(self.settings['pore.invasion_sequence']),
                          np.max(self.settings['throat.invasion_sequence'])])
        start=max_seq//Snw_num
        stop=max_seq
        step=max_seq//Snw_num
        self.settings['sat']=[]
        self.settings['relperm_wp']=[]
        self.settings['relperm_nwp']=[]
        for i in range(start, stop, step):
            sat=self._sat_occ_update(i)
            self.settings['sat'].append(sat)
            [Kewp,Kenwp]=self.rel_perm_calc(B_pores=[self.settings['flow_inlet'],self.settings['flow_outlet']],
                            in_outlet_pores=[self.settings['flow_inlet'],self.settings['flow_outlet']])
            Krwp=Kewp/self.settings['perm_wp']
            Krnwp=Kenwp/self.settings['perm_nwp']
            self.settings['relperm_wp'].append(Krwp)
            self.settings['relperm_nwp'].append(Krnwp)
        
    def get_Kr_data(self):
        self.settings['results']['sat']=self.settings['sat']
        self.settings['results']['krw']=self.settings['relperm_wp']
        self.settings['results']['krnw']=self.settings['relperm_nwp']
        return self.settings['results']
        
    def plot_Kr_curve(self):
        f = plt.figure()
        sp = f.add_subplot(111)
        sp.plot(self.settings['sat'], self.settings['relperm_wp'], 'o-')
        sp.plot(self.settings['sat'], self.settings['relperm_nwp'], '*-')
        return f
        
        