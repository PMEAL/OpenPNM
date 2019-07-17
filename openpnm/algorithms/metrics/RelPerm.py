from openpnm.algorithms import GenericAlgorithm, StokesFlow
from openpnm.utils import logging
from openpnm import models
import numpy as np
import matplotlib.pyplot as plt
import openpnm
logger = logging.getLogger(__name__)


default_settings = {'sat': dict(),
                    'relperm_wp': dict(),
                    'relperm_nwp': dict(),
                    'perm_wp': dict(),
                    'perm_nwp': dict(),
                    'wp': [],
                    'nwp': [],
                    'pore.invasion_sequence': [],
                    'throat.invasion_sequence': [],
                    'flow_inlet': dict(),
                    'flow_outlet': dict(),
                    'pore_volume': [],
                    'throat_volume': [],
                    'BP_1': dict(),
                    'BP_2': dict(),
                    'results': {'sat': [], 'krw': [], 'krnw': []}
                    }

class RelPerm(GenericAlgorithm):
    r"""
    A subclass of Generic Algorithm to calculate relative permeabilities of
    fluids in a drainage process. The main roles of this subclass are to
    get invasion sequence and implement a method for calculating the relative
    permeabilities of the fluids flowing through a user-specified direction.
    """
    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        self.settings.update(settings)

    def setup(self, invading_phase=None, defending_phase=None,
              invasion_sequence=None, multiphase=None):
        self.settings['nwp']=invading_phase
        self.settings['wp']= defending_phase
        self.settings['throat_volume']='throat.volume'
        self.settings['pore_volume']='pore.volume'
        if (invasion_sequence=='invasion_sequence'):
            self.settings['pore.invasion_sequence']=self.settings['nwp']['pore.invasion_sequence']
            self.settings['throat.invasion_sequence']=self.settings['nwp']['throat.invasion_sequence']
        self.settings['BP_1']={'x': 'left', 'y': 'front', 'z': 'top'}
        self.settings['BP_2']={'x': 'right', 'y': 'back', 'z': 'bottom'}
        self.settings['flow_inlets']=self.settings['BP_1']
        self.settings['flow_outlets']=self.settings['BP_2']
        #provide an else statement later
        #provide statement for multiphase later
    def set_inlets(self, pores):
        self.settings['flow_inlet'] = pores

    def set_outlets(self, pores):
        self.settings['flow_outlet'] = pores

    def _regenerate_models(self):
        if self.settings['wp'] is not None:
            modelwp=models.physics.multiphase.conduit_conductance
            self.settings['wp'].add_model(model=modelwp,
                                          propname='throat.conduit_hydraulic_conductance',
                                          throat_conductance='throat.hydraulic_conductance')
        modelnwp=models.physics.multiphase.conduit_conductance
        self.settings['nwp'].add_model(model=modelnwp,
                                       propname='throat.conduit_hydraulic_conductance',
                                       throat_conductance='throat.hydraulic_conductance')

    def abs_perm_calc(self, B_pores, in_outlet_pores):
        if self.settings['wp'] is not None:
            network=self.project.network
            St_wp = StokesFlow(network=network, phase=self.settings['wp'])
            St_wp.set_value_BC(pores=B_pores[0], values=1)
            St_wp.set_value_BC(pores=B_pores[1], values=0)
            St_wp.run()
            val=St_wp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                                  outlets=in_outlet_pores[1])
            Kwp=val
            self.project.purge_object(obj=St_wp)
        else:
            Kwp=None
        St_nwp = StokesFlow(network=network, phase=self.settings['nwp'])
        St_nwp.set_value_BC(pores=B_pores[0], values=1)
        St_nwp.set_value_BC(pores=B_pores[1], values=0)
        St_nwp.run()
        val=St_nwp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                               outlets=in_outlet_pores[1])
        Knwp=val
        self.project.purge_object(obj=St_nwp)
        return [Kwp, Knwp]

    def rel_perm_calc(self, B_pores, in_outlet_pores):
        network=self.project.network
        self._regenerate_models()
        if self.settings['wp'] is not None:
            St_mp_wp = StokesFlow(network=network, phase=self.settings['wp'])
            St_mp_wp.setup(conductance='throat.conduit_hydraulic_conductance')
            St_mp_wp.set_value_BC(pores=B_pores[0], values=1)
            St_mp_wp.set_value_BC(pores=B_pores[1], values=0)
            St_mp_wp.run()
            Kewp=St_mp_wp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                                        outlets=in_outlet_pores[1])
            self.project.purge_object(obj=St_mp_wp)
        else:
            Kewp=None
        St_mp_nwp = StokesFlow(network=network, phase=self.settings['nwp'])
        St_mp_nwp.set_value_BC(pores=B_pores[0], values=1)
        St_mp_nwp.set_value_BC(pores=B_pores[1], values=0)
        St_mp_nwp.setup(conductance='throat.conduit_hydraulic_conductance')
        St_mp_nwp.run()
        Kenwp=St_mp_nwp.calc_effective_permeability(inlets=in_outlet_pores[0],
                                                    outlets=in_outlet_pores[1])
        self.project.purge_object(obj=St_mp_nwp)
        return [Kewp, Kenwp]

    def _sat_occ_update(self, i):
        network=self.project.network
        pore_mask=self.settings['pore.invasion_sequence']<i
        throat_mask=self.settings['throat.invasion_sequence']<i
        sat_p=np.sum(network['pore.volume'][pore_mask])
        sat_t=np.sum(network['throat.volume'][throat_mask])
        sat1=sat_p+sat_t
        bulk=(np.sum(network['pore.volume']) + np.sum(network['throat.volume']))
        sat=sat1/bulk
        self.settings['nwp']['pore.occupancy'] = pore_mask
        self.settings['nwp']['throat.occupancy'] = throat_mask
        if self.settings['wp'] is not None:
            self.settings['wp']['throat.occupancy'] = 1-throat_mask
            self.settings['wp']['pore.occupancy'] = 1-pore_mask
        return sat

    def run(self, Snw_num=None, IP_pores=None):
        net= self.project.network
        Foutlets_init=dict()
        for dim in self.settings['flow_outlets']:
            Foutlets_init.update({dim: net.pores(self.settings['flow_outlets'][dim])})
        Foutlets=dict()
        outl=[]
        for key in Foutlets_init.keys():
            outl=[Foutlets_init[key][x] for x in range(0, len(Foutlets_init[key]), 2)]
            Foutlets.update({key: outl})
        Finlets_init=dict()
        for dim in self.settings['flow_inlets']:
            Finlets_init.update({dim: net.pores(self.settings['flow_inlets'][dim])})
        Finlets=dict()
        inl=[]
        for key in Finlets_init.keys():
            inl=([Finlets_init[key][x] for x in range(0, len(Finlets_init[key]), 2)])
            Finlets.update({key: inl})
        K_dir=set(self.settings['flow_inlets'].keys())
        for dim in K_dir:
            B_pores=[net.pores(self.settings['BP_1'][dim]),
                     net.pores(self.settings['BP_2'][dim])]
            in_outlet_pores=[Finlets_init[dim], Foutlets_init[dim]]
            [Kw, Knw]=self.abs_perm_calc(B_pores, in_outlet_pores)
            if self.settings['wp'] is not None:
                self.settings['perm_wp'].update({dim: Kw})
            self.settings['perm_nwp'].update({dim: Knw})
        for dirs in self.settings['flow_inlets']:
            if self.settings['wp'] is not None:
                relperm_wp=[]
            else:
                relperm_wp=None
            relperm_nwp=[]
            if Snw_num is None:
                Snw_num=10
            max_seq = np.max([np.max(self.settings['pore.invasion_sequence']),
                              np.max(self.settings['throat.invasion_sequence'])])
            start=max_seq//Snw_num
            stop=max_seq
            step=max_seq//Snw_num
            Snwparr = []
            B_pores=[net.pores(self.settings['BP_1'][dirs]),
                     net.pores(self.settings['BP_2'][dirs])]
            in_outlet_pores=[Finlets_init[dirs], Foutlets_init[dirs]]
            for j in range(start, stop, step):
                sat=self._sat_occ_update(j)
                Snwparr.append(sat)
                [Kewp, Kenwp]=self.rel_perm_calc(B_pores, in_outlet_pores)
                if self.settings['wp'] is not None:
                    relperm_wp.append(Kewp/self.settings['perm_wp'][dirs])
                relperm_nwp.append(Kenwp/self.settings['perm_nwp'][dirs])
            if self.settings['wp'] is not None:
                self.settings['relperm_wp'].update({dirs: relperm_wp})
            self.settings['relperm_nwp'].update({dirs: relperm_nwp})
            self.settings['sat'].update({dirs:Snwparr})

    def plot_Kr_curve(self):
        f = plt.figure()
        sp = f.add_subplot(111)
        for inp in self.settings['flow_inlets']:
            if self.settings['wp'] is not None:
                sp.plot(self.settings['sat'][inp], self.settings['relperm_wp'][inp],
                        'o-', label='Krwp'+inp)
            sp.plot(self.settings['sat'][inp], self.settings['relperm_nwp'][inp],
                    '*-', label='Krnwp'+inp)
        sp.set_xlabel('Snw')
        sp.set_ylabel('Kr')
        sp.set_title('Relative Permability Curves')
        sp.legend()
        return f

    def get_Kr_data(self):
        self.settings['results']['sat']=self.settings['sat']
        if self.settings['wp'] is not None:
            self.settings['results']['krw']=self.settings['relperm_wp']
        else:
             self.settings['results']['krw']=None
        self.settings['results']['krnw']=self.settings['relperm_nwp']
        return self.settings['results']
