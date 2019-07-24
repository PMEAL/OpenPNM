from openpnm.algorithms import GenericAlgorithm, StokesFlow
from openpnm.utils import logging
from openpnm import models
import numpy as np
import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


default_settings = {'wp': None,
                    'nwp': None,
                    'conduit_hydraulic_conductance':
                        'throat.conduit_hydraulic_conductance',
                    'hydraulic_conductance':
                        'throat.hydraulic_conductance',
                    'pore.invasion_sequence': 'pore.invasion_sequence',
                    'throat.invasion_sequence': 'throat.invasion_sequence',
                    'flow_inlet': None,
                    'flow_outlet': None,
                    }


class RelativePermeability(GenericAlgorithm):
    r"""
    A subclass of Generic Algorithm to calculate relative permeabilities of
    fluids in a drainage process. The main roles of this subclass are to
    get invasion sequence and implement a method for calculating the relative
    permeabilities of the fluids flowing in three directions.

    Notes
    -----
    The results can be plotted using `plot_Kr_curves`, and numerical data
    can be obtained with `get_Kr_data`.
    Properties related to the invading phase have subscript 'nwp', while those 
    related to the defending phase (if there is any) are named by subscript 'wp'
    """
    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        self.settings.update(settings)
        self.Kr_values = {'sat': dict(),
                          'relperm_wp': dict(),
                          'relperm_nwp': dict(),
                          'perm_wp': dict(),
                          'perm_nwp': dict(),
                          'results': {'sat': [], 'krw': [], 'krnw': []}}

    def setup(self, invading_phase=None, defending_phase=None,
              invasion_sequence=None):
        r"""
        Assigns values to the algorithms ``settings``

        Parameters
        ----------
        invading_phase : string
            The invading or non-wetting phase
        defending_phase : string, optional
            If defending phase is specified, then it's permeability will
            also be calculated, otherwise only the invading phase is
            considered.
        invasion_sequence : string (default is 'invasion_sequence')
            The dictionary key on the invading phase object where the invasion
            sequence is stored.  The default from both the IP and OP algorithms
            is 'invasion_sequence', so this is the default here.
        """
        if invading_phase is not None:
            self.settings['nwp'] = invading_phase.name
        if defending_phase is not None:
            self.settings['wp'] = defending_phase.name
        if (invasion_sequence == 'invasion_sequence'):
            nwp = self.project[self.settings['nwp']]
            seq_p = nwp['pore.invasion_sequence']
            self.settings['pore.invasion_sequence'] = seq_p
            seq_t = nwp['throat.invasion_sequence']
            self.settings['throat.invasion_sequence'] = seq_t
        self.settings['flow_inlets'] = {'x': 'left',
                                        'y': 'front',
                                        'z': 'top'}
        self.settings['flow_outlets'] = {'x': 'right',
                                         'y': 'back',
                                         'z': 'bottom'}

    def _regenerate_models(self):
        r"""
        Updates the multiphase physics model for each saturation
        """
        prop = self.settings['conduit_hydraulic_conductance']
        prop_q = self.settings['hydraulic_conductance']
        try:
            wp = self.project[self.settings['wp']]
            modelwp = models.physics.multiphase.conduit_conductance
            wp.add_model(model=modelwp, propname=prop,
                         throat_conductance=prop_q)
        except:
            pass
        nwp = self.project[self.settings['nwp']]
        modelnwp = models.physics.multiphase.conduit_conductance
        nwp.add_model(model=modelnwp, propname=prop,
                      throat_conductance=prop_q)

    def _abs_perm_calc(self, flow_pores):
        r"""
        Calculates absolute permeability of the medium using StokesFlow algorithm.
        The direction of flow is defined by flow_pores. 
        
        Parameters
        ----------
        flow_pores: numpy array 
        Boundary pores that will have constant value boundary condition to in
        StokesFlow algorithm. First element is the inlet face (pores) for flow of 
        invading phase through porous media. Second element is the outlet face (pores).
        
        Output: array [Kwp, Knwp]
        The value of absolute permeability of defending (if there is any) and
        invadin phase in the direction that is defined by flow_pores.
        """
        network = self.project.network
        try:
            wp = self.project[self.settings['wp']]
            St_wp = StokesFlow(network=network, phase=wp)
            St_wp.set_value_BC(pores=flow_pores[0], values=1)
            St_wp.set_value_BC(pores=flow_pores[1], values=0)
            St_wp.run()
            val = St_wp.calc_effective_permeability(inlets=flow_pores[0],
                                                    outlets=flow_pores[1])
            Kwp = val
            self.project.purge_object(obj=St_wp)
        except:
            Kwp = None
            pass
        nwp = self.project[self.settings['nwp']]
        St_nwp = StokesFlow(network=network, phase=nwp)
        St_nwp.set_value_BC(pores=flow_pores[0], values=1)
        St_nwp.set_value_BC(pores=flow_pores[1], values=0)
        St_nwp.run()
        val = St_nwp.calc_effective_permeability(inlets=flow_pores[0],
                                               outlets=flow_pores[1])
        Knwp = val
        self.project.purge_object(obj=St_nwp)
        return [Kwp, Knwp]
    #shouldn't it be just invading phase? it's the same, right?############################################
    ######################################################################################################

    def _eff_perm_calc(self, flow_pores):
        r"""
        Calculates effective permeability of each phase using StokesFlow algorithm
        with updated multiphase physics models to account for the multiphase flow.
        The direction of flow is defined by flow_pores. 

        Parameters
        ----------
        flow_pores: numpy array 
        Boundary pores that will have constant value boundary condition to in
        StokesFlow algorithm. First element is the inlet face (pores) for flow of 
        invading phase through porous media. Second element is the outlet face (pores).
        
        Output: array [Kewp, Kenwp]
        The value of effective permeability of defending (if there is any) and
        invadin phase in the direction that is defined by flow_pores.
        """
        network=self.project.network
        self._regenerate_models()
        try:
            wp = self.project[self.settings['wp']]
            St_mp_wp = StokesFlow(network=network, phase=wp)
            St_mp_wp.setup(conductance='throat.conduit_hydraulic_conductance')
            St_mp_wp.set_value_BC(pores=flow_pores[0], values=1)
            St_mp_wp.set_value_BC(pores=flow_pores[1], values=0)
            St_mp_wp.run()
            Kewp=St_mp_wp.calc_effective_permeability(inlets=flow_pores[0],
                                                        outlets=flow_pores[1])
            self.project.purge_object(obj=St_mp_wp)
        except :
            Kewp = None
            pass
        nwp = self.project[self.settings['nwp']]
        St_mp_nwp = StokesFlow(network=network, phase=nwp)
        St_mp_nwp.set_value_BC(pores=flow_pores[0], values=1)
        St_mp_nwp.set_value_BC(pores=flow_pores[1], values=0)
        St_mp_nwp.setup(conductance='throat.conduit_hydraulic_conductance')
        St_mp_nwp.run()
        Kenwp=St_mp_nwp.calc_effective_permeability(inlets=flow_pores[0],
                                                    outlets=flow_pores[1])
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
        nwp = self.project[self.settings['nwp']]
        nwp['pore.occupancy'] = pore_mask
        nwp['throat.occupancy'] = throat_mask
        try:
            wp = self.project[self.settings['wp']]
            wp['throat.occupancy'] = 1-throat_mask
            wp['pore.occupancy'] = 1-pore_mask
        except :
            pass
        return sat

    def run(self, Snw_num=None, IP_pores=None):
        net= self.project.network
        # The following 1/2 of the inlet to ??? because, etc
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
            flow_pores=[net.pores(self.settings['flow_inlets'][dim]),
                     net.pores(self.settings['flow_outlets'][dim])]
            in_outlet_pores=[Finlets_init[dim], Foutlets_init[dim]]
            [Kw, Knw]=self._abs_perm_calc(flow_pores)
            try:
                wp = self.project[self.settings['wp']]
                self.Kr_values['perm_wp'].update({dim: Kw})
            except:
                pass
            self.Kr_values['perm_nwp'].update({dim: Knw})
        for dirs in self.settings['flow_inlets']:
            try:
                wp = self.project[self.settings['wp']]
                relperm_wp=[]
            except:
                relperm_wp=None
                pass
            relperm_nwp=[]
            if Snw_num is None:
                Snw_num=10
            max_seq = np.max([np.max(self.settings['pore.invasion_sequence']),
                              np.max(self.settings['throat.invasion_sequence'])])
            start=max_seq//Snw_num
            stop=max_seq
            step=max_seq//Snw_num
            Snwparr = []
            flow_pores=[net.pores(self.settings['flow_inlets'][dirs]),
                     net.pores(self.settings['flow_outlets'][dirs])]
            in_outlet_pores=[Finlets_init[dirs], Foutlets_init[dirs]]
            for j in range(start, stop, step):
                sat=self._sat_occ_update(j)
                Snwparr.append(sat)
                [Kewp, Kenwp]=self._eff_perm_calc(flow_pores)
                if self.settings['wp'] is not None:
                    relperm_wp.append(Kewp/self.Kr_values['perm_wp'][dirs])
                relperm_nwp.append(Kenwp/self.Kr_values['perm_nwp'][dirs])
            if self.settings['wp'] is not None:
                self.Kr_value['relperm_wp'].update({dirs: relperm_wp})
            self.Kr_values['relperm_nwp'].update({dirs: relperm_nwp})
            self.Kr_values['sat'].update({dirs: Snwparr})

    def plot_Kr_curves(self):
        f = plt.figure()
        sp = f.add_subplot(111)
        for inp in self.settings['flow_inlets']:
            if self.settings['wp'] is not None:
                sp.plot(self.Kr_values['sat'][inp], self.Kr_values['relperm_wp'][inp],
                        'o-', label='Krwp'+inp)
            sp.plot(self.Kr_values['sat'][inp], self.Kr_values['relperm_nwp'][inp],
                    '*-', label='Krnwp'+inp)
        sp.set_xlabel('Snw')
        sp.set_ylabel('Kr')
        sp.set_title('Relative Permability Curves')
        sp.legend()
        return f

    def get_Kr_data(self):
        self.Kr_values['results']['sat']=self.Kr_values['sat']
        if self.settings['wp'] is not None:
            self.Kr_values['results']['krw']=self.Kr_values['relperm_wp']
        else:
            self.Kr_values['results']['krw']=None
        self.Kr_values['results']['krnw']=self.Kr_values['relperm_nwp']
        return self.settings['results']