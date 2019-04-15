from openpnm.algorithms import GenericAlgorithm, StokesFlow, InvasionPercolation
from openpnm.utils import logging
from openpnm import models
import numpy as np
import openpnm
import scipy as sp
# import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


default_settings = {
                    'mode': 'strict',
                    'sat': [],
                    'relperm_wp': [],
                    'relperm_nwp':[],
                    'perm_wp':[],
                    'perm_nwp':[],
                    'wp': [],
                    'nwp': [],
                    'pore.invasion_sequence': [],
                    'throat.invasion_sequence': [],
                    'flow_inlet': [],
                    'flow_outlet':[],
                    'pore_volume':[],
                    'throat_volume':[]
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

    def setup(self, invading_phase=None, defending_phase=None, pore_invasion_sequence=None,throat_invasion_sequence=None):
        self.settings['nwp']=invading_phase
        self.settings['wp']= defending_phase
        self.settings['pore.invasion_sequence']=pore_invasion_sequence
        self.settings['throat.invasion_sequence']=throat_invasion_sequence
        self.settings['throat_volume']='pore.volume'
        self.settings['pore_volume']='throat.volume'
        #self._setup_rp_algs()
        #self.run()
               
#    def run(self):
#        self._setup_ip_algs()
#        return {'ky': self.settings['relperm_wp']['y'], 'Saty': self.settings['sat']['y']}
    #def _setup_ip_algs(self):


#    def _setup_rp_algs(self):

#    def IP_occ(self, Snwp):
#        net = self.project.network
#        P12 = net['throat.conns']
#            # Fetch void volume for pores and throats
#        Vp = net[self.settings['pore_volume']]
#        Vt = net[self.settings['throat_volume']]
#            # Fetch the order of filling
#        Np = self.settings['pore.invasion_seq']
#        Nt = self.settings['throat.invasion_seq']
#        # Create Nt-long mask of which pores were filled when throat was filled
#        Pinv = (Np[P12].T == Nt).T
#        # If a pore and throat filled together, find combined volume
#        Vinv = sp.vstack(((Pinv*Vp[P12]).T, Vt)).T
#        Vinv = sp.sum(Vinv, axis=1)
#        # Convert to cumulative volume filled as each throat is invaded
#        x = sp.argsort(Nt)  # Find order throats were invaded
#        Vinv_cum = np.cumsum(Vinv[x])
#        # Normalized cumulative volume filled into saturation
#        S = Vinv_cum/(Vp.sum() + Vt.sum())
#        # Find throat invasion step where Snwp was reached
#        try:
#            N = sp.where(S < Snwp)[0][-1]
#        except:
#            N = -np.inf
#        data = {'pore.occupancy': Np <= N, 'throat.occupancy': Nt <= N}
#        return data
        
    def set_inlets(self, pores):
        #network=self.project.network
        self.settings['flow_inlet'] = pores

    def set_outlets(self, pores):
        self.settings['flow_outlet'] = pores
#        pores_out=[]
#        for outlet_num in range(len(pores)):
#            self['pore.outlets'] = False
#            self['pore.outlets'][pores[outlet_num]] = True
#            pores_out.append(self['pore.outlets'])
#        # print('outlets are', pores_out)
#        self['pore.outlets']= False
#        return pores_out           
        
#    def IP_occ(self, Snwp_num):
#        pore_occ=[]
#        throat_occ=[]
#        pn=self.project.network
#        seq_sort = np.argsort(self.settings['throat.invasion_sequence'])
#        sat = np.sum(pn['pore.volume'][self.settings['IP_pores']])
#        #pressure = 0
#        sats = [sat]
#        print(self.settings['throat.invasion_sequence'])
#        print(self.settings['pore.invasion_sequence'])
#        #pressures = [pressure]
#        pores = pn.pores()
#        throats = pn.throats()
#        max_seq = np.max([np.max(self.settings['pore.invasion_sequence']),
#                  np.max(self.settings['throat.invasion_sequence'])])
#        print('max is',max_seq)
#        start=max_seq//Snwp_num
#        stop=max_seq
#        step=max_seq//Snwp_num
#        for seq in range(start,stop,step):
#            pore_o=np.zeros((1, pn.Np), dtype=bool)
#            throat_o=np.zeros((1, pn.Nt), dtype=bool)
#            if seq in self.settings['throat.invasion_sequence']:
#                throat = throats[self.settings['throat.invasion_sequence'] == seq]
#                print('throat',throat)
#                sat += pn['throat.volume'][throat]
#                #add occupancy
#                map(lambda x: x>=4, throats)
#                throat_o[throats[self.settings['throat.invasion_sequence'] <= seq]]=True
#                throat_occ.append(throat_o)
#                #pressure = ip['throat.entry_pressure'][throat][0]
#            if seq in self.settings['pore.invasion_sequence']:
#                pore = pores[ip['pore.invasion_sequence'] == seq]
#                sat += pn['pore.volume'][pore]
#                #add occupancy
#                pore_o[pores[self.settings['pore.invasion_sequence'] <= seq]]=True
#                pore_occ.append(pore_o)
#            else:
#                pore_occ.append(pore_o)
#            sats.append(sat)
#            #pressures.append(pressure)
#        dom_vol = np.sum(pn['pore.volume']) + np.sum(pn['throat.volume'])
#        sats = np.asarray(sats)/dom_vol
#        #pressures = np.asarray(pressures)
#        data = {'pore.occupancy': throat_occ , 'throat.occupancy': pore_occ, 'sats': sats }
#        return data
        
    def run(self,Snw_num=None,IP_pores=None):
        network=self.project.network
        Stokes_alg_wp = StokesFlow(network=network, phase=self.settings['wp'])
        Stokes_alg_wp.set_value_BC(pores=(self.settings['flow_inlet']),values=1)
        Stokes_alg_wp.set_value_BC(pores=(self.settings['flow_outlet']),values=0)
        Stokes_alg_wp.run()
        val=Stokes_alg_wp.calc_effective_permeability(inlets=self.settings['flow_inlet'], outlets=self.settings['flow_outlet'])
        self.settings['perm_wp']=val
        self.project.purge_object(obj=Stokes_alg_wp)
        Stokes_alg_nwp = StokesFlow(network=network, phase=self.settings['nwp'])
        Stokes_alg_nwp.set_value_BC(pores=self.settings['flow_inlet'],values=1)
        Stokes_alg_nwp.set_value_BC(pores=self.settings['flow_outlet'],values=0)
        Stokes_alg_nwp.run()
        val=Stokes_alg_nwp.calc_effective_permeability(inlets=self.settings['flow_inlet'], outlets=self.settings['flow_outlet'])
        self.settings['perm_nwp']=val
        self.project.purge_object(obj=Stokes_alg_nwp)
        self.settings['IP_pores']=IP_pores
        #network=self.project.network
        if Snw_num==None:
            Snw_num=10
        #res1=self.IP_occ(Snw_num)
        #Snw=res1['sats']
        max_seq = np.max([np.max(self.settings['pore.invasion_sequence']),
                          np.max(self.settings['throat.invasion_sequence'])])
        start=max_seq//Snw_num
        stop=max_seq
        step=max_seq//Snw_num
        for i in range(start,stop,step):
            pore_mask=self.settings['pore.invasion_sequence']<i
            throat_mask=self.settings['throat.invasion_sequence']<i
            sat=(np.sum(network['pore.volume'][pore_mask])+np.sum(network['throat.volume'][throat_mask]))/(np.sum(network['pore.volume']) + np.sum(network['throat.volume']))
            self.settings['sat'].append(sat)
            self.settings['nwp']['pore.occupancy'] = pore_mask
            self.settings['wp']['pore.occupancy'] = 1-pore_mask
            self.settings['nwp']['throat.occupancy'] = throat_mask
            self.settings['wp']['throat.occupancy'] = 1-throat_mask
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
            Stokes_alg_mp_wp.set_value_BC(pores=self.settings['flow_inlet'],values=1)
            Stokes_alg_mp_wp.set_value_BC(pores=self.settings['flow_outlet'],values=0)
            Stokes_alg_mp_nwp = StokesFlow(network=network,phase=self.settings['nwp'])
            Stokes_alg_mp_nwp.set_value_BC(pores=self.settings['flow_inlet'],values=1)
            Stokes_alg_mp_nwp.set_value_BC(pores=self.settings['flow_outlet'],values=0)
            Stokes_alg_mp_nwp.setup(conductance='throat.conduit_hydraulic_conductance')
            Stokes_alg_mp_wp.run()
            Stokes_alg_mp_nwp.run()
            Keff_mp_wp = Stokes_alg_mp_wp.calc_effective_permeability(inlets=self.settings['flow_inlet'], outlets=self.settings['flow_outlet'])
            Keff_mp_nwp = Stokes_alg_mp_nwp.calc_effective_permeability(inlets=self.settings['flow_inlet'], outlets=self.settings['flow_outlet'])
            self.settings['relperm_wp'].append(Keff_mp_wp/self.settings['perm_wp'])
            self.settings['relperm_nwp'].append(Keff_mp_nwp/self.settings['perm_nwp'])
            self.project.purge_object(obj=Stokes_alg_mp_wp)
            self.project.purge_object(obj=Stokes_alg_mp_nwp)
        print('relperm_wp',self.settings['relperm_wp'])
        print('relperm_nwp',self.settings['relperm_nwp'])
        print('sat',self.settings['sat'])
            
            
                
    

