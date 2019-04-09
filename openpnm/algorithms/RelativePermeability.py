from openpnm.algorithms import GenericAlgorithm, StokesFlow, InvasionPercolation
from openpnm.utils import logging
from openpnm import models
import numpy as np
import openpnm
# import matplotlib.pyplot as plt
logger = logging.getLogger(__name__)


default_settings = {'pore_inv_seq': [],
                    'throat_inv_seq': [],
                    'points': 20,
                    'gh': 'throat.hydraulic_conductance',
                    'mode': 'strict',
                    'sat': [],
                    'inv_results': [],
                    'inlets': [],
                    'outlets': [],
                    'user_inlets': [],
                    'Pcap': [],
                    'pore_occ': [],
                    'throat_occ': [],
                    'auto': [],
                    'inv_phase': [],
                    'def_phase': [],
                    'auto_def': [],
                    'auto_inv': [],
                    'BP_1': [],
                    'BP_2': [],
                    'BC_1': [],
                    'BC_2': [],
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
        # Apply default settings
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        # Apply any settings received during initialization
        self.settings.update(settings)

    def setup(self, inv_phase=None, def_phase=None, points=None,
              sats=None,
              pore_inv_seq=None,
              throat_inv_seq=None,
              inlets=None,
              outlets=None,
              pore_occ=None,
              throat_occ=None,
              b1_pores=None,
              b2_pores=None,
              auto=None):
        self.settings['auto']=auto
        network = self.project.network
        if self.settings['auto'] is not False:
            
            # define phases (in this case we ignore inv_phase and def_phase, as they are defined for general algorithm)
            oil = openpnm.phases.GenericPhase(network=network, name='oil')
            water = openpnm.phases.GenericPhase(network=network, name='water')
            oil['pore.viscosity']=0.547
            oil['throat.contact_angle'] =110
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
            self.settings['auto_inv']= oil
            self.settings['auto_def']= water
            
            # define inlet/outlets
            self.settings['user_inlets']=False
            Hinlets = [network.pores(['top']), network.pores(['left']),
                          network.pores(['front'])]
            inlets=[]
            for i in range(len(Hinlets)):
                inlets.append([Hinlets[i][x] for x in range(0, len(Hinlets[i]), 2)])
            self.settings['inlets']=self.set_inlets(inlets)
            Houtlets = [network.pores(['bottom']), network.pores(['right']),
                           network.pores(['back'])]
            outlets=[]
            for i in range(len(Houtlets)):
                outlets.append([Houtlets[i][x] for x in range(0, len(Houtlets[i]), 2)])
            self.settings['outlets']=self.set_outlets(outlets)
            
            # define BC and B_pores
            self.settings['BP_1']=Hinlets
            self.settings['BP_2']=Houtlets
            # By default gets values not type. type is 'value' by default
            self.settings['BC_1']=[1,1,1]
            self.settings['BC_2']=[0,0,0]
            
            # Apply 3D Invasion Percolation
            inv=InvasionPercolation(phase=self.settings['auto_inv'],network=network)
            inv.setup(phase=self.settings['auto_inv'],
                      entry_pressure='throat.entry_pressure',
                      pore_volume='pore.volume',
                      throat_volume='throat.volume')
            pore_inv_seq={'0': [],'1': [],'2': []} # each element is a vector
            throat_inv_seq= {'0': [],'1': [],'2': []} # each element is a vector
            self.settings['sat'] = {'0': [],'1': [],'2': []} # each element is a vector, which is an approved saturation sequence for IP
            pore_occ = {'0': [],'1': [],'2': []} # each element(dimension) is a matrix
            throat_occ = {'0': [],'1': [],'2': []} # each element(dimension) is a matrix
            for i in range(len(self.settings['inlets'])): # for each dimension
                inv.set_inlets(pores=self.settings['inlets'][i])
                inv.run()
                pore_inv_seq[str(i)]=inv['pore.invasion_sequence']
                throat_inv_seq[str(i)]=inv['throat.invasion_sequence']
                Snwparr =  []
                Pcarr =  []
                Sarr=np.linspace(0,1,num=20)
                for Snw in Sarr:
                    #print('Snw is', Snw)
                    res1=inv.results(Snwp=Snw)
                    occ_ts=res1['throat.occupancy']
                    if np.any(occ_ts):
                        max_pthroat=np.max(self.settings['auto_inv']['throat.entry_pressure'][occ_ts])
                        Pcarr.append(max_pthroat)
                        Snwparr.append(Snw)
                self.settings['sat'][str(i)].append(Snwparr)
                c=-1
                for Sp in Snwparr:
                    c=c+1
                    res=inv.results(Sp)
                    pore_occ[str(i)].append(res['pore.occupancy'])
                    throat_occ[str(i)].append(res['throat.occupancy'])
            self.settings['pore_occ']=pore_occ
            self.settings['throat_occ']=throat_occ
            
    def run(self,network):
            # find Kx,Ky,Kz
            single_perms_water = [None,None,None] # each element is a scalar
            single_perms_oil = [None,None,None] # each element is a scalar
            [Da,Dl]=self.domain_l_a() # Da and Dl are vectors of 3 values for 3 dimensions
            for i in range(len(self.settings['inlets'])):
                BC1_pores = self.settings['BP_1'][i]
                BC2_pores = self.settings['BP_2'][i]
                #self.autorun()
                [da,dl]=[Da[i],Dl[i]]
                Stokes_alg_single_phase_water = StokesFlow(network=network, phase=self.settings['auto_def'])
                Stokes_alg_single_phase_water.setup(conductance='throat.hydraulic_conductance')
                Stokes_alg_single_phase_water._set_BC(pores=BC1_pores, bctype='value', bcvalues=self.settings['BC_1'][i])
                Stokes_alg_single_phase_water._set_BC(pores=BC2_pores, bctype='value', bcvalues=self.settings['BC_2'][i])
                Stokes_alg_single_phase_water.run()
#                single_perms_water[i] = Stokes_alg_single_phase_water.calc_effective_permeability(domain_area=da, domain_length=dl,inlets=BC1_pores, outlets=BC2_pores)
                single_perms_water[i] = Stokes_alg_single_phase_water.calc_effective_permeability(inlets=BC1_pores, outlets=BC2_pores)
                self.project.purge_object(obj=Stokes_alg_single_phase_water)
                Stokes_alg_single_phase_oil = StokesFlow(network=network, phase=self.settings['auto_inv'])
                Stokes_alg_single_phase_oil.setup(conductance='throat.hydraulic_conductance')
                Stokes_alg_single_phase_oil._set_BC(pores=BC1_pores, bctype='value', bcvalues=self.settings['BC_1'][i])
                Stokes_alg_single_phase_oil._set_BC(pores=BC2_pores, bctype='value', bcvalues=self.settings['BC_2'][i])
                Stokes_alg_single_phase_oil.run()
                single_perms_oil[i] = Stokes_alg_single_phase_oil.calc_effective_permeability(inlets=BC1_pores, outlets=BC2_pores)
#                single_perms_oil[i] = Stokes_alg_single_phase_oil.calc_effective_permeability(domain_area=da, domain_length=dl,inlets=BC1_pores, outlets=BC2_pores)
                self.project.purge_object(obj=Stokes_alg_single_phase_oil)
            # find Krx,Kry,Krz
            perm_water = {'0': [],'1': [],'2': []} # each element is a vector
            perm_oil = {'0': [],'1': [],'2': []} # each element is a vector
            # print(self.settings['sat'])
            for i in range(len(self.settings['sat'])):
                S_vec=self.settings['sat'][str(i)]
                [da,dl]=[Da[i],Dl[i]]
                for j in range(len(S_vec[0])):
                    self.update_phase_and_phys(self.settings['pore_occ'][str(i)][j], self.settings['throat_occ'][str(i)][j])
                    #print('sat is equal to', S_vec[0][j])
                    Stokes_alg_multi_phase_water = StokesFlow(network=network,phase=self.settings['auto_def'])
                    Stokes_alg_multi_phase_water.setup(conductance='throat.conduit_hydraulic_conductance')
                    Stokes_alg_multi_phase_water.set_value_BC(values=self.settings['BC_1'][i], pores=self.settings['BP_1'][i])
                    Stokes_alg_multi_phase_water.set_value_BC(values=self.settings['BC_2'][i], pores=self.settings['BP_2'][i])
                    #oil
                    Stokes_alg_multi_phase_oil = StokesFlow(network=network,phase=self.settings['auto_inv'])
                    Stokes_alg_multi_phase_oil.setup(conductance='throat.conduit_hydraulic_conductance')
                    Stokes_alg_multi_phase_oil.set_value_BC(values=self.settings['BC_1'][i], pores=self.settings['BP_1'][i])
                    Stokes_alg_multi_phase_oil.set_value_BC(values=self.settings['BC_2'][i], pores=self.settings['BP_2'][i])
                    # Run Multiphase algs
                    Stokes_alg_multi_phase_water.run()
                    Stokes_alg_multi_phase_oil.run()
                    effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_effective_permeability(inlets=BC1_pores, outlets=BC2_pores)
                    effective_permeability_oil_multi = Stokes_alg_multi_phase_oil.calc_effective_permeability(inlets=BC1_pores, outlets=BC2_pores)
                    #effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_effective_permeability(domain_area=da, domain_length=dl)
                    #effective_permeability_oil_multi = Stokes_alg_multi_phase_oil.calc_effective_permeability(domain_area=da, domain_length=dl)
                    relative_eff_perm_water = effective_permeability_water_multi/single_perms_water[i]
                    relative_eff_perm_oil = effective_permeability_oil_multi/single_perms_oil[i]
                    perm_water[str(i)].append(relative_eff_perm_water)
                    perm_oil[str(i)].append(relative_eff_perm_oil)
                    self.project.purge_object(obj=Stokes_alg_multi_phase_water)
                    self.project.purge_object(obj=Stokes_alg_multi_phase_oil)
                    #print('krw is for direction',i, 'kr', perm_water[str(i)])
#    def result(self,perm_water,perm_oil):
            return [self.settings['sat'],perm_water, perm_oil]
            
    def update_phase_and_phys(self, pore_occ, throat_occ):
        # inv_p=self.project.phases(self.settings['inv_phase'])
        # def_p=self.project.phases(self.settings['def_phase'])
        if self.settings['auto'] is not False:
            self.settings['auto_inv']['pore.occupancy'] = pore_occ
            self.settings['auto_def']['pore.occupancy'] = 1-pore_occ
            self.settings['auto_inv']['throat.occupancy'] = throat_occ
            self.settings['auto_def']['throat.occupancy'] = 1-throat_occ
            mode=self.settings['mode']
            self.settings['auto_def'].add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)
            self.settings['auto_inv'].add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)

    def domain_l_a(self):
        # for now we end up with defining default domain length and area
        if self.settings['user_inlets'] is not True:
            da=[]
            dl=[]
            network = self.project.network
            [amax, bmax, cmax] = np.max(network['pore.coords'], axis=0)
            [amin, bmin, cmin] = np.min(network['pore.coords'], axis=0)
            lx = amax-amin
            ly = bmax-bmin
            lz = cmax-cmin
            options = {0: self.top_b(lx, ly, lz),
                       1: self.left_r(lx, ly, lz),
                       2: self.front_b(lx, ly, lz)}
            for i in range(len(options)):
                [Da, Dl]=options[i]
                da.append(Da)
                dl.append(Dl)
        return [da, dl]
    
    def top_b(self, lx, ly, lz):
        da = lx*ly
        dl = lz
        res_2=[da, dl]
        return res_2

    def left_r(self, lx, ly, lz):
        da = lx*lz
        dl = ly
        res_2=[da, dl]
        return res_2

    def front_b(self, lx, ly, lz):
        da = ly*lz
        dl = lx
        res_2=[da, dl]
        return res_2
    
    def set_inlets(self, pores):
        r"""
        """
        pores_in=[]
        for inlet_num in range(len(pores)):
            self['pore.inlets'] = False
            self['pore.inlets'][pores[inlet_num]] = True
            pores_in.append(self['pore.inlets'])
        # print('inlets are',pores_in)
        self['pore.inlets'] = False
        return pores_in

    def set_outlets(self, pores):
        r"""
        """
        pores_out=[]
        for outlet_num in range(len(pores)):
            self['pore.outlets'] = False
            self['pore.outlets'][pores[outlet_num]] = True
            pores_out.append(self['pore.outlets'])
        # print('outlets are', pores_out)
        self['pore.outlets']= False
        return pores_out
        