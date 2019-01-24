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
              auto=None):
        r"""
         Set up the required parameters for the algorithm

        Parameters
        ----------
        inv_phase: OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        def_phase: OpenPNM Phase object
            The phase which will be withdrawn from the Network during the
            injection of invading phase.

        points : Scalar, integer
            Number of saturation points for which the relative permeabilities
            are calculated.

        pore_inv_seq : list of boolean array
            Arrays of the pores where invaded pores are assigned as True. If
            the data is not provided by user, the setup will call IP method
            in order to get the results by implementing Invasion Percolation.

        throat_inv_seq: list of boolean array
            Arrays of the throats where invaded throats are assigned as True.
            If the data is not provided by user, the setup will call IP method
            in order to get the results by implementing Invasion Percolation.
        inlets: list of inlets (list of pores)
            A list of lists, in which each element is a list of pores as an
            inlet defined for each simulation.
        outlets: list of inlets (list of pores)
            A list of list,l in which each element is a list of pores as an
            outlet defined for each simulation.
        """
        self.settings['auto']=auto
        network = self.project.network
        if inlets:
            self.settings['user_inlets']=True
            # self.settings['inlets']=self.set_inlets(inlets)
            self.settings['inlets']=inlets
            print('1line', self.settings['inlets'])
        else:
            print('not specified')
            self.settings['user_inlets']=False
            Hinlets = [network.pores(['top']), network.pores(['top']),
                          network.pores(['top'])]
            inlets=[]
            for i in range(len(Hinlets)):
                inlets.append([Hinlets[i][x] for x in range(0, len(Hinlets[i]), 2)])
            self.settings['inlets']=self.set_inlets(inlets)
            outlets = [network.pores(['bottom']), network.pores(['bottom']),
                           network.pores(['bottom'])]
            self.settings['outlets']=self.set_outlets(outlets)   # should
            # be changed with output of outlet calc because we need pore
            # itself (object) not just the indice
            # print('inlets are',self.settings['inlets'])
            # print('outlets are',self.settings['outlets'])
#        if outlets.any():
#            # self.settings['outlets']=self.set_outlets(outlets)
#            self.settings['outlets']= outlets
#            print('1line', self.settings['outlets'])
        if inv_phase:
            self.settings['inv_phase'] = inv_phase.name
        else:
            oil = openpnm.phases.GenericPhase(network=network)
            water = openpnm.phases.GenericPhase(network=network)
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
#            phys_water= openpnm.physics.GenericPhysics(network=network, phase=water, geometry=self.project.geometry)
#            phys_oil = openpnm.physics.GenericPhysics(network=network, phase=oil, geometry=self.project.geometry)
            mod = openpnm.models.physics.hydraulic_conductance.hagen_poiseuille
            oil.add_model(propname='throat.hydraulic_conductance',
                              model=mod)
            oil.add_model(propname='throat.entry_pressure',
                              model=openpnm.models.physics.capillary_pressure.washburn)
            water.add_model(propname='throat.hydraulic_conductance',
                              model=mod)
            water.add_model(propname='throat.entry_pressure',
                              model=openpnm.models.physics.capillary_pressure.washburn)
            self.settings['inv_phase']= oil.name
            self.settings['def_phase']= water.name
        if points:
            self.settings['points'] = points
        else:
            self.settings['points']= 15
        if pore_inv_seq:
            self.settings['sat']=sats
            self.settings['pore_inv_seq'] = pore_inv_seq
            self.settings['throat_inv_seq'] = throat_inv_seq
            self.settings['pore_occ']=pore_occ
            self.settings['throat_occ']=throat_occ
            print('1lineis', self.settings['pore_inv_seq'])
            print('1lineis', self.settings['throat_inv_seq'])
            print('1lineis', self.settings['pore_occ'])
            print('1lineis', self.settings['throat_occ'])
        elif self.settings['user_inlets'] is not True:   # ## will be uncommented later on
            inv=InvasionPercolation(phase=self.project.phases(self.settings['inv_phase']),network=network)
            inv.setup(phase=self.project.phases(self.settings['inv_phase']),
                      entry_pressure='throat.entry_pressure',
                      pore_volume='pore.volume',
                      throat_volume='throat.volume')
            pore_inv_seq=[]
            throat_inv_seq=[]
            sat = {'0': [],'1': [],'2': []}
            perm_water = {'0': [],'1': [],'2': []}
            perm_oil = {'0': [],'1': [],'2': []}
            pore_occ = {'0': [],'1': [],'2': []}
            throat_occ = {'0': [],'1': [],'2': []}
            single_perms_water = [None,None,None]
            single_perms_oil = [None,None,None]
            [Da,Dl]=self.domain_l_a()
            for i in range(len(self.settings['inlets'])):
                inv.set_inlets(pores=self.settings['inlets'][i])
                inv.run()
                pore_inv_seq=inv['pore.invasion_sequence']
                throat_inv_seq=inv['throat.invasion_sequence']
                Snwparr =  []
                Pcarr =  []
                Sarr=np.linspace(0,1,num=15)
                for Snw in Sarr:
                    res1=inv.results(Snwp=Snw)
                    occ_ts=res1['throat.occupancy']
                    if np.any(occ_ts):
                        max_pthroat=np.max(oil['throat.entry_pressure'][occ_ts])
                        Pcarr.append(max_pthroat)
                        Snwparr.append(Snw)
                sat[str(i)].append(Snwparr)
                c=-1
                for Sp in Snwparr:
                    c=c+1
                    res=inv.results(Sp)
                    pore_occ[str(i)].append(res['pore.occupancy'])
                    throat_occ[str(i)].append(res['throat.occupancy'])
            for i in range(len(self.settings['inlets'])):
                BC1_pores = Hinlets[i]
                BC2_pores = outlets[i]
                #self.autorun()
                [da,dl]=[Da[i],Dl[i]]
                Stokes_alg_single_phase_water = StokesFlow(network=network, phase=water)
                Stokes_alg_single_phase_water.setup(conductance='throat.hydraulic_conductance')
                Stokes_alg_single_phase_water._set_BC(pores=BC1_pores, bctype='value', bcvalues=100000)
                Stokes_alg_single_phase_water._set_BC(pores=BC2_pores, bctype='value', bcvalues=1000)
                Stokes_alg_single_phase_water.run()
                single_perms_water[i] = Stokes_alg_single_phase_water.calc_effective_permeability(domain_area=da, domain_length=dl,inlets=BC1_pores, outlets=BC2_pores)
                self.project.purge_object(obj=Stokes_alg_single_phase_water)
                Stokes_alg_single_phase_oil = StokesFlow(network=network, phase=oil)
                Stokes_alg_single_phase_oil.setup(conductance='throat.hydraulic_conductance')
                Stokes_alg_single_phase_oil._set_BC(pores=BC1_pores, bctype='value', bcvalues=10000)
                Stokes_alg_single_phase_oil._set_BC(pores=BC2_pores, bctype='value', bcvalues=1000)
                Stokes_alg_single_phase_oil.run()
                single_perms_oil[i] = Stokes_alg_single_phase_oil.calc_effective_permeability(domain_area=da, domain_length=dl,inlets=BC1_pores, outlets=BC2_pores)
                self.project.purge_object(obj=Stokes_alg_single_phase_oil)
                c=-1
                for Sp in range(len(sat[str(i)])):
                    c=c+1
                    self.update_phase_and_phys(c, i, pore_occ, throat_occ)
                    print('sat is equal to', Sp)
                    Stokes_alg_multi_phase_water = StokesFlow(network=network,phase=water)
                    Stokes_alg_multi_phase_water.setup(conductance='throat.conduit_hydraulic_conductance')
                    Stokes_alg_multi_phase_water.set_value_BC(values=100000, pores=BC1_pores)
                    Stokes_alg_multi_phase_water.set_value_BC(values=1000, pores=BC2_pores)
                    #oil
                    Stokes_alg_multi_phase_oil = StokesFlow(network=network,phase=oil)
                    Stokes_alg_multi_phase_oil.setup(conductance='throat.conduit_hydraulic_conductance')
                    Stokes_alg_multi_phase_oil.set_value_BC(values=100000, pores=BC1_pores)
                    Stokes_alg_multi_phase_oil.set_value_BC(values=1000, pores=BC2_pores)
                    # Run Multiphase algs
                    Stokes_alg_multi_phase_water.run()
                    Stokes_alg_multi_phase_oil.run()
                    effective_permeability_water_multi = Stokes_alg_multi_phase_water.calc_effective_permeability(domain_area=da, domain_length=dl)
                    effective_permeability_oil_multi = Stokes_alg_multi_phase_oil.calc_effective_permeability(domain_area=da, domain_length=dl)
                    relative_eff_perm_water = effective_permeability_water_multi/single_perms_water[i]
                    relative_eff_perm_oil = effective_permeability_oil_multi/single_perms_oil[i]
                    perm_water[str(i)].append(relative_eff_perm_water)
                    perm_oil[str(i)].append(relative_eff_perm_oil)
                    self.project.purge_object(obj=Stokes_alg_multi_phase_water)
                    self.project.purge_object(obj=Stokes_alg_multi_phase_oil)
                
                
            
        else:
            print('not specified')
#            ins=self.settings['inlets']
#            pore_inv=[]
#            throat_inv=[]
#            pore_occ=[]
#            throat_occ=[]
#            sat=[]
#            for inlet_num in range(len(ins)):
#                results=self.IP(inlets=ins[inlet_num],
#                                sim_num=inlet_num)
#                pore_inv.append(results['pore_inv'])
#                throat_inv.append(results['throat_inv'])
#                pore_occ.append(results['pore_occ'])
#                throat_occ.append(results['throat_occ'])
#                sat.append(results['sat'])
#            self.settings['pore_occ']=pore_occ
#            self.settings['throat_occ']=throat_occ
#            self.settings['pore_inv_seq']=pore_inv
#            self.settings['thorat_inv_seq']=throat_inv
#            self.settings['sat']=sat
#            print('1',self.settings['pore_occ'])
#            print('2',self.settings['throat_occ'])
#            print('3',self.settings['pore_inv_seq'])
#            print('4',self.settings['thorat_inv_seq'])
#            print('5',self.settings['sat'])
# the following lines are ignored assumming that once we have
# the pore_inv_seq we also have throat_inv_seq as long as
# both of them are produced as a result of running IP.
#        if throat_inv_seq:
#            self.settings['thorat_inv_seq'] = throat_inv_seq
#       else:
#            self.IP()
        
    def IP(self, inlets=None, sim_num=1):
        network = self.project.network
        phase = self.project.phases(self.settings['inv_phase'])
        inv=InvasionPercolation(phase=phase, network=network, project=self.project)
        inv.setup(phase=phase, pore_volume='pore.volume', throat_volume='throat.volume')
        inv.set_inlets(pores=inlets)
        inv.run()
        Snwparr =  []
        Pcarr =  []
        Sarr=np.linspace(0, 1, num=self.settings['points'])
        for Snw in Sarr:
            res1=inv.results(Snwp=Snw)
            print('sw is  :', Snw, 'result is', res1 )
            occ_ts=res1['throat.occupancy']
            if np.any(occ_ts):
                max_pthroat=np.max(phase['throat.entry_pressure'][occ_ts])
                Pcarr.append(max_pthroat)
                Snwparr.append(Snw)
#        x=1.0-np.array(Snwparr[:])
#        plt.xticks(np.arange(x.min(), x.max(), 0.05))
#        plt.yticks(np.arange(y.min(), y.max(),0.1))
#        plt.plot(x, y)
#        plt.xlabel('Invading Phase Saturation')
#        plt.ylabel('Capillary Pressure')
#        plt.grid(True)
        sat=np.array(Snwparr[:])
        print('sat is', sat)
        p_occ=[]
        t_occ=[]
        for Sp in sat:
            # inv_res.append(inv.results(Sp)) # gives pore and throat occupancy at Sp
            res=inv.results(Sp)
            p_occ.append(res['pore.occupancy'])
            t_occ.append(res['throat.occupancy'])
        inv_seq=[inv['pore.invasion_sequence'], inv['throat.invasion_sequence']]
        results = {'pore_inv': [], 'throat_inv': [], 'pore_occ': [], 'throat_occ': [], 'sat': []}
        # 'pore_inv' states invasion sequence in pores in an IP run
        # 'throat_inv' states invasion sequence in throats in an IP run
        # 'pore_occ' states pore occupancies for each saturation
        # 'throat_occ' states throat occupancy for each saturation
        # assumming the last array is corresponding to the Capillary pressure
        # we did not include saturations in the results
        # saturations can be taken from self.settings['sat']
        # ## self.settings['inv_results'][sim_num].append(Pcarr)
        results['pore_inv']=inv_seq[0]
        results['throat_inv']=inv_seq[1]
        results['pore_occ']=p_occ
        results['throat_occ']=t_occ
        results['sat']=sat
        return results

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

    def update_phase_and_phys(self, satnum, invnum, pore_occ, throat_occ):
        # inv_p=self.project.phases(self.settings['inv_phase'])
        # def_p=self.project.phases(self.settings['def_phase'])
        self.project.phases(self.settings['inv_phase'])['pore.occupancy'] = pore_occ[str(invnum)][satnum]
        self.project.phases(self.settings['def_phase'])['pore.occupancy'] = 1-pore_occ[str(invnum)][satnum]
        self.project.phases(self.settings['inv_phase'])['throat.occupancy'] = throat_occ[str(invnum)][satnum]
        self.project.phases(self.settings['def_phase'])['throat.occupancy'] = 1-throat_occ[str(invnum)][satnum]
        mode=self.settings['mode']
        self.project.phases(self.settings['def_phase']).add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)
        self.project.phases(self.settings['inv_phase']).add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)

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

    def run(self):
        r"""
        """
        Results = {'sat': [], 'k_inv': [], 'k_def': [], 'K_rel_inv': [], 'K_rel_def': []}
        print('sat length', len(Results['sat']))
        print('k_inv', len(Results['k_inv']))
        print('k_def', len(Results['k_def']))
        print('krelinv length', len(Results['K_rel_inv']))
        print('K_rel_def', len(Results['K_rel_def']))
        K_rel_def=[]
        K_rel_inv=[]
        network = self.project.network
        inlets=self.settings['inlets']
        outlets=self.settings['outlets']
        print('ins are', inlets)
        print('outs are', outlets)
        St_def = StokesFlow(network=network,
                            phase=self.project.phases(self.settings['def_phase']))
        St_def.setup(conductance='throat.hydraulic_conductance')
        St_def._set_BC(pores=inlets, bctype='value', bcvalues=1)
        St_def._set_BC(pores=outlets, bctype='value', bcvalues=0)
        St_def.run()
        if self.settings['user_inlets'] is not True:
            print('not specified')
#            [da, dl]=self.domain_l_a()
#            K_def = St_def.calc_effective_permeability(domain_area=da[bound_num],
#                                                       domain_length=dl[bound_num],
#                                                       inlets=inlets[bound_num],
#                                                       outlets=outlets[bound_num])
        else:
            K_def = St_def.calc_effective_permeability(inlets=inlets,
                                                          outlets=outlets)
        Results['k_def'].append(K_def)
        self.project.purge_object(obj=St_def)
        St_inv = StokesFlow(network=network,
                            phase=self.project.phases(self.settings['inv_phase']))
        St_inv.setup(conductance='throat.hydraulic_conductance')
        St_inv._set_BC(pores=inlets, bctype='value', bcvalues=1)
        St_inv._set_BC(pores=outlets, bctype='value', bcvalues=0)
        St_inv.run()
        if self.settings['user_inlets'] is not True:
            print('not specified')
#            K_inv = St_inv.calc_effective_permeability(domain_area=da[bound_num],
#                                                          domain_length=dl[bound_num],
#                                                          inlets=inlets[bound_num],
#                                                          outlets=outlets[bound_num])
        else:
            K_inv = St_inv.calc_effective_permeability(inlets=inlets,
                                                           outlets=outlets)
        Results['k_inv'].append(K_inv)
        self.project.purge_object(obj=St_inv)
        cn=-1
        print('saturation is',self.settings['sat'])
        for Sp in self.settings['sat']:
            cn=cn+1
            self.update_phase_and_phys(sat_num=cn)
            # inv_p=self.project.phases(self.settings['inv_phase'])
            # def_p=self.project.phases(self.settings['def_phase'])
            St_def_tp= StokesFlow(network=network,
                                      phase=self.project.phases(self.settings['def_phase']))
            St_def_tp.setup(conductance='throat.conduit_hydraulic_conductance')
            St_def_tp.set_value_BC(pores=inlets, values=1)
            St_def_tp.set_value_BC(pores=outlets, values=0)
            St_def_tp.run()
            St_inv_tp = StokesFlow(network=network,
                                       phase=self.project.phases(self.settings['inv_phase']))
            St_inv_tp.setup(conductance='throat.conduit_hydraulic_conductance')
            St_inv_tp.set_value_BC(pores=inlets, values=1)
            St_inv_tp.set_value_BC(pores=outlets, values=0)
            St_inv_tp.run()
            if self.settings['user_inlets'] is not True:
                print('not specified')
#                def_ef=St_def_tp.calc_effective_permeability
#                K_def_tp = def_ef(domain_area=da[bound_num],
#                                      domain_length=dl[bound_num],
#                                      inlets=inlets[bound_num],
#                                      outlets=outlets[bound_num])
#                inv_ef=St_inv_tp.calc_effective_permeability
#                K_inv_tp = inv_ef(domain_area=da[bound_num],
#                                      domain_length=dl[bound_num],
#                                      inlets=inlets[bound_num],
#                                      outlets=outlets[bound_num])
#                print('keffinv',K_inv_tp)
#                print('keffdef',K_def_tp)
            else:
                def_ef= St_def_tp.calc_effective_permeability
                K_def_tp = St_def_tp.calc_effective_permeability(inlets=inlets,
                                      outlets=outlets)
                inv_ef=St_inv_tp.calc_effective_permeability
                K_inv_tp = St_inv_tp.calc_effective_permeability(inlets=inlets,
                                      outlets=outlets)
            krel_def =K_def_tp/Results['k_def']
            krel_inv= K_inv_tp /Results['k_inv']
            K_rel_def.append(krel_def)
            K_rel_inv.append(krel_inv)
            self.project.purge_object(obj=St_def_tp)
            self.project.purge_object(obj=St_inv_tp)
        Results['K_rel_inv']=K_rel_inv
        Results['K_rel_def']=K_rel_def
        Results['sat']=self.settings['sat']
        return Results
