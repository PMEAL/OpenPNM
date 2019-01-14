from openpnm.algorithms import GenericAlgorithm, StokesFlow, InvasionPercolation
from openpnm.utils import logging
from openpnm import models
import numpy as np
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
              pore_inv_seq=None,
              throat_inv_seq=None,
              inlets=None,
              outlets=None):
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
        network = self.project.network
        if inlets:
            self.settings['user_inlets']=True
            self.settings['inlets']=self.set_inlets(inlets)
        else:
            self.settings['user_inlets']=False
            pores=[]
            inlets = [network.pores(['top']), network.pores(['front']),
                      network.pores(['left'])]
            for inlet_num in range(len(inlets)):
                used_inlets = [inlets[x] for x in range(0, len(inlets[inlet_num]), 2)]
                pores.append[used_inlets]
            self.settings['inlets']=self.set_inlets(pores)
            pores=[]
            outlets = [network.pores(['bottom']), network.pores(['back']),
                      network.pores(['right'])]
            for outlet_num in range(len(outlets)):
                used_outlets = [outlets[x] for x in range(0, len(outlets[outlet_num]), 2)]
                pores.append[used_outlets]
            self.settings['outlets']=self.set_outlets(pores)
            #######################calc domain l and a
        if outlets:
            self.settings['outlets']=self.set_outlets(outlets)
        if inv_phase:
            self.settings['inv_phase'] = inv_phase.name
        if def_phase:
            self.settings['def_phase'] = def_phase.name
        if points:
            self.settings['points'] = points
        if pore_inv_seq:
            self.settings['pore_inv_seq'] = pore_inv_seq
        else:
            for inlet_num in range(len(pores)):
                inv_seq=self.IP(inlets=self.settings['inlets'][inlet_num],
                                outlets=self.settings['outlets'][inlet_num],
                                sim_num=inlet_num)
                self.settings['pore_inv_seq'].append[inv_seq[inlet_num][0]]
                self.settings['thorat_inv_seq'].append[inv_seq[inlet_num][1]]
                # the following lines are ignored assumming that once we have
                # the pore_inv_seq we also have throat_inv_seq as long as
                # both of them are produced as a result of running IP.
        if throat_inv_seq:
            self.settings['thorat_inv_seq'] = throat_inv_seq
#       else:
#            self.IP()

    def IP(self, inlets=None, outlets=None, sim_num=1):
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
            occ_ts=res1['throat.occupancy']
            if np.any(occ_ts):
                max_pthroat=np.max(phase['throat.entry_pressure'][occ_ts])
                Pcarr.append(max_pthroat)
                Snwparr.append(Snw)
#        plt.figure(1)
#        y=np.array(Pcarr[:])
#        x=1.0-np.array(Snwparr[:])
#        plt.xticks(np.arange(x.min(), x.max(), 0.05))
#        plt.yticks(np.arange(y.min(), y.max(),0.1))
#        plt.plot(x, y)
#        plt.xlabel('Invading Phase Saturation')
#        plt.ylabel('Capillary Pressure')
#        plt.grid(True)
        self.settings['sat']=np.array(Snwparr[:])
        for Sp in self.settings['sat']:
            self.settings['inv_results'][sim_num].append(inv.results(Sp))
        inv_seq=[inv['pore.invasion_sequence'], inv['throat.invasion_sequence']]
        # assumming the last array is corresponding to the Capillary pressure
        # we did not include saturations in the results
        # saturations can be taken from self.settings['sat']
        self.settings['inv_results'][sim_num].append(Pcarr)
        return inv_seq
    
    def domain_l_a(self):
        # for now we end up with defining default domain length and area
        if self.settings['user_inlets'] =='False':
            da=[]
            dl=[]
            network = self.project.network
            [amax, bmax, cmax] = np.max(network['pore.coords'], axis=0)
            [amin, bmin, cmin] = np.min(network['pore.coords'], axis=0)
            lx = amax-amin
            ly = bmax-bmin
            lz = cmax-cmin
            options = {0 : self.top_b(lx,ly,lz), 1 : self.left_r(lx,ly,lz), 2 : self.front_b(lx,ly,lz)}
            for i in range(len(options)):
                [da[i], dl[i]]=options[i]
        return [da,dl]

    def set_inlets(self, pores):
        r"""
        """
        for inlet_num in range(len(pores)):
            self['pore.inlets'] = False
            self['pore.inlets'][pores[inlet_num]] = True
            self.settings['inlets'][inlet_num]=self['pore.inlets'][pores[inlet_num]]

    def set_outlets(self, pores):
        r"""
        """
        for outlet_num in range(len(pores)):
            self['pore.outlets'] = False
            self['pore.outlets'][pores[outlet_num]] = True
            self.settimgs['outlets'][outlet_num]=self['pore.outlets'][pores[outlet_num]]

    def update_phase_and_phys(self, results):
        inv_p=self.project.phases(self.settings['inv_phase'])
        def_p=self.project.phases(self.settings['def_phase'])
        inv_p['pore.occupancy'] = results['pore.occupancy'] #################################### check
        def_p['pore.occupancy'] = 1-results['pore.occupancy']
        inv_p['throat.occupancy'] = results['throat.occupancy']
        def_p['throat.occupancy'] = 1-results['throat.occupancy']
        # adding multiphase conductances
        mode=self.settings['mode']
        def_p.add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)
        inv_p.add_model(model=models.physics.multiphase.conduit_conductance,
                        propname='throat.conduit_hydraulic_conductance',
                        throat_conductance='throat.hydraulic_conductance',
                        mode=mode)

    def top_b(self, lx, ly, lz):
        da = lx*ly
        dl = lz
        res_2=[da, dl]
        return res_2
    def left_r(lx,ly,lz):
        da = lx*lz
        dl = ly
        res_2=[da,dl]
        return res_2

    def front_b(lx,ly,lz):
        da = ly*lz
        dl = lx
        res_2=[da,dl]
        return res_2

    def run(self):
        r"""
        """
        Results = {'k_inv': [], 'k_def': [], 'K_rel_inv': [], 'K_rel_def': []}
        # Retrieve phase and network
        K_rel_def=[]
        K_rel_inv=[]
        for i in range(len(self.settings['inlets'])):
            K_rel_def[i]=[]
            K_rel_inv[i]=[]
        network = self.project.network
        # K_def=1
        # K_inv=[]
        # # apply single phase flow
        inlets=self.settings['inlets']
        outlets=self.settings['outlets']
        for bound_num in range(len(self.settings['inlets'])):
            # Run Single phase algs effective properties
            # Kw
            St_def = StokesFlow(network=network,
                                phase=self.project.phases(self.settings['def_phase']))
            St_def.setup(conductance=self.settings['gh'])
            St_def._set_BC(pores=inlets[bound_num], bctype='value', bcvalues=1)
            St_def._set_BC(pores=outlets[bound_num], bctype='value', bcvalues=0)
            St_def.run()
            if self.settings['user_inlets']==False:
                [da,dl]=self.domain_l_a()
                K_def = St_def.calc_effective_permeability(domain_area=da[bound_num],
                                                           domain_length=dl[bound_num],
                                                           inlets=inlets[bound_num],
                                                           outlets=outlets[bound_num])
            else:
                K_def = St_def.calc_effective_permeability(inlets=inlets[bound_num],
                                                           outlets=outlets[bound_num])
            # proj.purge_object(obj=St_def)
            # Ko
            Results['k_def'].append(K_def)
            St_inv = StokesFlow(network=network,
                                phase=self.project.phases(self.settings['inv_phase']))
            St_inv.setup(conductance=self.settings['gh'])
            St_inv._set_BC(pores=inlets[bound_num], bctype='value', bcvalues=1)
            St_inv._set_BC(pores=outlets[bound_num], bctype='value', bcvalues=0)
            St_inv.run()
            if self.settings['user_inlets']==False:
                K_inv = St_inv.calc_effective_permeability(domain_area=da[bound_num],
                                                           domain_length=dl[bound_num],
                                                      inlets=inlets[bound_num],
                                                           outlets=outlets[bound_num])     
            else:
                K_inv = St_inv.calc_effective_permeability(inlets=inlets[bound_num],
                                                           outlets=outlets[bound_num])
            Results['k_inv'].append(K_inv)
            self.project.purge_object(obj=St_inv)
        # apply two phase effective perm calculation
        for bound_num in range(len(self.settings['inlets'])):
            cn=-1
            for Sp in self.settings['sat']:
                cn=cn+1
                self.update_phase_and_phys(self.settings['inv_results'][bound_num][cn])
                # print('sat is equal to', Sp)
                inv_p=self.project.phases(self.settings['inv_phase'])
                def_p=self.project.phases(self.settings['def_phase'])
                # water
                St_def_tp= StokesFlow(network=network,
                                      phase=def_p)
                St_def_tp.setup(conductance='throat.conduit_hydraulic_conductance')
                St_def_tp.set_value_BC(pores=inlets[bound_num], values=1)
                St_def_tp.set_value_BC(pores=outlets[bound_num], values=0)
                # oil
                St_inv_tp = StokesFlow(network=network,
                                       phase=inv_p)
                St_inv_tp.setup(conductance='throat.conduit_hydraulic_conductance')
                St_inv_tp.set_value_BC(pores=inlets[bound_num], values=1)
                St_inv_tp.set_value_BC(pores=outlets[bound_num], values=0)
                # Run Multiphase algs
                St_def_tp.run()
                St_inv_tp.run()
                if self.settings['user_inlets']==False:
                    K_def_tp = St_def_tp.calc_effective_permeability(domain_area=da[bound_num],
                                                                 domain_length=dl[bound_num],
                                                                 inlets=inlets[bound_num],
                                                                 outlets=outlets[bound_num])
                    K_inv_tp = St_inv_tp.calc_effective_permeability(domain_area=da[bound_num],
                                                                 domain_length=dl[bound_num],
                                                                 inlets=inlets[bound_num],
                                                                 outlets=outlets[bound_num])
                else:
                    K_def_tp = St_def_tp.calc_effective_permeability(inlets=inlets[bound_num],
                                                                 outlets=outlets[bound_num])
                    K_inv_tp = St_inv_tp.calc_effective_permeability(inlets=inlets[bound_num],
                                                                 outlets=outlets[bound_num])
                krel_def =K_def_tp/Results['k_def'][bound_num]
                krel_inv= K_inv_tp /Results['k_inv'][bound_num]
                K_rel_def[bound_num].append(krel_def)
                K_rel_inv[bound_num].append(krel_inv)
                self.project.purge_object(obj=St_def_tp)
                self.project.purge_object(obj=St_inv_tp)
            Results['K_rel_inv']=K_rel_inv
            Results['K_rel_def']=K_rel_def
        return Results
