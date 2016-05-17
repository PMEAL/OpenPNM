# -*- coding: utf-8 -*-
"""
===============================================================================
module __ViscousDrainage__: Viscous fluid flow with capillary pressure
===============================================================================

"""
import scipy as sp
from OpenPNM.Algorithms import GenericLinearTransport
from OpenPNM.Algorithms import Drainage as _drainage
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class ViscousDrainage(GenericLinearTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous drainage
    taking into account capillary pressures.

    This class adds several functions from __Drainage__ nescesary for the
    simultation.

    References
    ----------
    .. [1] Ferer, M., Bromnhal, G.S., Duane, H.S.
           Pore-level modeling of immiscible drainage: validation in the
           invasion percolation and DLA limits. Physica A319, 11-35 (2003)
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        logger.info('Create ' + self.__class__.__name__ + ' Object')

    def setup(self,
              wetting_phase,
              nonwetting_phase,
              injection_rate = None,
              conductance='hydraulic_conductance',
              entry_pressure='throat.capillary_pressure',
              quantity='pressure',
              pore_volume='pore.volume',
              throat_volume='throat.volume',
              super_pore_conductance=None,
              sat_tol=1.0E-6,
              max_steps=1E5, **params):
        r"""
        This setup provides the initial requirements for the solver setup
        and additional parameters for the drainage simulation.

        Parameters
        ----------
        nonwetting_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        wetting_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the defending
            fluid.

        injection_rate : float
            The bulk rate in m^3/sec the invading phase is injected into the
            network through the inlets. Used to set the Nuemann_group BC
            for the inlet.

        conductance : string (optional)
            The dictionary key on the Physics object where the throat conductance
            values are found.  The default is 'hydraulic_conductance'.

        entry_pressure : string (optional)
            The dictionary key on the Phase object where the throat entry
            pressure values can be found.  The default is
            'throat.capillary_pressure'.

        pore_volume and throat_volume : string (optional)
            The dictionary key on the Geometry object where the pore or throat
            volume data is located.  The defaults is 'pore.volume' and
            'throat.volume'.
        sat_tol : float (optional)
            Sets the maximum or minimum saturation value for a pore or
            throat to be rounded up to 1.0 or down to 0.0 during the simulation.
        """
        logger.info('Setup ' + self.__class__.__name__)
        #
        if (injection_rate is None):
            raise Exception('Error - injection rate must be specified')
        #
        if (nonwetting_phase is None):
            raise Exception('Error - Non-wetting phase must be specified')
        #
        if (wetting_phase is None):
            raise Exception('Error - Wetting phase must be specified')
        #
        # setting argument based values
        self._inv_phase = nonwetting_phase
        self._def_phase = wetting_phase
        self._inj_rate = injection_rate
        self['throat.entry_pressure'] = nonwetting_phase[entry_pressure]
        self._max_pc = nonwetting_phase[entry_pressure]
        self._throat_volume = throat_volume
        self._pore_volume = pore_volume
        self._sat_tol = sat_tol
        self._max_steps = max_steps
        # initializing simulation variables
        self._net_vol  = sp.sum(self._net['pore.volume'])
        self._net_vol += sp.sum(self._net['throat.volume'])
        self._trapping = True
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        # fraction of invading phase occupying pore/throat
        self._pore_inv_frac = sp.zeros(self.Np,dtype=float)
        self._throat_inv_frac = sp.zeros(self.Nt,dtype=float)
        # used in advance interface to specify how to change saturation fraction
        self._throat_sup_fact = sp.ones(self.Nt,dtype=float)*-1.0
        # bool identifying pores/throats that are fully invaded
        self._pore_invaded = sp.zeros(self.Np,dtype=bool)
        self._throat_invaded = sp.zeros(self.Nt,dtype=bool)
        # bool idenifying pores/throats that are being invaded
        self._pore_contested = sp.zeros(self.Np,dtype=bool)
        self._throat_contested = sp.zeros(self.Nt,dtype=bool)
        # menisci positions are set relative to lower indexed pore.
        self._menisci = [[] for i in range(self.Nt)]
        # function to scale capillary pressure based on miniscus location
        self._pc_func = lambda x: sp.sin(sp.pi * x)
        #
        self._log_file = open('VD-Log.txt','w');
        #
        super().setup(conductance=conductance, quantity=quantity,
                      super_pore_conductance=super_pore_conductance)

    def set_inlets(self, pores=None, mode='add'):
        r"""
        Sets inlet pores as well as inital menicsi to start simulation.
        """
        _drainage.set_inlets(self,pores=pores,mode=mode)
        #
        # inlets start off as fully invaded pores
        inlets = sp.where(self['pore.inlets'])[0]
        self._pore_inv_frac[inlets] = 1.0
        self._pore_invaded[inlets] = True
        #
        # throats between two inlets are set as filled to prevent initial
        # fluid plugs
        filled_throats = self._net.find_neighbor_throats(pores=inlets,
                                                         mode='intersection')
        self._throat_inv_frac[filled_throats] = 1.0
        self._throat_invaded[filled_throats] = True
        #
        # adding menisci
        interface = self._net.find_neighbor_throats(pores=inlets,
                                                    mode='not_intersection')
        self._throat_contested[interface] = True
        self._message('Setting menisci for throats on the inlet pores')
        self._q = sp.zeros(self.Nt)
        for th in interface:
            pore1,pore2 = self._net['throat.conns'][th]
            if (self._pore_invaded[pore1]):
                self._set_menisci(pore1,[th])
            else:
                self._set_menisci(pore2,[th])
        self._message('')

    def set_outlets(self, pores=None, mode='add'):
        _drainage.set_outlets(self,pores=pores,mode=mode)

    def run(self,**kwargs):
        #
        # Ensure inlets are set
        if sp.size(self['pore.inlets']) == 0:
            raise Exception('Inlet pores have not been specified')

        # Ensure outlet pores are set
        if sp.size(self['pore.outlets']) == 0:
            raise Exception('Outlet pores have not been specified')
        #
        # setting initial time to maintain mass balance if inlet pores and
        # throats have a non-zero volume
        self._total_time = 0.0
        tot_vol = sp.sum(sp.multiply(self._net['pore.volume'],self._pore_inv_frac))
        tot_vol += sp.sum(sp.multiply(self._net['throat.volume'],self._throat_inv_frac))
        tot_sat = tot_vol/self._net_vol
        self._total_time = tot_vol/self._inj_rate
        self._message('Initial Saturation of Invading Phase: ',tot_sat)
        self._message('')
        #
        # beginning simulation
        self._gdef = sp.copy(self['throat.conductance'])
        self._total_inv_out = 0.0
        self._total_def_out = 0.0
        self._do_outer_iteration_stage(**kwargs)
        #
        self._log_file.close()

    def _do_outer_iteration_stage(self, **kwargs):
        r"""
        This calls the solve method in the algorithm.
        Handles the tracking and movement of phases throat the network.
        """
        #
        self._inlets = sp.where(self['pore.inlets'])[0]
        self._outlets = sp.where(self['pore.outlets'])[0]
        #
        self.set_boundary_conditions(bctype='Dirichlet',
                                       mode='overwrite',
                                       bcvalue=0.0,
                                       pores=self._outlets)
        self.set_boundary_conditions(bctype='Neumann_group',
                                       mode='merge',
                                       bcvalue=self._inj_rate,
                                       pores=self._inlets)
        #
        break_through_time = -1
        break_through_steps = 0
        self._i = 0
        while True:
            self._modify_conductance()
            self._update_RHS_pcap_data()
            self.A = self._build_coefficient_matrix()
            self.b = self._build_RHS_matrix(self._net.pores(),self._RHS_pcap_data)
            self.b = sp.negative(self.b)
            self.solve(iterative_solver='cg')
            dt = self._calculate_dt()
            #
            self._message('Time Step: ',self._i,' size: {:0.3E} '.format(dt))
            #pres = ['{:0.3f}'.format(p) for p in self['pore.pressure']]
            #self._message('Pore Pressures: {}'.format(','.join(pres)))
            #
            self._advance_interface(dt)
            self._total_time += dt
            self._calc_fluid_out(dt)
            self._print_step_stats(self._i,dt)
            #
            test = sp.where(self._pore_inv_frac[self._outlets] > 1-self._sat_tol)[0]
            if (sp.size(test) > 0 and break_through_time < 0):
                break_through_time = self._total_time
                break_through_steps = self._i
                #break
            #
            if (self._i > self._max_steps):
                break
            #
            self._i += 1
        #
        # checking overall mass balance
        q_inj = self._total_time * self._inj_rate
        tot_vol = sp.sum(sp.multiply(self._net['pore.volume'],self._pore_inv_frac))
        tot_vol += sp.sum(sp.multiply(self._net['throat.volume'],self._throat_inv_frac))
        tot_sat = tot_vol/self._net_vol
        mass_bal = (q_inj - tot_vol - self._total_inv_out)/self._net_vol
        #
        self._message('Total Simulation Time Until Break Through: ',break_through_time,' Steps:',break_through_steps)
        self._message('Total Volume: ',tot_vol)
        self._message('Total Inv Fluid Out: ',self._total_inv_out)
        self._message('Total saturation: ',tot_sat)
        self._message('Total injection: ',q_inj)
        self._message('Mass Difference / Total Vol: {:15.9e}'.format(mass_bal))

    def _modify_conductance(self):
        r"""
        Updates the conductance based on the viscosity ratio of the fluids
        and the fractional occupancy of the throats.
        """
        #
        for th in self._net.throats():
            pores = self._net['throat.conns'][th]
            dvisc = sp.average(self._def_phase['pore.viscosity'][pores])
            ivisc = sp.average(self._inv_phase['pore.viscosity'][pores])
            M = dvisc/ivisc
            f = self._throat_inv_frac[th]
            #
            f = 1 - f + f*M
            self['throat.conductance'][th] = f * self._gdef[th]

    def _update_RHS_pcap_data(self):
        r"""
        Adds f * g * pcap to RHS for pores containing menisci
        """
        self._RHS_pcap_data = sp.zeros(self.Np,dtype=float)
        #
        for th in sp.where(self._throat_contested)[0]:
            #
            for p in self._net['throat.conns'][th]:
                g = self['throat.conductance'][th]
                fpc = self._sum_fpcap(th,p)
                self._RHS_pcap_data[p] +=  -g * fpc  #negative because cap fact was subtracted over to RHS

    def _calculate_dt(self):
        r"""
        Calculates the maximum timestep that would not advance a meniscus
        out of a throat or overfill a pore
        """
        #
        dt = sp.inf
        self._q = sp.zeros(self.Nt,dtype=sp.float64)
        self._pore_qsum = sp.zeros(self.Np,dtype=sp.float64)
        #
        # calculating q for contested throats
        for th in sp.where(self._throat_contested)[0]:
            if (self._net['throat.volume'][th] == 0.0):
                # if zero vol throats exist, dt must be 0.0 to maintain
                # proper mass balance, otherise injected fluid is 'lost'
                dt = 0.0
            #
            p1,p2 = self._net['throat.conns'][th]
            pr1,pr2 = self['pore.pressure'][[p1,p2]]
            g = self['throat.conductance'][th]
            fpc = self._sum_fpcap(th,p1)
            #
            # negative dir is moving away from lower index pore
            self._q[th] = -g * (pr1 - pr2 + fpc)
        #
        # setting dt values based on maximum allowed throat travel distance
        for th in sp.where(self._throat_contested)[0]:
            dx_max = self._set_dx_max(th)
            v = self._q[th]/self._net['throat.area'][th]
            if (v == 0.0):
                continue
            dt_new = dx_max * self._net['throat.length'][th]/abs(v)
            if (dt_new < dt):
                dt = dt_new
        #
        # estimating dt for either phase to reach dv_max
        for p in sp.where(self._pore_contested)[0]:
            if (self._net['pore.volume'][p] == 0.0):
                dt= 0.0
            con_ts = self._net.find_neighbor_throats(p)
            con_ps = self._net.find_connected_pores(con_ts)
            con_ps = con_ps[con_ps != p]
            con_ts_sf = self._get_supply_facts(con_ts,p)
            #
            # throats supplying pore
            qsum = 0.0
            for i in range(len(con_ts)):
                th = con_ts[i]
                p1,p2 = self._net['throat.conns'][th]
                pr1,pr2 = self['pore.pressure'][[p1,p2]]
                g = self['throat.conductance'][th]
                fpc = self._sum_fpcap(th,p1)
                q = -g * (pr1 - pr2 + fpc) # neg value is flowing out of pore 1
                self._q[th] = q
                if (p == p2):
                    q = -1.0 * q #reversing sign b/c we're looking at pore 2
                # only accounting for the invading phase
                if (con_ts_sf[i] > 0):
                    qsum += q
            #
            self._pore_qsum[p] = qsum
            dv_max = self._set_dv_max(p,qsum)
            if (qsum == 0.0):
                continue
            dt_new = dv_max * self._net['pore.volume'][p]/abs(qsum)
            if (dt_new < dt):
                dt = dt_new
        return(dt)

    def _advance_interface(self,dt):
        r"""
        Updates simulation based on the chosen time step
        """
        #
        contested_pores = sp.where(self._pore_contested)[0]
        #
        # moving mensici and setting up a new contested pore if necessary
        for th in sp.where(self._throat_contested)[0]:
            v = self._q[th]/self._net['throat.area'][th]
            # neg v is moving away from lowest index pore
            dx = (-v * dt)/self._net['throat.length'][th] #positive dx is moving away from lower index pore
            self._menisci[th] = [m + dx for m in self._menisci[th]]
            ph_frac = dx * self._throat_sup_fact[th]
            # if even number of mensici then phase is same on both ends of throat
            if (len(self._menisci[th]) % 2 == 0):
                ph_frac = 0.0
            self._throat_inv_frac[th] += ph_frac
            #
            # checking if throat has zero-volume
            if (self._net['throat.volume'][th] == 0.0):
                m = self._menisci[th][0]
                self._advance_zero_vol_throat(th)
                dx = self._menisci[th][0] - m
                v = -dx
            #
            mens = ['{:0.5f}'.format(m) for m in self._menisci[th]]
            fmt_str = 'Throat {:2d}: inv_frac: {:0.5f} menisci advanced by {:0.5f} new positions: {}'
            self._message(fmt_str.format(th,self._throat_inv_frac[th],dx,', '.join(mens)))
            sat_adj = 0.0
            pore = -1
            if ((self._menisci[th][-1] > (1.0 - self._sat_tol)) and (v < 0.0)): #mensicus being pushed away
                pore = self._net['throat.conns'][th][1] #pore meniscus moved into
                sat_adj = (1.0-self._menisci[th][-1])*self._throat_sup_fact[th]
                del self._menisci[th][-1] # menisci are ordered 0 -> 1
            elif ((self._menisci[th][0] < self._sat_tol) and (v > 0.0)): #meniscus being pulled inward
                pore = self._net['throat.conns'][th][0]
                sat_adj = self._menisci[th][0]*self._throat_sup_fact[th]
                del self._menisci[th][0]
                # needs flipped because fluid supplying throat changed as miniscus moves into the pore
                self._throat_sup_fact[th] *= -1.0
            #
            # updating saturations if rounding was performed
            if (pore > -1):
                # changing throat saturation based on rounding to pore
                self._throat_inv_frac[th] += sat_adj
                #negative b/c it's the fluid opposite the meniscus
                self._pore_inv_frac[pore] += -sat_adj
                self._pore_contested[pore] = True
                self._message('New contested pore: ',pore)
            #
            # removing contested flag if no mensici exist in throat
            if (len(self._menisci[th]) == 0):
                self._throat_contested[th] = False
        #
        # updating contested pores phase fraction
        for p in contested_pores:
            # qsum is always in terms of invading phase
            qsum = self._pore_qsum[p]
            if (self._net['pore.volume'][p] == 0.0):
                print('zero vol pore qsum: ',self._pore_qsum[p])
                if (qsum > 0):
                    self._pore_inv_frac[p] = 1.0
                else:
                    self._pore_inv_frac[p] = 0.0
            else:
                self._pore_inv_frac[p] += dt*qsum/self._net['pore.volume'][p]
            #
            #
            frac = dt*qsum
            fmt_str = 'Pore {0:2d} filled to: {1:10.6f}, ph frac change: {2:10.6f}, overall change: {3:10.9f}'
            self._message(fmt_str.format(p,self._pore_inv_frac[p],frac/self._net['pore.volume'][p],frac/self._net_vol))
            if (self._pore_inv_frac[p] > (1 - self._sat_tol)):
                if (qsum >= 0):
                    self._fill_pore(p)
            elif (self._pore_inv_frac[p] < self._sat_tol):
                if (qsum <= 0):
                    self._fill_pore(p)

    def _calc_fluid_out(self,dt):
        r"""
        Calculates the total amount of each phase leaving the network.
        """
        self._def_out_rate = 0.0
        self._inv_out_rate = 0.0
        self._def_out = 0.0
        self._inv_out = 0.0
        #
        for p in self._outlets:
            con_ts = self._net.find_neighbor_throats(p)
            con_ps = self._net.find_connected_pores(con_ts)
            con_ps = con_ps[con_ps != p]
            #
            # throats supplying pore
            for i in range(len(con_ts)):
                th = con_ts[i]
                p1,p2 = self._net['throat.conns'][th]
                pr1,pr2 = self['pore.pressure'][[p1,p2]]
                g = self['throat.conductance'][th]
                fpc = self._sum_fpcap(th,p1)
                q = -g * (pr1 - pr2 + fpc) # neg value is flowing out of pore 1
                self._q[th] = q
                if (p == p2):
                    q = -1.0 * q #reversing sign b/c we're looking at pore 2
                # only accounting for the invading phase
                if self._pore_invaded[p]:
                    self._inv_out_rate += q
                else:
                    self._def_out_rate += q
        self._inv_out = self._inv_out_rate * dt
        self._def_out = self._def_out_rate * dt
        self._total_inv_out += self._inv_out_rate * dt
        self._total_def_out += self._def_out_rate * dt

    def _print_step_stats(self,*args):
        #
        # getting average pressure drop (outlet is set to 0)
        inlet_p = sp.average(self['pore.pressure'][self._inlets])
        #
        q_inj = self._total_time * self._inj_rate
        vp = sp.multiply(self._net['pore.volume'],self._pore_inv_frac)
        vt = sp.multiply(self._net['throat.volume'],self._throat_inv_frac)
        tot_vol  = sp.sum(vp) + sp.sum(vt)
        tot_sat = tot_vol/self._net_vol
        mass_bal = (q_inj - tot_vol - self._total_inv_out)/self._net_vol
        fmt_str = 'Tot Sat Frac: {:0.5f}, Norm Mass Diff: {:0.15F}'
        #
        #
        #print(args[0],'  diff: {:15.9e}'.format((q_inj - tot_vol)/self._net_vol), ' dt: ',args[1])
        self._message('Net Def Fluid Out: {:10.6e}'.format(self._def_out))
        self._message('Net Inv Fluid Out: {:10.6e}'.format(self._inv_out))
        self._message('Net Fluid In: {:10.6e}'.format(self._inj_rate*args[1]))
        self._message('Net Fluid Out: {:10.6e}'.format(self._def_out+self._inv_out))
        self._message('Average Pressure Drop: {:10.4f}'.format(inlet_p))
        self._message(fmt_str.format(tot_sat,mass_bal))
        self._message('-'*25)
        self._message('')


#
# Helper functions below here
#
    def _sum_fpcap(self,th,ref_pore):
        r"""
        Sums the capillary forces from minisci alternating the sign with
        the fluid type.
        """
        # determining loop order
        p1,p2 = self._net['throat.conns'][th]
        step = 1
        if (ref_pore == p2):
            step = -1
        # initializing capillary factor
        fpc = 0.0
        f = -1.0*self._get_supply_facts([th],ref_pore)[0]# needs reversed b/c 1.0 is invading phase
        for x in self._menisci[th][::step]:
            fpc += f * self._pc_func(x)*self._max_pc[th]
            f = f * -1.0
        #
        return(fpc)

    def _get_supply_facts(self,throats,ref_pore):
        Ts_sf = sp.zeros(sp.size(throats))
        for i in range(sp.size(throats)):
            th = throats[i]
            p1,p2 = self._net['throat.conns'][th]
            if (ref_pore == p2):
                # sup facts are based on lower indexed pore, needs flipped
                # based on number of mensici preset for upper pore sf
                Ts_sf[i] = self._throat_sup_fact[th] * (-1)**len(self._menisci[th])
            else:
                Ts_sf[i] = self._throat_sup_fact[th]
        #
        return(Ts_sf)

    def _set_dx_max(self,th):
        q = self._q[th]
        if (q < 0.0):
            x = self._menisci[th][-1]
        else:
            x = self._menisci[th][0]
        dx_max = 0.03
        #
        if (q < 0.0) and (x > 0.50): #pushing meniscus away from pore past halfway point (0.51 ->1)
            dx_max = 0.30
            if (self._menisci[th][0] < 0.50):
                dx_max = 0.03
            if (dx_max > 1.0-x):
                dx_max = 1.0 - x
        elif (q > 0.0) and (x < 0.50): #pulling meniscus towards from pore past halfway point (0 -> 0.49)
            dx_max = 0.30
            if (self._menisci[th][-1] > 0.50):
                dx_max = 0.03
            if (dx_max > x):
                dx_max =  x
        #
        return(dx_max)

    def _set_dv_max(self,pore,q):
        #
        dv_max = 0.25
        #filling pore
        if ((q > 0) and ((1 - self._pore_inv_frac[pore]) < dv_max)):
            dv_max = 1 - self._pore_inv_frac[pore]
        #emptying pore
        elif ((q < 0) and (self._pore_inv_frac[pore] < dv_max)):
             dv_max = self._pore_inv_frac[pore]
        #
        return(dv_max)

    def _advance_zero_vol_throat(self,th):
        r"""
        Fills the throat with matching pore fluid based on the flow
        through it.
        """
        p1,p2 = self._net['throat.conns'][th];
        phase = 1.0
        # fluid flowing from p1 into p2
        if (self._q[th] < 0.0):
            if not (self._pore_invaded[p1]):
                phase = 0.0
            self._menisci[th] = [1.0]
        # fluid flowing from p2 into p1
        else:
            if not (self._pore_invaded[p2]):
                phase = 0.0
            self._menisci[th] = [0.0]
        self._throat_inv_frac[th] = phase

    def _fill_pore(self,pore):
        r"""
        Handles filling of pores and creation of new menisci in throats.
        """
        #
        self._pore_inv_frac[pore] = round(self._pore_inv_frac[pore])
        if (int(self._pore_inv_frac[pore]) == 1.0):
            self._pore_invaded[pore] = True
            sf = -1
        else:
            self._pore_invaded[pore] = False
            sf = 1
        #
        # creating a meniscus in all throats that have a supply factor matching
        # the pores previous status
        Ts = self._net.find_neighbor_throats(pore)
        Ts_sf = self._get_supply_facts(Ts,pore)
        self._set_menisci(pore,Ts[Ts_sf == sf])
        #
        # testing if all throats have the same sf, if so then contested is false
        Ts_sf = self._get_supply_facts(Ts,pore)
        self._pore_contested[pore] = not sp.all(Ts_sf == Ts_sf[0])

    def _set_menisci(self,base_pore,Ts):
        for th in Ts:
            p1,p2 = list(self._net['throat.conns'][th])
            if (base_pore == p2):
                # checking if flow in throat is moving away from the pore
                # - Q means flow from p1 into p2
                if (self._q[th] < 0.0):
                    continue
                self._menisci[th].append(1.0)
            else:
                # checking if flow in throat is moving away from the pore
                if (self._q[th] > 0.0):
                    continue
                self._menisci[th].insert(0,0.0)
                # needs flipped because fluid supplying throat changed
                self._throat_sup_fact[th] *= -1.0
            #
            self._throat_contested[th] = True


    def _message(self,*args):
        #
        string = [str(a) for a in args]
        string = ' '.join(string)
        print(string)
        self._log_file.write(string+'\n')



    def return_results(self,**kwargs):
        self._net['pore.inv_phase_frac'] = self._pore_inv_frac
        self._net['throat.inv_phase_frac'] = self._throat_inv_frac
        super().return_results(**kwargs)







