# -*- coding: utf-8 -*-
"""
===============================================================================
module __ViscousDrainage__: Viscous fluid flow with capillary pressure
===============================================================================

"""
import scipy as sp
import OpenPNM.Utilities.IO as io
from OpenPNM.Algorithms import GenericLinearTransport
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
              invading_phase,
              defending_phase=None,
              injection_rate=None,
              conductance='hydraulic_conductance',
              entry_pressure='throat.capillary_pressure',
              pore_volume='pore.volume',
              throat_volume='throat.volume',
              super_pore_conductance=None,
              sat_tol=1.0E-6,
              max_steps=1E5):
        r"""
        This setup provides the initial requirements for the solver setup
        and additional parameters for the drainage simulation.

        Parameters
        ----------
        invading_phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        defending_phase : OpenPNM Phase object
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
        max_steps : int (optional)
            Sets the overall hard limit on number of time steps to perform
            during the simulation if other exit criteria are not met.
        """
        logger.info('Setup ' + self.__class__.__name__)
        #
        if injection_rate is None:
            raise Exception('Error - injection rate must be specified')
        #
        if invading_phase is None:
            raise Exception('Error - Invading phase phase must be specified')
        #
        if defending_phase is None:
            defending_phase = self._phase
        #
        self['pore.inlets'] = False
        self['pore.outlets'] = False
        self['throat.entry_pressure'] = invading_phase[entry_pressure]
        self['pore.inv_frac'] = sp.zeros(self.Np, dtype=float)
        self['throat.inv_frac'] = sp.zeros(self.Nt, dtype=float)
        self['pore.contested'] = sp.zeros(self.Np, dtype=bool)
        self['throat.contested'] = sp.zeros(self.Nt, dtype=bool)
        self['pore.invaded'] = sp.zeros(self.Np, dtype=bool)
        self['pore.pressure'] = sp.zeros(self.Np)
        self._inv_phase = invading_phase
        self._def_phase = defending_phase
        self._inj_rate = injection_rate
        self._th_q = sp.zeros(self.Nt)
        self._pore_qsum = sp.zeros(self.Np)
        self._menisci = [[] for i in range(self.Nt)]
        self._throat_volume = throat_volume
        self._pore_volume = pore_volume
        # used in advance interface to specify how to change saturation fraction
        self._throat_sup_fact = sp.ones(self.Nt, dtype=float)*-1.0
        #
        self._sat_tol = sat_tol
        self._max_steps = max_steps
        self._ts_num = 0
        self._total_time = 0.0
        self.break_through_time = -1
        self.break_through_steps = 0
        self._net_vol = sp.sum(self._net[pore_volume])
        self._net_vol += sp.sum(self._net[throat_volume])
        self._total_def_out = 0.0
        self._total_inv_out = 0.0
        self._def_out_rate = 0.0
        self._inv_out_rate = 0.0
        #
        self._log_fname = 'VD-Log-2.txt'
        super().setup(conductance=conductance, quantity='pressure',
                      super_pore_conductance=super_pore_conductance)

    def set_inlets(self, pores=None, mode='add'):
        r"""
        Sets inlet pores as well as inital menicsi to start simulation.
        """
        #
        Ps = self._parse_locations(pores)
        #
        if mode in ['clear', 'overwrite']:
            self['pore.inlets'] = False
        #
        if sum(self['pore.outlets'][Ps]) > 0:
            raise Exception('Some inlets are already defined as outlets')
        #
        bool_val = True
        if mode is 'remove':
            bool_val = False
        #
        self['pore.inlets'][Ps] = bool_val
        inlets = sp.where(self['pore.inlets'])[0]
        #
        # setting invasion status and BC on inlet pores
        self['pore.inv_frac'][inlets] = 1.0
        self['pore.invaded'][inlets] = True
        self.set_boundary_conditions(bctype='Neumann_group',
                                     mode='merge',
                                     bcvalue=-self._inj_rate,
                                     pores=inlets)
        logger.debug('Inlet pores set as invaded and Nuemann BC defined')
        #
        # throats between two inlets are set as filled to prevent plugs
        filled_throats = self._net.find_neighbor_throats(pores=inlets,
                                                         mode='intersection')
        if filled_throats:
            self['throat.inv_frac'][filled_throats] = 1.0
            logger.debug('Throats between inlet pores have been set as filled')
        #
        # adding menisci
        interface = self._net.find_neighbor_throats(pores=inlets,
                                                    mode='not_intersection',
                                                    flatten=False)
        self['throat.contested'][sp.ravel(interface)] = True
        self._set_menisci(inlets, interface)
        logger.info('Set menisci for throats connected to inlet pores')

    def set_outlets(self, pores=None, mode='add'):
        r"""
        Defines outlets for invading and defending phases
        """
        #
        Ps = self._parse_locations(pores)
        if mode in ['clear', 'overwrite']:
            self['pore.outlets'] = False
        #
        if sum(self['pore.inlets'][Ps]) > 0:
            raise Exception('Some outlets are already defined as inlets')
        #
        bool_val = True
        if mode is 'remove':
            bool_val = False
        self['pore.outlets'][Ps] = bool_val
        #
        self._outlets = sp.where(self['pore.outlets'])[0]
        self.set_boundary_conditions(bctype='Dirichlet',
                                     mode='overwrite',
                                     bcvalue=0.0,
                                     pores=self._outlets)

    def run(self, **kwargs):
        r"""
        Starts the simulation
        """
        #
        # Ensure inlets are set
        if sp.size(self['pore.inlets']) == 0:
            raise Exception('Inlet pores have not been specified')

        # Ensure outlet pores are set
        if sp.size(self['pore.outlets']) == 0:
            raise Exception('Outlet pores have not been specified')
        #
        # gdef is reused to calculate mixed throat conductances
        self._gdef = sp.copy(self['throat.conductance'])
        #
        # setting initial time to maintain mass balance if inlet pores and
        # throats have a non-zero volume
        tot_vol = sp.sum(sp.multiply(self._net[self._pore_volume],
                                     self['pore.inv_frac']))
        tot_vol += sp.sum(sp.multiply(self._net[self._throat_volume],
                                      self['throat.inv_frac']))
        tot_sat = tot_vol/self._net_vol
        self._total_time = tot_vol/self._inj_rate
        logger.info('Initial Saturation of Invading Phase: '+str(tot_sat))
        #
        # beginning simulation
        with open(self._log_fname, 'w') as self._log_file:
            self._do_outer_iteration_stage(**kwargs)

    def restart_simulation(self, max_steps=2e5, **kwargs):
        r"""
        Restarts a simulation to run until an exit condition is met
        """
        #
        self._max_steps = max_steps
        logger.debug('Simulation restarted')
        #
        raise NotImplementedError

    def _do_outer_iteration_stage(self, **kwargs):
        r"""
        This calls the solve method in the algorithm.
        Handles the tracking and movement of phases throat the network.
        """
        #
        self._zero_dt = 0
        while True:
            A = self._update_coefficient_matrix()
            b = self._update_rhs()
            self.solve(A, b)
            dt = self._calculate_dt()
            #
            if self._ts_num % 10000 == 0:
                suf = str(self._ts_num)[0]
                self.return_results()
                phases = [self._inv_phase, self._def_phase]
                fname = 'temp_files/'+self._net.name+'-vd-temp-file-'+suf
                io.VTK.save(self._net, fname, phases)
            elif self._ts_num % 2000 == 0:
                self.return_results()
                phases = [self._inv_phase, self._def_phase]
                fname = 'temp_files/'+self._net.name+'-vd-temp-file'
                io.VTK.save(self._net, fname, phases)
            #
            if dt == 0.0:
                self._zero_dt += 1
            #
            test = sp.where(self['pore.inv_frac'][self._outlets] > 1-self._sat_tol)[0]
            if sp.size(test) > 0 and self.break_through_time < 0:
                self.break_through_time = self._total_time
                self.break_through_steps = self._ts_num
                #break
            #
            if (abs(self._inv_out_rate - self._inj_rate)/self._inj_rate < 1e-9):
                break
            #
            if self._ts_num > self._max_steps:
                logger.info('Maximum step exit condition triggered')
                break
            #
            #
            self._message('Time Step: ', self._ts_num, ' size: {:0.3E} '.format(dt))
            #
            self._advance_interface(dt)
            self._total_time += dt
            self._calc_fluid_out(dt)
            self._print_step_stats(self._ts_num, dt)
            self._ts_num += 1
        #
        # checking overall mass balance
        q_inj = self._total_time * self._inj_rate
        tot_vol = sp.sum(sp.multiply(self._net[self._pore_volume], self['pore.inv_frac']))
        tot_vol += sp.sum(sp.multiply(self._net[self._throat_volume], self['throat.inv_frac']))
        tot_sat = tot_vol/self._net_vol
        mass_bal = (q_inj - tot_vol - self._total_inv_out)/self._net_vol
        #
        self._message('Total Simulation Time Until Break Through: ',
                      self.break_through_time, ' Steps:', self.break_through_steps)
        self._message('Total Simulation Time: ', self._total_time, ' Steps:', self._ts_num)
        self._message('Total Volume: ', tot_vol)
        self._message('Total Inv Fluid Out: ', self._total_inv_out)
        self._message('Total saturation: ', tot_sat)
        self._message('Total injection: ', q_inj)
        self._message('Mass Difference / Total Vol: {:15.9e}'.format(mass_bal))

    def _update_coefficient_matrix(self):
        r"""
        Updates the conductance based on the viscosity ratio of the fluids
        and the fractional occupancy of the throats.
        """
        #
        conns = self._net['throat.conns']
        dvisc = sp.average(self._def_phase['pore.viscosity'][conns],axis=1)
        ivisc = sp.average(self._inv_phase['pore.viscosity'][conns],axis=1)
        M = sp.divide(dvisc,ivisc)
        #
        fact = 1.0 - self['throat.inv_frac']
        fact += sp.multiply(self['throat.inv_frac'],M)
        self['throat.conductance'] = sp.multiply(fact,self._gdef)
        #
        return self._build_coefficient_matrix()

    def _update_rhs(self):
        r"""
        Adds f * g * pcap to RHS for pores containing menisci
        """
        #
        Ts = sp.where(self['throat.contested'])[0]
        Ps = sp.ravel(self._net['throat.conns'][Ts], order='F')
        Ts = sp.append(Ts,Ts)
        g = self['throat.conductance'][Ts]
        fpc = self._sum_fpcap(Ts,Ps)
        fpc = sp.multiply(g,fpc)
        #
        rhs_pcap_data = sp.zeros(self.Np, dtype=float)
        for p,f in zip(Ps,fpc):
            rhs_pcap_data[p] -= f
        #
        # just tried adding this term for kicks, I don't believe it actually should exist
        #rhs_pcap_data -= self.rate(self._net.pores(),mode='single')
        #
        return self._build_RHS_matrix(self._net.pores(), rhs_pcap_data)

    def _calculate_dt(self):
        r"""
        Calculates the maximum timestep that would not advance a meniscus
        out of a throat or overfill a pore
        """
        #
        # initial time step is the time to fill 1% of the total network volume
        dt = self._net_vol/self._inj_rate*0.01
        self._th_q = sp.zeros(self.Nt)
        self._pore_qsum = sp.zeros(self.Np)
        #
        # calculating flow for contested throats
        Ts = sp.where(self['throat.contested'])[0]
        if sp.size(sp.where(self._net[self._throat_volume][Ts] == 0.0)[0]):
            dt = 0.0
        else:
            Ps = self._net['throat.conns'][Ts]
            Ps_pr = self['pore.pressure'][Ps]
            g = self['throat.conductance'][Ts]
            fpc = self._sum_fpcap(Ts, Ps[:,0])
            # negative q is flowing away from lower index pore
            self._th_q[Ts] = -sp.multiply(g,(Ps_pr[:,0] - Ps_pr[:,1] - fpc))
        #
        # setting dt based on maximum allowed meniscus travel distance
        dx_max = self._set_dx_max(Ts)
        v = sp.divide(self._th_q[Ts],self._net['throat.area'][Ts])
        v = sp.absolute(v)
        locs = sp.where(v > 0.0)[0]
        if sp.size(locs):
            Ts = Ts[locs]
            dx_max = dx_max[locs]
            v = v[locs]
            dt_new = sp.multiply(dx_max,(self._net['throat.length'][Ts]/v))
            dt_new = sp.amin(dt_new)
            if dt_new < dt:
                dt = dt_new
        #
        # calculating flow for contested pores
        Ps = self.pores('contested')
        if sp.size(sp.where(self._net[self._pore_volume][Ps] == 0.0)[0]):
            dt = 0.0
        else:
            self._pore_qsum[Ps] = self.rate(Ps,mode='single',phase='invading')
        #
        # setting dt based on maximum allowed pore volume change
        dv_max = self._set_dv_max(Ps, self._pore_qsum[Ps])
        qsum = sp.absolute(self._pore_qsum[Ps])
        locs = sp.where(qsum > 0.0)[0]
        if sp.size(locs):
            Ps = Ps[locs]
            dv_max = dv_max[locs]
            qsum = qsum[locs]
            dt_new = sp.multiply(dv_max,(self._net[self._pore_volume][Ps]/qsum))
            dt_new = sp.amin(dt_new)
            if dt_new < dt:
                dt = dt_new
        #
        return dt

    def _advance_interface(self, dt):
        r"""
        Updates simulation based on the chosen time step
        """
        #
        contested_pores = sp.where(self['pore.contested'])[0]
        #
        # moving mensici and setting up a new contested pore if necessary
        for th in sp.where(self['throat.contested'])[0]:
            # neg v is moving away from lowest index pore
            v = self._th_q[th]/self._net['throat.area'][th]
            # positive dx is moving away from lower index pore
            dx = (-v * dt)/self._net['throat.length'][th]
            self._menisci[th] = [m + dx for m in self._menisci[th]]
            ph_frac = dx * self._throat_sup_fact[th]
            # if even number of mensici then phase is same on both ends of throat
            if (len(self._menisci[th]) % 2 == 0):
                ph_frac = 0.0
            self['throat.inv_frac'][th] += ph_frac
            #
            # checking if throat has zero-volume
            if self._net[self._throat_volume][th] == 0.0:
                m = self._menisci[th][0]
                self._advance_zero_vol_throat(th)
                dx = self._menisci[th][0] - m
                v = -dx
            #
            #mens = ['{:0.5f}'.format(m) for m in self._menisci[th]]
            #fmt_str = 'Throat {:2d}: inv_frac: {:0.5f} menisci '
            #fmt_str += 'advanced by {:0.5f} new positions: {}'
            #self._message(fmt_str.format(th, self['throat.inv_frac'][th], dx, ', '.join(mens)))
            sat_adj = 0.0
            pore = -1
             #mensicus being pushed away from p1
            if ((self._menisci[th][-1] > (1.0 - self._sat_tol)) and (v < 0.0)):
                pore = self._net['throat.conns'][th][1] #pore meniscus moved into
                sat_adj = (1.0-self._menisci[th][-1])*self._throat_sup_fact[th]
                del self._menisci[th][-1] # menisci are ordered 0 -> 1
            #meniscus being pulled towards p1
            elif ((self._menisci[th][0] < self._sat_tol) and (v > 0.0)):
                pore = self._net['throat.conns'][th][0]
                sat_adj = self._menisci[th][0]*self._throat_sup_fact[th]
                del self._menisci[th][0]
                # needs flipped because fluid supplying throat changed as miniscus moves into p1
                self._throat_sup_fact[th] *= -1.0
            #
            # updating saturations if rounding was performed
            if pore > -1:
                # changing throat saturation based on rounding to pore
                self['throat.inv_frac'][th] += sat_adj
                #negative b/c it's the fluid opposite the meniscus
                self['pore.inv_frac'][pore] += -sat_adj
                self['pore.contested'][pore] = True
                #self._message('New contested pore: ', pore)
            #
            # removing contested flag if no mensici exist in throat
            if len(self._menisci[th]) == 0:
                self['throat.contested'][th] = False
        #
        # updating contested pores phase fraction
        for p in contested_pores:
            # qsum is always in terms of invading phase
            qsum = self._pore_qsum[p]
            if self._net[self._pore_volume][p] == 0.0:
                if qsum > 0:
                    self['pore.inv_frac'][p] = 1.0
                else:
                    self['pore.inv_frac'][p] = 0.0
            else:
                self['pore.inv_frac'][p] += dt*qsum/self._net[self._pore_volume][p]
            #
            #
            #frac = dt*qsum
            #fmt_str = 'Pore {0:2d} filled to: {1:10.6f}, ph frac change: '
            #fmt_str +='{2:10.6f}, overall change: {3:10.9f}'
            #self._message(fmt_str.format(p, self['pore.inv_frac'][p],
            #                             frac/self._net[self._pore_volume][p],
            #                             frac/self._net_vol))
            if (self['pore.inv_frac'][p] > (1 - self._sat_tol)):
                if qsum >= 0:
                    self._fill_pore(p)
            elif self['pore.inv_frac'][p] < self._sat_tol:
                if qsum <= 0:
                    self._fill_pore(p)

    def _calc_fluid_out(self, dt):
        r"""
        Calculates the total amount of each phase leaving the network.
        """
        def_out_rate = self.rate(pores=self._outlets, phase='defending')[0]
        self._inv_out_rate = self.rate(pores=self._outlets, phase='invading')[0]
        self._inv_out = self._inv_out_rate * dt
        self._def_out = def_out_rate * dt
        self._total_inv_out += self._inv_out_rate * dt
        self._total_def_out += def_out_rate * dt

    def _print_step_stats(self, *args):
        #
        # getting average pressure drop (outlet is set to 0)
        inlets = sp.where(self['pore.inlets'])[0]
        inlet_p = sp.average(self['pore.pressure'][inlets])
        #
        q_inj = self._total_time * self._inj_rate
        pore_vol = sp.multiply(self._net[self._pore_volume], self['pore.inv_frac'])
        throat_vol = sp.multiply(self._net[self._throat_volume], self['throat.inv_frac'])
        tot_vol = sp.sum(pore_vol) + sp.sum(throat_vol)
        tot_sat = tot_vol/self._net_vol
        mass_bal = (q_inj - tot_vol - self._total_inv_out)/self._net_vol
        fmt_str = 'Tot Sat Frac: {:0.5f}, Norm Mass Diff: {:0.15F}'
        #
        chk_val = abs(self._inv_out_rate - self._inj_rate)/self._inj_rate
        strg = 'inv fluid out: {:15.6e}, normed value: {:15.6e}'.format(self._inv_out_rate,
                                                                        chk_val)
        #
        #print(args[0], strg, 'num zero steps: ', self._zero_dt)
        self._message(strg, 'num zero steps: ', self._zero_dt)
        self._message('Net Def Fluid Out: {:10.6e}'.format(self._def_out))
        self._message('Net Inv Fluid Out: {:10.6e}'.format(self._inv_out))
        self._message('Net Fluid In: {:10.6e}'.format(self._inj_rate*args[1]))
        self._message('Net Fluid Out: {:10.6e}'.format(self._def_out+self._inv_out))
        self._message('Average Pressure Drop: {:10.4f}'.format(inlet_p))
        self._message(fmt_str.format(tot_sat, mass_bal))
        self._message('-'*25)
        self._message('')

#
# Helper functions below here
#
    def _sum_fpcap(self, throats, ref_pores):
        r"""
        Sums the capillary forces from minisci alternating the sign with
        the fluid type.
        """
        #
        fpc = sp.zeros(sp.size(throats))
        for i,(th,ref_pore) in enumerate(zip(throats,ref_pores)):
            # determining loop order
            ps = self._net['throat.conns'][th]
            step = 1
            if ref_pore == ps[1]:
                step = -1
            # needs reversed b/c 1.0 is invading phase
            f = self._get_supply_facts([[th]], [ref_pore])[0][0]
            for x in self._menisci[th][::step]:
                fpc[i] += f * abs(sp.sin(sp.pi*x)) * self['throat.entry_pressure'][th]
                f = f * -1.0
        #
        return fpc

    def _get_supply_facts(self, throat_list, ref_pores):
        r"""
        Handles calculation of supply factors
        """
        sup_facts = []
        for throats,ref_pore in zip(throat_list,ref_pores):
            #
            Ts_sf = sp.zeros(sp.size(throats))
            for i,th in enumerate(throats):
                ps = self._net['throat.conns'][th]
                if ref_pore == ps[1]:
                    # sup facts are based on lower indexed pore, needs flipped
                    # based on number of mensici preset for upper pore sf
                    Ts_sf[i] = self._throat_sup_fact[th] * (-1)**len(self._menisci[th])
                else:
                    Ts_sf[i] = self._throat_sup_fact[th]
            sup_facts.append(Ts_sf)
        #
        return sp.array(sup_facts,ndmin=2)

    def _set_dx_max(self, throats):
        r"""
        Calculates maximum mensicus travel distance for an array of throats
        """
        #
        dists = []
        for th in throats:
            q = self._th_q[th]
            if q < 0.0:
                x = self._menisci[th][-1]
            else:
                x = self._menisci[th][0]
            dx_max = 0.03
            #
            #pushing meniscus away from pore past halfway point (0.51 ->1)
            if (q < 0.0) and (x > 0.50):
                dx_max = 0.30
                if self._menisci[th][0] < 0.50:
                    dx_max = 0.03
                if (dx_max > 1.0 - x):
                    dx_max = 1.0 - x
            #pulling meniscus towards from pore past halfway point (0 -> 0.49)
            elif (q > 0.0) and (x < 0.50):
                dx_max = 0.30
                if self._menisci[th][-1] > 0.50:
                    dx_max = 0.03
                if dx_max > x:
                    dx_max = x
            dists.append(dx_max)
        #
        return sp.array(dists)

    def _set_dv_max(self, pores, qs):
        #
        vols = []
        for pore,q in zip(pores,qs):
            dv_max = 0.25
            #filling pore
            if ((q > 0) and ((1 - self['pore.inv_frac'][pore]) < dv_max)):
                dv_max = 1 - self['pore.inv_frac'][pore]
            #emptying pore
            elif ((q < 0) and (self['pore.inv_frac'][pore] < dv_max)):
                dv_max = self['pore.inv_frac'][pore]
            #
            vols.append(dv_max)
        #
        return sp.array(vols)

    def _advance_zero_vol_throat(self, th):
        r"""
        Fills the throat with matching pore fluid based on the flow
        through it.
        """
        p1, p2 = self._net['throat.conns'][th]
        phase = 1.0
        # fluid flowing from p1 into p2
        if self._th_q[th] < 0.0:
            if not self['pore.invaded'][p1]:
                phase = 0.0
            self._menisci[th] = [1.0]
        # fluid flowing from p2 into p1
        else:
            if not self['pore.invaded'][p2]:
                phase = 0.0
            self._menisci[th] = [0.0]
        self['throat.inv_frac'][th] = phase

    def _fill_pore(self, pore):
        r"""
        Handles filling of pores and creation of new menisci in throats.
        """
        #
        self['pore.inv_frac'][pore] = round(self['pore.inv_frac'][pore])
        if int(self['pore.inv_frac'][pore]) == 1.0:
            self['pore.invaded'][pore] = True
            sf = -1
        else:
            self['pore.invaded'][pore] = False
            sf = 1
        #
        # creating a meniscus in all throats that have a supply factor matching
        # the pores previous status
        Ts = self._net.find_neighbor_throats(pore)
        Ts_sf = self._get_supply_facts([Ts], [pore])[0]
        self._set_menisci([pore], [Ts[Ts_sf == sf]])
        #
        # testing if all throats have the same sf, if so then contested is false
        Ts_sf = self._get_supply_facts([Ts], [pore])[0]
        self['pore.contested'][pore] = not sp.all(Ts_sf == Ts_sf[0])

    def _set_menisci(self, pores, throats):
        r"""
        Sets the menisci for a list of pores and corresponding throats
        """
        for i, base_pore in enumerate(pores):
            Ts = throats[i]
            for th in Ts:
                ps = list(self._net['throat.conns'][th])
                if base_pore == ps[1]:
                    # checking if flow in throat is moving away from the pore
                    # - Q means flow from p1 into p2
                    if self._th_q[th] < 0.0:
                        continue
                    self._menisci[th].append(1.0)
                else:
                    # checking if flow in throat is moving away from the pore
                    if self._th_q[th] > 0.0:
                        continue
                    self._menisci[th].insert(0, 0.0)
                    # needs flipped because fluid supplying throat changed
                    self._throat_sup_fact[th] *= -1.0
                #
                self['throat.contested'][th] = True


    def _message(self, *args):
        #
        string = [str(a) for a in args]
        string = ' '.join(string)
        #print(string)
        self._log_file.write(string+'\n')

    def rate(self, pores=None, mode='group', phase='both'):
        r"""
        Send a list of pores and receive the net rate
        of material moving into them.

        Parameters
        ----------
        pores : array_like
            The pores where the net rate will be calculated
        mode : string, optional
            Controls how to return the rate.  Options are:
            - 'group'(default): It returns the cumulative rate moving into them
            - 'single': It calculates the rate for each pore individually.
        phase : string, optional
            Specifies a specific phase to calculate rates for.  Options are:
            - 'both'(default): It returns the cumulative rate of both phases
            - 'defending': Calculates the cumlative rate of only defending phase.
            - 'invading': Calculates the cumlative rate of only invading phase.
        """
        network = self._net
        conductance = self['throat.conductance']
        X_value = self[self._quantity]
        pores = sp.array(pores, ndmin=1)
        R = []
        if mode == 'group':
            t = network.find_neighbor_throats(pores, flatten=True,
                                              mode='not_intersection')
            throat_group_num = 1
        elif mode == 'single':
            t = network.find_neighbor_throats(pores, flatten=False,
                                              mode='not_intersection')
            #fixes odd scipy error where I get 2-D array instead of an array of arrays
            throat_group_num = sp.size(t,0)
        for i in sp.r_[0: throat_group_num]:
            if mode == 'group':
                throats = t
                P = pores
            elif mode == 'single':
                throats = t[i]
                P = pores[i]
            p1 = network.find_connected_pores(throats)[:, 0]
            p2 = network.find_connected_pores(throats)[:, 1]
            pores1 = sp.copy(p1)
            pores2 = sp.copy(p2)
            #
            # Changes to pores1 and pores2 to make them as inner/outer pores
            pores1[~sp.in1d(p1, P)] = p2[~sp.in1d(p1, P)]
            pores2[~sp.in1d(p1, P)] = p1[~sp.in1d(p1, P)]
            X1 = X_value[pores1]
            X2 = X_value[pores2]
            g = conductance[throats]
            fpc = self._sum_fpcap(throats,pores1)
            throat_q = -sp.multiply(g, (X1 - X2 - fpc))
            #
            # adding flow into the throat q array
            # needs flipped here b/c throats always assume flow direction
            # from lower index to higher
            flip_q = sp.where(pores1 > pores2)[0]
            throat_q[flip_q] = sp.negative(throat_q[flip_q])
            self._th_q[throats] = throat_q
            # reflipping q to reference flow into pores
            throat_q[flip_q] = sp.negative(throat_q[flip_q])
            #
            # determining if a speific phase was requested
            if phase == 'invading':
                Ts_sf = self._get_supply_facts(sp.array(throats,ndmin=2).T,pores1)
                Ts_sf = sp.ravel(Ts_sf)
                throat_q = throat_q[sp.where(Ts_sf > 0)[0]]
            elif (phase == 'defending'):
                Ts_sf = self._get_supply_facts(sp.array(throats,ndmin=2).T,pores1)
                Ts_sf = sp.ravel(Ts_sf)
                throat_q = throat_q[sp.where(Ts_sf < 0)[0]]
            else:
                raise Exception('Unsupported phase use "invading" or "defending"')
            #
            # summing throats for overall pore rate
            R.append(sp.sum(throat_q))
        return(sp.array(R, ndmin=1))

    def return_results(self, **kwargs):
        #
        for element in ['pore', 'throat']:
            prop_name = element+'.volume_fraction'
            self._inv_phase[prop_name] = self[element+'.inv_frac']
            self._def_phase[prop_name] = 1.0 - self[element+'.inv_frac']
        #
        self._net['pore.contested'] = self['pore.contested']
        self._net['pore.invaded'] = self['pore.invaded']
        self._net['throat.menisci'] = [len(men) for men in self._menisci]
        super().return_results(**kwargs)







