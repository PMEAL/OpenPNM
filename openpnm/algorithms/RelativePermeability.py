from openpnm.algorithms import GenericAlgorithm, StokesFlow
from openpnm.utils import logging
from openpnm import models
import numpy as np
import openpnm as op
logger = logging.getLogger(__name__)


default_settings = {'pore_inv_seq': 'pore.invasion_sequence',
                    'throat_inv_seq': 'throat.invasion_sequence',
                    'points': 10,
                    'gh': 'throat.hydraulic_conductance',
                    'mode': 'strict',
                    }

class RelativePermeability(GenericAlgorithm):
    r"""
    Parameters
    ----------


    Notes
    -----


    """
    def __init__(self, settings={}, **kwargs):
        # Apply default settings
        self.settings.update(default_settings)
        # Apply any settings received during initialization
        self.settings.update(settings)
        super().__init__(**kwargs)
    def setup(self, inv_phase=None, def_phase=None,points=None,
              pore_inv_seq=None,
              throat_inv_seq=None):
        r"""
        """
        if inv_phase:
            self.settings['inv_phase'] = inv_phase
        if def_phase:
            self.settings['def_phase'] = def_phase
        if points:
            self.settings['points'] = points
        if pore_inv_seq:
            self.settings['pore_inv_seq'] = pore_inv_seq
        else:
                res=self.IP(self)
                self.settings['pore_inv_seq'] =res[0]
        if throat_inv_seq:
            self.settings['thorat_inv_seq'] = throat_inv_seq
    def IP(self,):
        inv=op.algorithms.InvasionPercolation(phase=self.settings['inv_phase']),network=
        self.project.network,project=self.project)
        inv.setup(phase=oil,entry_pressure='throat.entry_pressure',pore_volume='pore.volume', throat_volume='throat.volume')
        inlets = pn.pores(['top'])
        used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
        outlets = pn.pores(['bottom'])
        used_outlets = [outlets[x] for x in range(0, len(outlets), 2)]
        inv.set_inlets(pores=used_inlets)
        inv.run()
        Snwparr =  []
        Pcarr =  []
        Sarr=np.linspace(0,1,num=60)
        for Snw in Sarr:
            res1=inv.results(Snwp=Snw)
            occ_ts=res1['throat.occupancy']
            if np.any(occ_ts):
                max_pthroat=np.max(phys_oil['throat.entry_pressure'][occ_ts])
                Pcarr.append(max_pthroat)
                Snwparr.append(Snw)
        self.settings['pore_inv_seq'] = pore_inv_seq
        self.settings['thorat_inv_seq'] = throat_inv_seq
        

    def set_inlets(self, pores):
        r"""
        """
        self['pore.inlets'] = False
        self['pore.inlets'][pores] = True

    def set_outlets(self, pores):
        r"""
        """
        self['pore.outlets'] = False
        self['pore.outlets'][pores] = True

    def run(self, inlets=None, outlets=None):
        r"""
        """
        if inlets is not None:
            self.set_inlets(pores=inlets)
        if outlets is not None:
            self.set_outlets(pores=outlets)
        # Retrieve phase and network
        network = self.project.network
        phase = self.project.phases(self.settings['phase'])
        sf = StokesFlow(network=network)
        sf.setup(phase=phase,
                 quantity='pore.pressure',
                 conductance=self.settings['gh'])
        sf.set_value_BC(pores=self['pore.inlets'], values=1)
        sf.set_value_BC(pores=self['pore.outlets'], values=0)
        phase['pore.occupancy'] = 0
        phase['throat.occupancy'] = 0
        phase.add_model(propname='throat.multiphase_hydraulic_conductance',
                        model=models.physics.multiphase.conduit_conductance,
                        throat_conductance=self.settings['gh'],
                        throat_occupancy='throat.occupancy',
                        pore_occupancy='pore.occupancy',
                        mode=self.settings['mode'],
                        )
        max_inv = np.amax(phase['pore.invasion_sequence'])
        invs = np.linspace(0, max_inv, self.settings['points'])
        results = []
        for s in invs:
            phase['pore.occupancy'] = 1.0*(phase[self.settings['pore_inv_seq']] <= s)
            phase['throat.occupancy'] = 1.0*(phase[self.settings['throat_inv_seq']] <= s)
            phase.regenerate_models(deep=True)
            sf.run()
            results.append([[s, sf.rate(pores=self['pore.inlets'])]])
        return results
