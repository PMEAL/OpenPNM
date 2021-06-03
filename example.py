import openpnm as op
import numpy as np


def theta(target,
          alg_name,
          theta_max,
          m=3,
          diameter='pore.diameter',
          volume='pore.volume',
          area='pore.area'):
    r"""
    """
    from scipy.constants import Avogadro as N_A
    from scipy.constants import pi
    alg = target.project[alg_name]
    net = target.network
    try:
        t = alg.settings['t_solns'][-1]
        Cads = alg[f'pore.net_rate@{t}']
        r = net[diameter]
        Vp = net[volume]
        A = net[area]
        B = (1-pi*r**2*Vp/(A*theta_max) * Cads)**m
    except:
        B = np.zeros(alg.Np)
    return B


pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phases.Air(network=pn, name='air')
phys_air = op.physics.Standard(network=pn, phase=air, geometry=geo)

# Add reaction to phys_air
phys_air['pore.n'] = 1
phys_air['pore.k'] = -1e-5
phys_air.add_model(propname='pore.B',
                   model=theta,
                   alg_name=None,  # Will update later
                   theta_max=0.91,
                   m=3, regen_mode='deferred')
phys_air.add_model(propname='pore.A',
                   model=op.models.misc.product,
                   prop1='pore.k', prop2='pore.B',
                   regen_mode='deferred')
mod = op.models.physics.source_terms.standard_kinetics
phys_air.add_model(propname='pore.rxn', model=mod,
                   X='pore.concentration',
                   prefactor='pore.A', exponent='pore.n',
                   regen_mode='deferred')

rxn = op.algorithms.TransientFickianDiffusion(network=pn)
phys_air.models['pore.B']['alg_name'] = rxn.name  # This is kludgy for now
rxn.setup(phase=air)
Ps = pn.find_nearby_pores(pores=50, r=5e-4, flatten=True)
rxn.set_source(propname='pore.rxn', pores=Ps)
rxn.set_value_BC(pores=pn.pores('top'), values=1)
settings = {'store_rate': True,
            't_final': 100.0,
            't_initial': 0,
            't_output':10.0}
rxn.settings.update(settings)
rxn.run()
air.update(rxn.results())



