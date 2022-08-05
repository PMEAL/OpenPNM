import logging
import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator


__all__ = ['AdvectionDiffusion']


docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='AdvectionDiffusionSettings',
                     sections=['Parameters'])
@docstr.dedent
class AdvectionDiffusionSettings:
    r"""
    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s
    diffusive_conductance : str
        The name of the diffusive conductance values to be used by the
        specified ``conductance`` model to find the advective-diffusive
        conductance.
    hydraulic_conductance : str, optional
        The name of the hydraulic conductance values to be used by the
        specified ``conductance`` model to find the advective-diffusive
        conductance.
    pressure : str, optional
        The name of the pressure values calculated by the ``StokesFlow``
        algorithm.

    """
    quantity = 'pore.concentration'
    conductance = 'throat.ad_dif_conductance'
    diffusive_conductance = 'throat.diffusive_conductance'
    hydraulic_conductance = 'throat.hydraulic_conductance'
    pressure = 'pore.pressure'
    cache = False


class AdvectionDiffusion(ReactiveTransport):
    r"""
    A subclass of ReactiveTransport to simulate advection-diffusion
    """

    def __init__(self, name='ad_?', **kwargs):
        super().__init__(name=name, **kwargs)
        self.settings._update(AdvectionDiffusionSettings())
        self['pore.bc.outflow'] = np.nan

    def set_outflow_BC(self, pores, mode='add'):
        r"""
        Adds outflow boundary condition to the selected pores

        Parameters
        ----------
        pores : array_like
            The pore indices where the condition should be applied
        mode : str, optional
            Controls how the boundary conditions are applied. The default value
            is 'merge'. For definition of various modes, see the
            docstring for ``set_BC``.
        force : bool, optional
            If ``True`` then the ``'mode'`` is applied to all other bctypes as
            well. The default is ``False``.

        Notes
        -----
        Outflow condition means that the gradient of the solved quantity
        does not change, i.e. is 0.

        """
        pores = self._parse_indices(pores)

        # Calculating A[i,i] values to ensure the outflow condition
        network = self.project.network
        phase = self.project[self.settings['phase']]
        throats = network.find_neighbor_throats(pores=pores)
        C12 = network.conns[throats]
        P12 = phase[self.settings['pressure']][C12]
        gh = phase[self.settings['hydraulic_conductance']][throats]
        Q12 = -gh * np.diff(P12, axis=1).squeeze()
        Qp = np.zeros(self.Np)
        np.add.at(Qp, C12[:, 0], -Q12)
        np.add.at(Qp, C12[:, 1], Q12)

        self.set_BC(pores=pores, bcvalues=Qp[pores], bctype='outflow',
                    mode=mode)

    def _apply_BCs(self):
        r"""
        Applies Dirichlet, Neumann, and outflow BCs in order
        """
        # Apply Dirichlet and rate BCs
        super()._apply_BCs()
        if 'pore.bc.outflow' not in self.keys():
            return
        # Apply outflow BC
        diag = self.A.diagonal()
        ind = np.isfinite(self['pore.bc.outflow'])
        diag[ind] += self['pore.bc.outflow'][ind]
        self.A.setdiag(diag)


if __name__ == "__main__":
    import openpnm as op
    pn = op.network.Cubic(shape=[10, 10, 1])
    pn.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
    pn.regenerate_models()

    air = op.phase.Air(network=pn)
    air.add_model_collection(op.models.collections.physics.standard)
    air.regenerate_models()
    del air.models['throat.diffusive_conductance']
    del air.models['throat.hydraulic_conductance']
    air['throat.diffusive_conductance'] = 1e-15
    air['throat.hydraulic_conductance'] = 1e-15

    flow = op.algorithms.StokesFlow(network=pn, phase=air)
    flow.set_value_BC(pores=pn.pores('left'), values=1)
    flow.set_value_BC(pores=pn.pores('right'), values=0)
    flow.run()
    ad = op.algorithms.AdvectionDiffusion(network=pn, phase=air)
    ad.settings['cache'] = False
    ad.set_value_BC(pores=pn.pores('front'), values=1)
    ad.set_outflow_BC(pores=pn.pores('back'), mode='overwrite', force=False)
    ad.run()
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1)
    ax.imshow(ad['pore.concentration'].reshape([10, 10]))
