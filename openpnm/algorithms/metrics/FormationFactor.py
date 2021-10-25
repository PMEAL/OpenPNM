import numpy as np
from openpnm.utils import logging, Project, Workspace, PrintableDict
from openpnm.phases import GenericPhase
from openpnm.physics import GenericPhysics
from openpnm.algorithms import FickianDiffusion
from openpnm.algorithms.metrics import GenericMetric
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)
ws = Workspace()


class FormationFactor(GenericMetric):
    r"""
    This class works by applying 'value' boundary conditions across the
    domain to find molar flow, then using Fick's law to back-calculate
    the effective diffusivity of the domain.  The formation factor is
    defined as:

    .. math::

        F = \frac{D_{AB}}{D_{eff}} > 1

    where

    .. math::

        D_{eff} = \frac{n_{A} L }{\Delta C_{A} A }

    and

    .. math::

        D_{eff} = D_{AB} \frac{\varepsilon}{\tau}

    The formation factor is a convenient metric to compare diffusion in
    different pore networks since it does not require knowledge of the network
    porosity, unlike tortuosity. The porosity of a pore network is difficult
    to determine, mainly because the bulk volume of a network is not
    well known.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-5)
    >>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

    Now find the formation factor of the network:

    >>> F = op.algorithms.metrics.FormationFactor(network=pn)
    >>> F.run()
    >>> print(F.results)
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    Direction                           Formation Factor
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    x                                   ...
    y                                   ...
    z                                   ...
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.update({'inlets': {'x': 'left',
                                         'y': 'front',
                                         'z': 'top'},
                              'outlets': {'x': 'right',
                                          'y': 'back',
                                          'z': 'bottom'},
                              'areas': {'x': None,
                                        'y': None,
                                        'z': None},
                              'lengths': {'x': None,
                                          'y': None,
                                          'z': None}})

        self.results = PrintableDict()
        self.results._value = "Formation Factor"
        self.results._key = "Direction"
        if network is None:
            network = project.network
        if project is None:
            project = network.project
#        project = ws.copy_project(network.project)
#        keep = [network] + list(project.geometries().values())
#        for i in project:
#            if i not in keep:
#                project.purge_object(i)
        super().__init__(network=network, project=project, **kwargs)

    def run(self):
        r"""
        Execute the diffusion simulations in the principle directions.

        """
        phase = GenericPhase(network=self.network)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        mod = models.physics.diffusive_conductance.ordinary_diffusion
        for geom in self.project.geometries().values():
            phys = GenericPhysics(network=self.network,
                                  phase=phase, geometry=geom)
            phys.add_model(propname='throat.diffusive_conductance', model=mod)
        for bcs in self.settings['inlets'].keys():
            Diff = FickianDiffusion(network=self.project.network, phase=phase)
            Pin = self.network.pores(self.settings['inlets'][bcs])
            Diff.set_value_BC(pores=Pin, values=1.0)
            Pout = self.network.pores(self.settings['outlets'][bcs])
            Diff.set_value_BC(pores=Pout, values=0.0)
            Diff.run()
            A = self.settings['areas'][bcs]
            if A is None:
                A = Diff._get_domain_area()
                self.settings['areas'][bcs] = A
            L = self.settings['lengths'][bcs]
            if L is None:
                L = Diff._get_domain_length()
                self.settings['lengths'][bcs] = L
            R = Diff.rate(pores=Pin)
            Deff = R*L/A  # Conc gradient and diffusivity were both unity
            self.results[bcs] = 1/Deff[0]

    def set_inlets(self, direction, label):
        r"""
        Specifies which pore labels indicate the inlet face for a given
        direction

        Parameters
        ----------
        direction : string
            Either 'x', 'y', or 'z', indicating the principle direction
        label : string
            The pore label indicating which pores define the inlet face

        Notes
        -----
        By default inlets are taken as pores labelled 'left', 'top' and
        'front'. If your network does not have these labels you should
        either add them, or change which labels this algorithm looks for
        using this function, or manual changing the values in the
        ``settings['inlets']``.
        """
        self.settings['inlets'][direction] = label

    def set_outlets(self, direction, label):
        r"""
        Specifies which pore labels indicate the outlet face for a given
        direction

        Parameters
        ----------
        direction : string
            Either 'x', 'y', or 'z', indicating the principle direction
        label : string
            The pore label indicating which pores define the outlet face

        Notes
        -----
        By default outlets are taken as pores labelled 'right', 'bottom' and
        'back'. If your network does not have these labels you should
        either add them, or change which labels this algorithm looks for
        using this function, or manual changing the values in the
        ``settings['outlets']``.
        """
        self.settings['outlets'][direction] = label

    def set_area(self, direction, area):
        r"""
        Specifies the area of the inlet faces normal to the given direction

        Parameters
        ----------
        direction : string
            Either 'x', 'y', or 'z', indicating the principle direction
        area : scalar
            The numerical value of the inlet face area

        Notes
        -----
        If this value is not provided then an attempt is made to calculate it
        automatically, though this estimate is not perfectly accurate,
        especially for small networks.  This value is placed into the object's
        setting dictionary under 'areas'.
        """
        self.settings['areas'][direction] = area

    def set_length(self, direction, length):
        r"""
        Specifies the length of the domain along the given direction

        Parameters
        ----------
        direction : string
            Either 'x', 'y', or 'z', indicating the principle direction
        area : scalar
            The numerical value of the domain length

        Notes
        -----
        If this value is not provided then an attempt is made to calculate it
        automatically, though this estimate is not perfectly accurate,
        especially for small networks.  This value is placed into the object's
        setting dictionary under 'lengths'.
        """
        self.settings['lengths'][direction] = length
