# -*- coding: utf-8 -*-
"""
===============================================================================
Boundary -- Subclass of GenericGeometry for Boundary Pores
===============================================================================

"""

from openpnm.models import geometry as gm
from openpnm.models import misc as mm
from openpnm.geometry import GenericGeometry


class Boundary(GenericGeometry):
    r"""
    Boundary subclass of GenericGeometry.

    The aim of this class is to assign geometric properties to boundary pores
    that are appropriate, such as 0 volume pores (for invasion simulations and
    porosity calculations) and very large throats (to offer negligible
    resistance in transport calculations).  When necessary, this class assumes
    that the pores are spherical and the throats are cylindrical.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network to which the Geometry object should be associated

    pores : array_like
        The pores indice where the Geometry should be applied

    throats : array_like
        The throat indices where the Geometry should be applied

    name : string, optional
        A name to help identify the object, if not given one is automatically
        generated

    project : OpenPNM Project, optional
        Either a Network or a Project must be specified

    Examples
    --------
    >>> import openpnm
    >>> pn = openpnm.network.Cubic(shape=[3, 3, 3])
    >>> Ps_int = pn.pores(labels=['top', 'bottom'], mode='not')
    >>> Ps_boun = pn.pores(labels=['top', 'bottom'], mode='union')
    >>> Ts_int = pn.throats(labels=['internal'])
    >>> Ts_boun = pn.throats(labels=['internal'], mode='not')
    >>> geo = openpnm.geometry.GenericGeometry(network=pn,
    ...                                        pores=Ps_int, throats=Ts_int)
    >>> boun = openpnm.geometry.Boundary(network=pn, pores=Ps_boun,
    ...                                  throats=Ts_boun)

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Temporary value to make function below work
        self.regenerate_models()

    def regenerate_models(self):
        self['pore.diameter'] = 0.0
        self.add_model(propname='throat.diameter',
                       model=mm.from_neighbor_pores,
                       pore_prop='pore.diameter',
                       mode='max')
        # Make slightly smaller than pore
        self['throat.diameter'] /= 1.1
        # Copy value from throat model just added
        self.add_model(propname='pore.diameter',
                       model=mm.from_neighbor_throats,
                       throat_prop='throat.diameter',
                       mode='max')
        # Make slightly bigger than throat
        self['pore.diameter'] *= 1.1
        del self.models['throat.diameter']
        del self.models['pore.diameter']
        self.add_model(propname='pore.area',
                       model=gm.pore_area.sphere)

        self['pore.volume'] = 0.0
        self['pore.seed'] = 1.0
        self['throat.seed'] = 1.0

        self.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        self.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        self.add_model(propname='throat.volume',
                       model=gm.throat_volume.extrusion)
        self.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)
        self.add_model(propname='throat.conduit_lengths',
                       model=gm.throat_length.boundary)
        self.add_model(propname='throat.equivalent_area',
                       model=gm.throat_equivalent_area.boundary)
