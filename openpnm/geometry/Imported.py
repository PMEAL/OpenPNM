import openpnm.models as mods
from openpnm.geometry import GenericGeometry
from openpnm.utils import logging, GenericSettings
logger = logging.getLogger(__name__)


class ImportedSettings(GenericSettings):
    r"""

    Parameters
    ----------

    pore_diameter : str (default = 'pore.extended_diameter')
        Key into the extracted data array to use as pore diameter in other
        geometry calculations. The default is .  Use of 'pore.' is not
        required.
    throat_diameter : str (default = 'throat.equivalent_diameter')
        Key into the extracted data array to use as throat diameter in other
        geometry calculations. Use of 'throat.' is not required.
    pore_shape : string {'sphere' (default), 'cube'}
        Specifies which shape to assume when calculating dependent properties
        such as volume and surface area.
    throat_shape : string {'cylinder' (default), 'cuboid'}
        Specifies which shape to assume when calculating dependent properties
        such as volume and surface area.

    """
    pore_diameter = 'pore.extended_diameter'
    throat_diameter = 'throat.equivalent_diameter'
    pore_shape = 'sphere'
    throat_shape = 'cylinder'


class Imported(GenericGeometry):
    r"""
    This geometry class extracts all numerical properites from the received
    network object and moves them to itself.

    This class is intended for use with networks imported from network
    extraction codes, where the geometry properties are included on the
    network itself.

    Parameters
    ----------
    network : OpenPNM Network object
        The network with which this Geometry should be associated

    exclude : list of strings
        A list of which network properties should *not* be transferred to
        new geometry object.  'pore.coords' and 'throat.conns' are *always*
        excluded.  Note that labels are not transferred, only properties.

    project : OpenPNM Project object, optional
        Can be supplied in addition to ``network`` but is inferred from the
        network's project if not given.

    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    Notes
    -----
    An error occurs when adding other geometries to a network that has
    geometrical properties such as 'pore.diameter'.  This can occur when
    adding boundary pores or in more elaborate scenarios such as stitching
    networks together.  The issue arises because OpenPNM prevents a property,
    such as 'pore.volume', from existing on both the network and also a
    geometry.  Thus it is necessary to move the extracted network properties
    to this ``Imported`` class, then create new geometry objects for any
    added pores as needed.

    """

    def __init__(self, network, exclude=[], settings={}, **kwargs):
        self.settings._update_settings_and_docs(ImportedSettings())
        self.settings.update(settings)
        super().__init__(network=network, **kwargs)
        # Transfer all geometrical properties off of network
        exclude.extend(['pore.coords', 'throat.conns'])
        for item in network.props():
            if item not in exclude:
                self[item] = network.pop(item)

        # If the following 'essential' props are not already defined, then
        # they should be added using the specified values or models

        if 'pore.diameter' not in self.keys():
            pdia = 'pore.'+self.settings['pore_diameter'].split('pore.')[-1]
            try:
                self['pore.diameter'] = self[pdia]
            except KeyError:
                logger.error(pdia + " not found, can't assign 'pore.diameter'")

        if 'pore.volume' not in self.keys():
            pore_shape = self.settings['pore_shape']
            m = getattr(mods.geometry.pore_volume, pore_shape)
            self.add_model(propname='pore.volume',
                           model=m, pore_diameter='pore.diameter')

        if 'pore.area' not in self.keys():
            pore_shape = self.settings['pore_shape']
            m = getattr(mods.geometry.pore_area, pore_shape)
            self.add_model(propname='pore.area',
                           model=m)

        if 'throat.diameter' not in self.keys():
            tdia = 'throat.'+self.settings['throat_diameter'].split('throat.')[-1]
            try:
                self['throat.diameter'] = self[tdia]
            except KeyError:
                logger.error(tdia + " not found, can't assign 'throat.diameter'")

        if 'throat.endpoints' not in self.keys():
            self.add_model(propname='throat.endpoints',
                           model=mods.geometry.throat_endpoints.spherical_pores,
                           pore_diameter='pore.diameter',
                           throat_diameter='throat.diameter')

        if 'throat.length' not in self.keys():
            self.add_model(propname='throat.length',
                           model=mods.geometry.throat_length.piecewise,
                           throat_endpoints='throat.endpoints')

        if 'throat.volume' not in self.keys():
            shape = self.settings['throat_shape']
            m = getattr(mods.geometry.throat_volume, shape)
            self.add_model(propname='throat.volume',
                           model=m,
                           throat_length='throat.length',
                           throat_diameter='throat.diameter')

        if 'throat.conduit_lengths' not in self.keys():
            self.add_model(propname='throat.conduit_lengths',
                           model=mods.geometry.throat_length.conduit_lengths,
                           throat_endpoints='throat.endpoints',
                           throat_length='throat.length')
