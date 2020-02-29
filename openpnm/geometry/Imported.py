import openpnm.models as mods
from openpnm.geometry import GenericGeometry


defset = {'pore_diameter': 'equivalent_diameter',
          'throat_diameter': 'equivalent_diameter'}

# The following will appear as the "help" docstring for the settings attribute
s = r"""
    The following table lists the various settings on this object and
    provides a brief description of their meaning.

    ================  =========================================================
    pore_diameter     Key into the extracted data array to use as pore
                      diameter in other geometry calculations. The default is
                      'pore.equivalent_diameter'.  Use of 'pore.' is not
                      required.
    ----------------  ---------------------------------------------------------
    throat_diameter   Key into the extracted data array to use as throat
                      diameter in other geometry calculations. The default is
                      'throat.equivalent_diameter'.  Use of 'throat.' is not
                      required.
    ================  =========================================================

    """


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
        super().__init__(network=network, **kwargs)
        self.settings.update(defset)
        self.settings.__doc__ = s
        self.settings.update(settings)
        exclude.extend(['pore.coords', 'throat.conns'])
        for item in network.props():
            if item not in exclude:
                self[item] = network.pop(item)

        if 'pore.diameter' not in self.keys():
            pdia = 'pore.'+self.settings['pore_diameter'].split('pore.')[-1]
            self['pore.diameter'] = self[pdia]

        if 'throat.diameter' not in self.keys():
            tdia = 'throat.'+self.settings['throat_diameter'].split('thraot.')[-1]
            self['throat.diameter'] = self[tdia]

        if 'throat.endpoints' not in self.keys():
            self.add_model(propname='throat.endpoints',
                           model=mods.geometry.throat_endpoints.spherical_pores,
                           pore_diameter='pore.diameter',
                           throat_diameter='throat.diameter')

        if 'throat.length' not in self.keys():
            self.add_model(propname='throat.length',
                           model=mods.geometry.throat_length.piecewise,
                           throat_endpoints='throat.endpoints')

        self.add_model(propname='throat.conduit_lengths',
                       model=mods.geometry.throat_length.conduit_lengths,
                       throat_endpoints='throat.endpoints',
                       throat_length='throat.length')
