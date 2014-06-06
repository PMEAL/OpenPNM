r"""
===============================================================================
Submodule -- throat_volume
===============================================================================

"""
import scipy as sp

def constant(geometry,
             network,
             propname,
             value,
             **params):
    r"""
    Assigns specified constant value
    """
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cylinder(geometry,
             network,
             propname,
             throat_length='length',
             throat_diameter='diameter',
             **params):
    r"""
    Calculate throat diameter from seeds for a cylindrical throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    Tlen = network.get_data(prop=throat_length,throats=geometry.throats())
    Tdia = network.get_data(prop=throat_diameter,throats=geometry.throats())
    value = sp.pi/4*Tlen*Tdia**2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def cuboid(geometry,
           network,
           propname,
           throat_length='length',
           throat_diameter='diameter',
           **params):
    r"""
    Calculate throat volume of cuboidal throat
    - note: this will need to account for volume taken up by spherical pore bodies
    """
    Tlen = network.get_data(prop=throat_length,throats=geometry.throats())
    Tdia = network.get_data(prop=throat_diameter,throats=geometry.throats())
    value = Tlen*Tdia**2
    network.set_data(prop=propname,throats=geometry.throats(),data=value)

def voronoi(geometry,
            network,
            propname,
            throat_length='length',
            throat_area='area',
            **params):
    r"""
    Calculate volume from the voronoi facet area and the throat length
    """
    value=network.get_throat_data(prop=throat_length,locations=geometry)*network.get_throat_data(prop=throat_area,locations=geometry)
    network.set_data(prop=propname,throats=geometry.throats(),data=value)