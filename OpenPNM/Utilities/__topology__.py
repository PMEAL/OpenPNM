# -*- coding: utf-8 -*-
"""
===============================================================================
Utilities.topology: Assorted topological manipulation methods
===============================================================================

"""
from OpenPNM.Network import tools
from OpenPNM.Base import logging as _logging
logger = _logging.getLogger(__name__)


class topology(object):

    @staticmethod
    def extend(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.extend(network=network, **kwargs)

    @staticmethod
    def trim(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.trim(network=network, **kwargs)

    @staticmethod
    def clone_pores(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.clone_pores(network=network, **kwargs)

    @staticmethod
    def stitch(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.stitch(network=network, **kwargs)

    @staticmethod
    def connect_pores(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.connect_pores(network=network, **kwargs)

    @staticmethod
    def find_centroid(coords, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        c = tools.find_centroid(coords, **kwargs)
        return c

    @staticmethod
    def find_pores_distance(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        dist = tools.find_pores_distance(network=network, **kwargs)
        return dist

    @staticmethod
    def subdivide(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.subdivide(network=network, **kwargs)

    @staticmethod
    def trim_occluded_throats(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.trim_occluded_throats(network=network, **kwargs)

    @staticmethod
    def merge_pores(network, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.merge_pores(network=network, **kwargs)

    @staticmethod
    def template_sphere_shell(**kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        img = tools.template_sphere_shell(**kwargs)
        return img
