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

    def extend(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.extend(**kwargs)

    def trim(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.trim(**kwargs)

    def clone_pores(**kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.clone_pores(**kwargs)

    def stitch(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.stitch(**kwargs)

    def connect_pores(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.connect_pores(**kwargs)

    def find_centroid(**kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.find_centroid(**kwargs)

    def find_pores_distance(**kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.find_pores_distance(**kwargs)

    def subdivide(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.subdivide(**kwargs)

    def trim_occluded_throats(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.trim_occluded_throats(**kwargs)

    def merge_pores(self, **kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.merge_pores(**kwargs)

    def template_sphere_shell(**kwargs):
        r"""
        This function as been moved to ``Network.tools`` and remains here for
        backward compatibility.
        """
        tools.template_sphere_shell(**kwargs)
