# -*- coding: utf-8 -*-

from collections import namedtuple as _namedtuple


class PhysicalConstants:
    r"""
    A generic object contain the intrinsic properties of the Phase
    """
    molecular_weight = None
    critical_pressure = None
    critical_temperature = None
    melting_point = None
    freezing_point = None
    antoines_coefficients = _namedtuple('antoines_coefficients',
                                        ('A', 'B', 'C'))
