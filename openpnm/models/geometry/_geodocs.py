try:
    from matplotlib._docstring import Substitution
except ModuleNotFoundError:
    from matplotlib.docstring import Substitution


__all__ = [
    '_geodocs',
]


_geodocs = Substitution(
    network=
    r"""network : OpenPNM Network object

    """,
    Dp=
    r"""pore_diameter : str
            Name of the dictionary key on ``network`` containing the ndarray of
            pore diameter values.""",
    Dt=
    r"""throat_diameter : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat diameter values.""",
    Vp=
    r"""pore_volume : str
            Name of the dictionary key on ``network`` containing the ndarray of
            pore volume values.""",
    Vt=
    r"""throat_volume : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat volume values.""",
    Lt=
    r"""throat_length : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat length values.""",
    Pcoords=
    r"""pore_coords : str
            Name of the dictionary key on ``network`` containing the ndarray of
            pore coordinate values.""",
    Tcoords=
    r"""throat_coords : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat centroid coordinate values.""",
    At=
    r"""throat_area : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat surface values.""",
    Pt=
    r"""throat_perimeter : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat perimeter values.""",
    Tcen=
    r"""throat_centroid : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat centroid coordinate values.""",
    Act=
    r"""throat_cross_sectional_area : str
            Name of the dictionary key on ``network`` containing the ndarray of
            throat cross-sectional area values.""",
)
