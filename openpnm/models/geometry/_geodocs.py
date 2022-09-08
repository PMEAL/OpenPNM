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
            Name of the dictionary key on ``target`` where the array containing
            pore diameter values is stored""",
    Dt=
    r"""throat_diameter : str
            Name of the dictionary key on ``target`` where the array containing
            pore diameter values is stored""",
    Vp=
    r"""pore_volume : str
            Name of the dictionary key on ``target`` where the array containing
            pore volume values is stored""",
    Vt=
    r"""throat_volume : str
            Name of the dictionary key on ``target`` where the array containing
            throat volume values is stored""",
    Lt=
    r"""throat_length : str
            Name of the dictionary key on ``target`` where the array containing
            throat length values is stored""",
    Pcoords=
    r"""pore_coords : str

    """,
    Tcoords=
    r"""throat_coords : str

    """,
    At=
    r"""throat_area : str

    """,
    Pt=
    r"""throat_perimeter : str

    """,
    Dcen=
    r"""throat_centroid : str

    """,
)
