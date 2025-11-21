try:
    from matplotlib._docstring import Substitution
except ModuleNotFoundError:
    from matplotlib.docstring import Substitution



__all__ = [
    '_doctxt',
]


_doctxt = Substitution(
    phase=
    r"""phase : OpenPNM Phase object
            The phase object to which this model is associated (i.e. attached).
            This controls the length of the calculated array(s), and also
            provides access to other necessary properties.""",
    network=
    r"""network : OpenPNM Network object
            The network object to which this model is associated (i.e.
            attached). This controls the length of the calculated array(s),
            and also provides access to other necessary properties.""",
    dict_blurb=
    r"""Name of the dictionary key on the target object pointing to the
    array containing values of """,
    return_arr=
    r"""values : ndarray
            A numpy ndarray containing the computed values of """,
)
