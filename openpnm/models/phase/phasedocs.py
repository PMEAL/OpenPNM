from matplotlib.docstring import Substitution


__all__ = [
    '_phasedocs',
]


_phasedocs = Substitution(
    phase=
    r"""phase : OpenPNM Phase object

    """,
    T=
    r"""T : str (dict key) or scalar
            Name of the dictionary key on ``phase`` containing the array of
            temperature values. If a scalar is passed it get used directly.""",
    P=
    r"""P : str (dict key) or scalar
            Name of the dictionary key on ``phase`` containing the array of
            pressure values. If a scalar is passed it get used directly.""",
    MWs=
    r"""MWs : dict key (str) or scalar

    """,
    MW=
    r"""MW : dict key (str) or scalar

    """,
    Tcs=
    r"""Tcs : str (dict key) or scalar

    """,
    Tc=
    r"""Tc : str (dict key) or scalar

    """,
    Pcs=
    r"""Pcs : str (dict key) or scalar

    """,
    Pc=
    r"""Pc : str (dict key) or scalar

    """,
    Tb=
    r"""Tb : str (dict key) or scalar

    """,
    Vcs=
    r"""Vcs : str (dict key) or scalar

    """,
    Vc=
    r"""Vc : str (dict key) or scalar

    """,
    omegas=
    r"""omegas : str (dict key) or scalar

    """,
    omega=
    r"""omega : str (dict key) or scalar

    """,
    epsilons=
    r"""epsilons : str (dict key) or scalar

    """,
    epsilon=
    r"""epsilon : str (dict key) or scalar

    """,
    sigmas=
    r"""sigmas : str (dict key) or scalar

    """,
    sigma=
    r"""sigma : str (dict key) or scalar

    """,
    Cpg=
    r"""Cpg : str (dict key) or scalar

    """,
    Cps=
    r"""Cps : str (dict key) or scalar

    """,
    ks=
    r"""ks : str (dict key) or scalar

    """,
    Vdms=
    r"""

    """,
    rhos=
    r"""rhos : str
            The dictionary key containing the density in kg/m3""",
    rho=
    r"""rho : str
            The dictionary key containing the density in kg/m3""",
    mu=
    r"""mu : str

    """,
    SInote=
    r"""Since numpy arrays don't natively support units, OpenPNM cannot enforce
    units on values and arrays, however the use of SI is assumed.""",
    salinity=
    r"""salinity : str
            The dictionary key containing the salinity values in g/kg""",
)
