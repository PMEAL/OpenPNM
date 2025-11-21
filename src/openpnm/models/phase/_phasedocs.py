try:
    from matplotlib._docstring import Substitution
except ModuleNotFoundError:
    from matplotlib.docstring import Substitution


__all__ = [
    '_phasedocs',
]


_phasedocs = Substitution(
    phase=
    r"""phase : OpenPNM Phase object
            The phase object to which this model is associated (i.e. attached).
            This controls the length of the calculated array(s), and also
            provides access to other necessary properties.""",
    T=
    r"""T : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            temperature values, in units of [K].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    P=
    r"""P : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            pressure values, in units of [Pa].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    MWs=
    r"""MWs : str (dict key) or dict of numeric values for each species

    """,
    MW=
    r"""MW : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            molecular weight values, in units of [g/mol].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    Tcs=
    r"""Tcs : str (dict key) or scalar

    """,
    Tc=
    r"""Tc : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            critical temperature values, in units of [K].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    Pcs=
    r"""Pcs : str (dict key) or scalar

    """,
    Pc=
    r"""Pc : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            critical pressure values, in units of [Pa].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    Tb=
    r"""Tb : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            boiling temperature values, in units of [K].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    Vcs=
    r"""Vcs : str (dict key) or scalar

    """,
    Vc=
    r"""Vc : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            critical volume values, in units of [m^3].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    omegas=
    r"""omegas : str (dict key) or scalar

    """,
    omega=
    r"""omega : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            accentric factor values, which is dimensionless.
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    epsilons=
    r"""epsilons : str (dict key) or scalar

    """,
    epsilon=
    r"""epsilon : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            Lennard-Jones epsilon values, in units of [K].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    sigmas=
    r"""sigmas : str (dict key) or scalar

    """,
    sigma=
    r"""sigma : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            Lennard-Jones molecular diameter values, in units of [A].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    Cpg=
    r"""Cpg : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            gas phase heat capacity values, in units of [J/mol].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    Cps=
    r"""Cps : str (dict key) or scalar

    """,
    ks=
    r"""ks : str (dict key) or numeric
    """,
    Vdms=
    r"""

    """,
    Vm=
    r"""Vm : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            molar volume, in units of [m3/mol].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    n_V=
    r"""n_V : str (dict key) or scalar

    """,
    rhos=
    r"""rhos : str
            The dictionary key containing the density in kg/m3""",
    rho=
    r"""rho : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            mass density values, in units of [kg/m^3].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    mus=
    r"""mus : str

    """,
    mu=
    r"""mu : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            kinematic viscosity values, in units of [Pa.s].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    SInote=
    r"""Since numpy arrays don't natively support units, OpenPNM cannot enforce
    units on values and arrays, however the use of SI is assumed.""",
    salinity=
    r"""salinity : str
            The dictionary key containing the salinity values in g/kg""",
    conc=
    r"""conc : str (dict key) or numeric
            Name of the dictionary key on ``phase`` containing the ndarray of
            concentration values, in units of [mol/L].
            If a numerical value is passed (i.e. a scalar or ndarray)
            it gets used directly.""",
    zs=
    r"""T : str (dict key)
            Name of the dictionary key on ``phase`` containing the ndarray of
            mole fraction values.""",
)
