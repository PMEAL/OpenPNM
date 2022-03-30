import numpy as _np
import scipy as _sp
from openpnm.models import _doctxt


__all__ = [
    "charge_conservation",
    "standard_kinetics",
    "linear",
    "power_law",
    "exponential",
    "natural_exponential",
    "logarithm",
    "natural_logarithm",
    "general_symbolic",
    "butler_volmer_conc",
    "butler_volmer_voltage"
]


@_doctxt
def charge_conservation(target, phase, p_alg, e_alg, assumption):
    r"""
    Applies the source term on the charge conservation equation when solving
    for ions transport.

    Parameters
    ----------
    %(target_blurb)s
    p_alg : GenericAlgorithm
        The algorithm used to enforce charge conservation.
    e_alg : list
        The list of algorithms used to solve for transport of different
        ionic species of the mixture phase.
    assumption : str
        The assumption adopted to enforce charge conservation. Options are:

        ================= ====================================================
        Options           Description
        ================= ====================================================
        poisson           ?
        electroneutrality ?
        laplace           ?
        ================= ====================================================

    Returns
    -------
    rate_info : dict
        A dictionary containing the following three items:

        ======= ==============================================================
        Item    Description
        ======= ==============================================================
        rate    The value of the source term function for the given list
                of algortihms under the provided assumption.
        S1      A placeholder (zeros array)
        S2      The value of the source term function for the given list of
                algortihms under the provided assumption (same as 'rate').
        ======= ==============================================================

    """
    assumption = assumption.lower()
    import scipy.sparse.csgraph as _spgr

    F = 96485.3321233100184
    rhs = _np.zeros(shape=(p_alg.Np, ), dtype=float)
    network = p_alg.project.network
    if assumption == 'poisson':
        v = network['pore.volume']
        for e in e_alg:
            rhs += (v * F * phase['pore.valence.'+e.settings['ion']]
                    * e[e.settings['quantity']])
    elif assumption == 'poisson_2d':
        s = network['pore.cross_sectional_area']
        for e in e_alg:
            rhs += (s * F * phase['pore.valence.'+e.settings['ion']]
                    * e[e.settings['quantity']])
    elif assumption in ['electroneutrality', 'electroneutrality_2d']:
        for e in e_alg:
            try:
                c = e[e.settings['quantity']]
            except KeyError:
                c = _np.zeros(shape=(e.Np, ), dtype=float)
            network = e.project.network
            g = phase['throat.diffusive_conductance.'+e.settings['ion']]
            am = network.create_adjacency_matrix(weights=g, fmt='coo')
            A = _spgr.laplacian(am)
            rhs += - F * phase['pore.valence.'+e.settings['ion']] * (A * c)
    elif assumption in ['laplace', 'laplace_2d']:
        pass  # rhs should remain 0
    else:
        raise Exception('Unknown keyword for charge_conservation, pick from:'
                        + ' poisson, poisson_2d, laplace, laplace2d,'
                        + ' electroneutrality or electroneutrality_2d')
    S1 = _np.zeros(shape=(p_alg.Np, ), dtype=float)
    values = {'S1': S1, 'S2': rhs, 'rate': rhs}
    return values


@_doctxt
def standard_kinetics(target, X, prefactor, exponent):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of ``X``:

        .. math::
            r = A X^b

    Parameters
    ----------
    %(target_blurb)s
    X : str
        %(dict_blurb)s quantity of interest
    prefactor : str
        %(dict_blurb)s the prefactor to be used in the source term model
    exponent : str
        %(dict_blurb)s the exponent to be used in the source term model

    Returns
    -------
    rate_info : dict
        A dictionary containing the following three items:

        ======= ==============================================================
        Item    Description
        ======= ==============================================================
        rate    The value of the source term function at the given X.
        S1      The slope of the source term function at the given X.
        S2      The intercept of the source term function at the given X.
        ======= ==============================================================

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    X = target[X]
    A = target[prefactor]
    b = target[exponent]

    r = A*(X**b)
    S1 = A*b*(X**(b - 1))
    S2 = A*(1 - b)*(X**b)
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def _parse_args(target, key, default):
    r"""
    Parses the passed ``key``: returns ``default`` if ``key`` is empty,
    otherwise fetches and returns the values from ``target``.
    """
    if key == '':
        val = default
    else:
        val = target[key]
    return val


@_doctxt
def linear(target, X, A1='', A2=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of ``X``:

        .. math::
            r = A_{1} X + A_{2}

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key on the target object containing the the quantity
        of interest
    A1 -> A2 : str
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    Returns
    -------
    dict
        A dictionary containing the following three items:

            'rate'
                The value of the source term function at the given X.
            'S1'
                The slope of the source term function at the given X.
            'S2'
                The intercept of the source term function at the given X.

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    A = _parse_args(target=target, key=A1, default=0.0)
    B = _parse_args(target=target, key=A2, default=0.0)
    X = target[X]

    r = A * X + B
    S1 = A
    S2 = B
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


@_doctxt
def power_law(target, X, A1='', A2='', A3=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *X*:

        .. math::
            r = A_{1} x^{A_{2}} + A_{3}

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key on the target object containing the the quantity
        of interest
    A1 -> A3 : str
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    Returns
    -------
    dict
        A dictionary containing the following three items:

            'rate'
                The value of the source term function at the given X.
            'S1'
                The slope of the source term function at the given X.
            'S2'
                The intercept of the source term function at the given X.

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    A = _parse_args(target=target, key=A1, default=0.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=0.0)
    X = target[X]

    r = A * X ** B + C
    S1 = A * B * X ** (B - 1)
    S2 = A * X ** B * (1 - B) + C
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


@_doctxt
def exponential(target, X, A1='', A2='', A3='', A4='', A5='', A6=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6}

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key on the target object containing the the quantity
        of interest
    A1 -> A6 : str
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    Returns
    -------
    dict
        A dictionary containing the following three items:

            'rate'
                The value of the source term function at the given X.
            'S1'
                The slope of the source term function at the given X.
            'S2'
                The intercept of the source term function at the given X.

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    A = _parse_args(target=target, key=A1, default=0.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    F = _parse_args(target=target, key=A6, default=0.0)
    X = target[X]

    r = A * B ** (C * X ** D + E) + F
    S1 = A * C * D * X ** (D - 1) * _np.log(B) * B ** (C * X ** D + E)
    S2 = A * B ** (C * X ** D + E) * (1 - C * D * _np.log(B) * X ** D) + F
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


@_doctxt
def natural_exponential(target, X, A1='', A2='', A3='', A4='', A5=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r = A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5}

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key on the target object containing the the quantity
        of interest
    A1 -> A5 : str
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    Returns
    -------
    dict
        A dictionary containing the following three items:

            'rate'
                The value of the source term function at the given X.
            'S1'
                The slope of the source term function at the given X.
            'S2'
                The intercept of the source term function at the given X.

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    A = _parse_args(target=target, key=A1, default=0.0)
    B = _parse_args(target=target, key=A2, default=0.0)
    C = _parse_args(target=target, key=A3, default=0.0)
    D = _parse_args(target=target, key=A4, default=0.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    X = target[X]

    r = A * _np.exp(B * X ** C + D) + E
    S1 = A * B * C * X ** (C - 1) * _np.exp(B * X ** C + D)
    S2 = A * (1 - B * C * X ** C) * _np.exp(B * X ** C + D) + E
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


@_doctxt
def logarithm(target, X, A1='', A2='', A3='', A4='', A5='', A6=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =  A_{1} Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key on the target object containing the the quantity
        of interest
    A1 -> A6 : str
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    Returns
    -------
    dict
        A dictionary containing the following three items:

            'rate'
                The value of the source term function at the given X.
            'S1'
                The slope of the source term function at the given X.
            'S2'
                The intercept of the source term function at the given X.

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    A = _parse_args(target=target, key=A1, default=0.0)
    B = _parse_args(target=target, key=A2, default=10.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    F = _parse_args(target=target, key=A6, default=0.0)
    X = target[X]

    r = (A * _np.log(C * X ** D + E)/_np.log(B) + F)
    S1 = A * C * D * X ** (D - 1) / (_np.log(B) * (C * X ** D + E))
    S2 = A * _np.log(C * X ** D + E) / _np.log(B) + F - A * C * D * X ** D / \
        (_np.log(B) * (C * X ** D + E))
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


@_doctxt
def natural_logarithm(target, X, A1='', A2='', A3='', A4='', A5=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =   A_{1} Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key on the target object containing the the quantity
        of interest
    A1 -> A5 : str
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    Returns
    -------
    dict
        A dictionary containing the following three items:

            'rate'
                The value of the source term function at the given X.
            'S1'
                The slope of the source term function at the given X.
            'S2'
                The intercept of the source term function at the given X.

    Notes
    -----
    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1} X + S_{2}

    """
    A = _parse_args(target=target, key=A1, default=0.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=0.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    X = target[X]

    r = A*_np.log(B*X**C + D) + E
    S1 = A*B*C*X**(C - 1) / (B * X ** C + D)
    S2 = A*_np.log(B*X**C + D) + E - A*B*C*X**C / (B*X**C + D)
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def _build_func(eq, **args):
    r"""
    Take a symbolic equation and return the lambdified version plus the
    linearization of form S1 * x + S2
    """
    from sympy import lambdify
    eq_prime = eq.diff(args['x'])
    s1 = eq_prime
    s2 = eq - eq_prime*args['x']
    EQ = lambdify(args.values(), expr=eq, modules='numpy')
    S1 = lambdify(args.values(), expr=s1, modules='numpy')
    S2 = lambdify(args.values(), expr=s2, modules='numpy')
    return EQ, S1, S2


@_doctxt
def general_symbolic(target, eqn, x, **kwargs):
    r"""
    A general function to interpret a sympy equation and evaluate the linear
    components of the source term.

    Parameters
    ----------
    %(target_blurb)s
    eqn : str
        The str representation of the equation to use.  This will be
        passed to sympy's ``sympify`` function to make a *live* sympy object.
    x : str
        The dictionary key of the independent variable
    kwargs
        All additional keyword arguments are converted to sympy variables
        using the ``symbols`` function.  Note that IF the arguments are
        strs, it is assumed they are dictionary keys pointing to arrays
        on the ``target`` object.  If they are numerical values they are
        used 'as is'.  Numpy arrays are not accepted.  These must be stored
        in the ``target`` dictionary and referenced by key.

    Examples
    --------
    >>> import openpnm as op
    >>> from openpnm.models.physics import generic_source_term as gst
    >>> import numpy as np
    >>> import sympy
    >>> pn = op.network.Cubic(shape=[5, 5, 5], spacing=0.0001)
    >>> water = op.phase.Water(network=pn)
    >>> water['pore.a'] = 1
    >>> water['pore.b'] = 2
    >>> water['pore.c'] = 3
    >>> water['pore.x'] = np.random.random(water.Np)
    >>> y = 'a*x**b + c'
    >>> arg_map = {'a':'pore.a', 'b':'pore.b', 'c':'pore.c'}
    >>> water.add_model(propname='pore.general',
    ...                 model=gst.general_symbolic,
    ...                 eqn=y, x='pore.x', **arg_map)

    """
    from sympy import symbols, sympify
    eqn = sympify(eqn)
    # Get the data
    data = {'x': target[x]}
    args = {'x': symbols('x')}
    for key in kwargs.keys():
        if isinstance(kwargs[key], str):
            data[key] = target[kwargs[key]]
        else:
            data[key] = kwargs[key]
        args[key] = symbols(key)
    r, s1, s2 = _build_func(eqn, **args)
    r_val = r(*data.values())
    s1_val = s1(*data.values())
    s2_val = s2(*data.values())
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


@_doctxt
def butler_volmer_conc(
    target, X, z, j0, c_ref, alpha_anode, alpha_cathode, z_rds=1,
    reaction_order=1,
    temperature="pore.temperature",
    reaction_area="pore.reaction_area",
    solid_voltage="pore.solid_voltage",
    electrolyte_voltage="pore.electrolyte_voltage",
    open_circuit_voltage="pore.open_circuit_voltage",
):
    r"""
    Calculates the rate, slope and intercept of the Butler-Volmer kinetic
    model based on **concentration** to be used in mass transfer
    algorithms.

    Parameters
    ----------
    %(target_blurb)s
    X : str
        The dictionary key of the quantity of interest (i.e. main variable
        to be solved; in this case, concentration).
    z : float
        Number of electrons transferred in the redox reaction.
    j0 : float
        Exchange current density [A/m^2].
    c_ref : float
        Reference concentration [mol/m^3].
    alpha_anode : float
        Anodic transfer coefficient.
    alpha_cathode : float
        Cathodic transfer coefficient.
    z_rds : float
        Number of electrons transferred in the rate determining step.
    reaction_order : float
        Reaction order, i.e. the exponent of the concentration term
    solid_voltage : str
        The dictionary key of solid phase voltages [V].
    electrolyte_voltage : str
        The dictionary key of electrolyte phase voltages [V].
    open_circuit_voltage : str
        The dictionary key of open-circuit voltage values [V].
    reaction_area : str
        The dictionary key of reaction area values [m^2].
    temperature : str
        The dictionary key of temperature values [K].

    Returns
    -------
    dict
        Dictionary containing the following key/value pairs:

        - rate : The value of the source term function at the given X.
        - S1 : The slope of the source term function at the given X.
        - S2 : The intercept of the source term function at the given X.

    Notes
    -----
    The difference between butler_volmer_conc and butler_volmer_voltage is
    that the former is linearized with respect to the electrolyte
    concentration whereas the latter is linearized with respect to the
    electrolyte voltage.

    Consequently, while the S1 and S2 produced by these model shouldn't be
    identical, they should both produce the same **rate** with the only
    difference that the rate generated by butler_volmer_conc has the units
    [mol/s] whereas that generated by butler_volmer_voltage has the units
    [C/s]. Therefore, the two rates will differ by z * F, where z is the
    number of electrons transferred and F is the Faraday's constant.

    .. math::
        r =  j_0 A_{rxn} (\frac{ 1 }{ z F })(\frac{ X }{ c_{ref} }) ^ {\nu}
        \Big(
            \exp(  \frac{\alpha_a z_{rds} F}{RT} \eta )
          - \exp( -\frac{\alpha_c z_{rds} F}{RT} \eta )
        \Big)

    where:

    .. math::
        \eta = V_s - V_{\ell} - V_{oc}

    where ``V_s`` is the solid voltage, ``V_l`` is the electrolyte voltage,
    and ``V_oc`` is the open-circuit voltage.

    The slope and intercept provide a linearized source term equation
    about the current value of X as follow:

    .. math::
        rate = S_{1} X + S_{2}
    
    Reference: Rahn, Christopher D. Wang, Chao-Yang. (2013). Battery Systems Engineering - 
    3.2.1 The Butler-Volmer Equation. (pp. 27). John Wiley & Sons.

    """
    network = target.project.network
    domain = target._domain
    pores = domain.pores(target.name)

    # Fetch model variables
    X = target[X]
    T = target[temperature]
    Vs = target[solid_voltage]
    Ve = target[electrolyte_voltage]
    Voc = target[open_circuit_voltage]
    A_rxn = network[reaction_area][pores]
    F = _sp.constants.physical_constants["Faraday constant"][0]
    R = _sp.constants.R
    nu = reaction_order

    # Linearize with respect to X (electrolyte concentration)
    eta = Vs - Ve - Voc
    cte = j0 * A_rxn / (z * F)
    m1 = alpha_anode * z_rds * F / (R * T)
    m2 = alpha_cathode * z_rds * F / (R * T)
    fV = _np.exp(m1 * eta) - _np.exp(-m2 * eta)
    fC = (X / c_ref)**nu
    r = cte * fC * fV
    drdC = cte * (X / c_ref)**(nu - 1) * (1 / c_ref) * fV
    S1 = drdC
    S2 = r - drdC * X

    values = {"pore.S1": S1, "pore.S2": S2, "pore.rate": r}
    return values


def butler_volmer_voltage(
    target, X, z, j0, c_ref, alpha_anode, alpha_cathode, z_rds=1,
    reaction_order=1,
    temperature="pore.temperature",
    reaction_area="pore.reaction_area",
    solid_voltage="pore.solid_voltage",
    open_circuit_voltage="pore.open_circuit_voltage",
    electrolyte_concentration="pore.electrolyte_concentration",
):
    r"""
    Calculates the rate, slope and intercept of the Butler-Volmer kinetic model
    based on **voltage** to be used in electron conduction algorithms.

    Parameters
    ----------
    target : GenericPhysics
        The Physics object where the result will be applied.
    X : str
        The dictionary key of the quantity of interest (i.e. main variable
        to be solved; in this case, electrolyte voltage).
    z : float
        Number of electrons transferred in the redox reaction.
    j0 : float
        Exchange current density [A/m^2].
    c_ref : float
        Reference concentration [mol/m^3].
    alpha_anode : float
        Anodic transfer coefficient.
    alpha_cathode : float
        Cathodic transfer coefficient.
    z_rds : float
        Number of electrons transferred in the rate determining step.
    electrolyte_concentration : str
        The dictionary key of the electrolyte concentrations [mol/m^3].
    reaction_order : float
        Reaction order, i.e. the exponent of the concentration term
    solid_voltage : str
        The dictionary key of solid phase voltages [V].
    electrolyte_voltage : str
        The dictionary key of electrolyte phase voltages [V].
    open_circuit_voltage : str
        The dictionary key of open-circuit voltage values [V].
    reaction_area : str
        The dictionary key of reaction area values [m^2].
    temperature : str
        The dictionary key of temperature values [K].

    Returns
    -------
    rate_info : dict
        Dictionary containing the following key/value pairs:

        - rate : The value of the source term function at the given X.
        - S1 : The slope of the source term function at the given X.
        - S2 : The intercept of the source term function at the given X.

    Notes
    -----
    The difference between butler_volmer_conc and butler_volmer_voltage is
    that the former is linearized with respect to the electrolyte
    concentration whereas the latter is linearized with respect to the
    electrolyte voltage.

    Consequently, while the S1 and S2 produced by these model shouldn't be
    identical, they should both produce the same **rate** with the only
    difference that the rate generated by butler_volmer_conc has the units
    [mol/s] whereas that generated by butler_volmer_voltage has the units
    [C/s]. Therefore, the two rates will differ by z * F, where z is the
    number of electrons transferred and F is the Faraday's constant.

    .. math::
        r = j_0 A_{rxn} (\frac{ c }{ c_{ref} }) ^ {\nu}
        \Big(
            \exp(  \frac{\alpha_a z_{rds} F}{RT} \eta )
          - \exp( -\frac{\alpha_c z_{rds} F}{RT} \eta )
        \Big)

    where:

    .. math::
        \eta = V_s - X - V_{oc}

    where ``V_s`` is the solid voltage, ``X`` is the electrolyte voltage,
    and ``V_oc`` is the open-circuit voltage.

    The slope and intercept provide a linearized source term equation
    about the current value of X as follow:

    .. math::
        rate = S_{1} X + S_{2}
    
    Reference: Rahn, Christopher D. Wang, Chao-Yang. (2013). Battery Systems Engineering - 
    3.2.1 The Butler-Volmer Equation. (pp. 27). John Wiley & Sons.
    """
    network = target.project.network
    domain = target._domain
    pores = domain.pores(target.name)

    # Fetch model variables
    A_rxn = network[reaction_area][pores]
    Vs = target[solid_voltage]
    Voc = target[open_circuit_voltage]
    c = target[electrolyte_concentration]
    T = target[temperature]
    X = target[X]
    F = _sp.constants.physical_constants["Faraday constant"][0]
    R = _sp.constants.R
    nu = reaction_order

    # Linearize with respect to X (electrolyte voltage)
    eta = Vs - X - Voc
    cte = j0 * A_rxn
    m1 = alpha_anode * z_rds * F / (R * T)
    m2 = alpha_cathode * z_rds * F / (R * T)
    fV = _np.exp(m1 * eta) - _np.exp(-m2 * eta)
    dfVdV = -(m1 * _np.exp(m1 * eta) + m2 * _np.exp(-m2 * eta))
    fC = (c / c_ref)**nu
    r = cte * fC * fV
    drdV = cte * fC * dfVdV
    S1 = drdV
    S2 = r - drdV * X

    values = {"pore.S1": S1, "pore.S2": S2, "pore.rate": r}
    return values
