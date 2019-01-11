r"""

.. autofunction:: openpnm.models.physics.generic_source_term.standard_kinetics
.. autofunction:: openpnm.models.physics.generic_source_term.linear
.. autofunction:: openpnm.models.physics.generic_source_term.power_law
.. autofunction:: openpnm.models.physics.generic_source_term.exponential
.. autofunction:: openpnm.models.physics.generic_source_term.natural_exponential
.. autofunction:: openpnm.models.physics.generic_source_term.logarithm
.. autofunction:: openpnm.models.physics.generic_source_term.natural_logarithm
.. autofunction:: openpnm.models.physics.generic_source_term.general_symbolic

"""

import scipy as _sp
import sympy as _syp


def standard_kinetics(target, quantity, prefactor, exponent):
    r"""

    """
    X = target[quantity]
    A = target[prefactor]
    b = target[exponent]

    r = A*(X**b)
    S1 = A*b*(X**(b - 1))
    S2 = A*(1 - b)*(X**b)
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def _parse_args(target, key, default):
    if key == '':
        val = default
    else:
        val = target[key]
    return val


def linear(target, X, A1='', A2=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r = A_{1}   X  +  A_{2}

    Parameters
    ----------
    A1 -> A2 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=0.0)
    X = target[X]

    r = A * X + B
    S1 = A
    S2 = B
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def power_law(target, X, A1='', A2='', A3=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *X*:

        .. math::
            r = A_{1}   x^{A_{2}}  +  A_{3}

    Parameters
    ----------
    A1 -> A3 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=0.0)
    X = target[X]

    r = A * X ** B + C
    S1 = A * B * X ** (B - 1)
    S2 = A * X ** B * (1 - B) + C
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def exponential(target, X, A1='', A2='', A3='', A4='', A5='', A6=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6}

    Parameters
    ----------
    A1 -> A6 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    F = _parse_args(target=target, key=A6, default=0.0)
    X = target[X]

    r = A * B ** (C * X ** D + E) + F
    S1 = A * C * D * X ** (D - 1) * _sp.log(B) * B ** (C * X ** D + E)
    S2 = A * B ** (C * X ** D + E) * (1 - C * D * _sp.log(B) * X ** D) + F
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def natural_exponential(target, X, A1='', A2='', A3='', A4='', A5=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =   A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5}

    Parameters
    ----------
    A1 -> A5 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    X = target[X]

    r = A * _sp.exp(B * X ** C + D) + E
    S1 = A * B * C * X ** (C - 1) * _sp.exp(B * X ** C + D)
    S2 = A * (1 - B * C * X ** C) * _sp.exp(B * X ** C + D) + E
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def logarithm(target, X, A1='', A2='', A3='', A4='', A5='', A6=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =  A_{1}   Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}

    Parameters
    ----------
    A1 -> A6 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    F = _parse_args(target=target, key=A6, default=0.0)
    X = target[X]

    r = (A * _sp.log(C * X ** D + E)/_sp.log(B) + F)
    S1 = A * C * D * X ** (D - 1) / (_sp.log(B) * (C * X ** D + E))
    S2 = A * _sp.log(C * X ** D + E) / _sp.log(B) + F - A * C * D * X ** D / \
        (_sp.log(B) * (C * X ** D + E))
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def natural_logarithm(target, X, A1='', A2='', A3='', A4='', A5=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of `X`:

        .. math::
            r =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}

    Parameters
    ----------
    A1 -> A5 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    X = target[X]

    r = A*_sp.log(B*X**C + D) + E
    S1 = A*B*C*X**(C - 1) / (B * X ** C + D)
    S2 = A*_sp.log(B*X**C + D) + E - A*B*C*X**C / (B*X**C + D)
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def _build_func(eq, **args):
    r'''
    Take a symbolic equation and return the lambdified version plus the
    linearization of form S1 * x + S2
    '''
    eq_prime = eq.diff(args['x'])
    s1 = eq_prime
    s2 = eq - eq_prime*args['x']
    EQ = _syp.lambdify(args.values(), expr=eq, modules='numpy')
    S1 = _syp.lambdify(args.values(), expr=s1, modules='numpy')
    S2 = _syp.lambdify(args.values(), expr=s2, modules='numpy')
    return EQ, S1, S2


def linear_sym(target, X, A1='', A2=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *x*:

        .. math::
            r = A_{1}   x  +  A_{2}

    Parameters
    ----------
    A1 -> A2 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target object containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=0.0)
    X = target[X]
    # Symbols used in symbolic function
    a, b, x = _syp.symbols('a,b,x')
    # Equation
    y = a*x + b
    # Callable functions
    r, s1, s2 = _build_func(eq=y, a=a, b=b,  x=x)
    # Values
    r_val = r(A, B, X)
    s1_val = s1(A, B, X)
    s2_val = s2(A, B, X)
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


def power_law_sym(target, X, A1='', A2='', A3=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *x*:

        .. math::
            r = A_{1}   x^{A_{2}}  +  A_{3}

    Parameters
    ----------
    A1 -> A3 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=0.0)
    X = target[X]
    # Symbols used in symbolic function
    a, b, c, x = _syp.symbols('a,b,c,x')
    # Equation
    y = a*x**b + c
    # Callable functions
    r, s1, s2 = _build_func(eq=y, a=a, b=b, c=c, x=x)
    # Values
    r_val = r(A, B, C, X)
    s1_val = s1(A, B, C, X)
    s2_val = s2(A, B, C, X)
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


def exponential_sym(target, X, A1='', A2='', A3='', A4='', A5='', A6=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *x*:

        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6}

    Parameters
    ----------
    A1 -> A6 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string or float/int or array/list
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    F = _parse_args(target=target, key=A6, default=0.0)
    X = target[X]
    # Symbols used in symbolic function
    a, b, c, d, e, f, x = _syp.symbols('a,b,c,d,e,f,x')
    # Equation
    y = a*b**(c*x**d + e) + f
    # Callable functions
    r, s1, s2 = _build_func(eq=y, a=a, b=b, c=c, d=d, e=e, f=f, x=x)
    # Values
    r_val = r(A, B, C, D, E, F, X)
    s1_val = s1(A, B, C, D, E, F, X)
    s2_val = s2(A, B, C, D, E, F, X)
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


def natural_exponential_sym(target, X, A1='', A2='', A3='', A4='', A5=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *x*:

        .. math::
            r =   A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5}

    Parameters
    ----------
    A1 -> A6 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string or float/int or array/list
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    X = target[X]
    # Symbols used in symbolic function
    a, b, c, d, e, x = _syp.symbols('a,b,c,d,e,x')
    # Equation
    y = a*_syp.exp(b*x**c + d) + e
    # Callable functions
    r, s1, s2 = _build_func(eq=y, a=a, b=b, c=c, d=d, e=e, x=x)
    # Values
    r_val = r(A, B, C, D, E, X)
    s1_val = s1(A, B, C, D, E, X)
    s2_val = s2(A, B, C, D, E, X)
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


def logarithm_sym(target, X, A1='', A2='', A3='', A4='', A5='', A6=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *x*:

        .. math::
            r =  A_{1}   Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}

    Parameters
    ----------
    A1 -> A6 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    x : string or float/int or array/list
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    F = _parse_args(target=target, key=A6, default=0.0)
    X = target[X]
    # Symbols used in symbolic function
    a, b, c, d, e, f, x = _syp.symbols('a,b,c,d,e,f,x')
    # Equation
    y = a*_syp.log((c*x**d + e), b) + f
    # Callable functions
    r, s1, s2 = _build_func(eq=y, a=a, b=b, c=c, d=d, e=e, f=f, x=x)
    # Values
    r_val = r(A, B, C, D, E, F, X)
    s1_val = s1(A, B, C, D, E, F, X)
    s2_val = s2(A, B, C, D, E, F, X)
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


def natural_logarithm_sym(target, X, A1='', A2='', A3='', A4='', A5=''):
    r"""
    Calculates the rate, as well as slope and intercept of the following
    function at the given value of *x*:

        .. math::
            rate =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}

    Parameters
    ----------
    A1 -> A5 : string
        The dictionary keys on the target object containing the coefficients
        values to be used in the source term model

    X : string or float/int or array/list
        The dictionary key on the target objecxt containing the the quantity
        of interest

    Returns
    -------
    A dictionary containing the following three items:

        **'rate'** - The value of the source term function at the given X.

        **'S1'** - The slope of the source term function at the given X.

        **'S2'** - The intercept of the source term function at the given X.

    The slope and intercept provide a linearized source term equation about the
    current value of X as follow:

        .. math::
            rate = S_{1}   X  +  S_{2}

    """
    A = _parse_args(target=target, key=A1, default=1.0)
    B = _parse_args(target=target, key=A2, default=1.0)
    C = _parse_args(target=target, key=A3, default=1.0)
    D = _parse_args(target=target, key=A4, default=1.0)
    E = _parse_args(target=target, key=A5, default=0.0)
    X = target[X]
    # Symbols used in symbolic function
    a, b, c, d, e, x = _syp.symbols('a,b,c,d,e,x')
    # Equation
    y = a*_syp.ln(b*x**c + d) + e
    # Callable functions
    r, s1, s2 = _build_func(eq=y, a=a, b=b, c=c, d=d, e=e, x=x)
    # Values
    r_val = r(A, B, C, D, E, X)
    s1_val = s1(A, B, C, D, E, X)
    s2_val = s2(A, B, C, D, E, X)
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values


def general_symbolic(target, eqn=None, arg_map=None):
    r'''
    A general function to interpret a sympy equation and evaluate the linear
    components of the source term.

    Parameters
    ----------
    target : OpenPNM object
        The OpenPNM object where the result will be applied.

    eqn : sympy symbolic expression for the source terms
        e.g. y = a*x**b + c

    arg_map : Dict mapping the symbols in the expression to OpenPNM data
        on the target. Must contain 'x' which is the independent variable.
        e.g. arg_map={'a':'pore.a', 'b':'pore.b', 'c':'pore.c', 'x':'pore.x'}

    Example
    ----------
    >>> import openpnm as op
    >>> from openpnm.models.physics import generic_source_term as gst
    >>> import scipy as sp
    >>> import sympy as _syp
    >>> pn = op.network.Cubic(shape=[5, 5, 5], spacing=0.0001)
    >>> water = op.phases.Water(network=pn)
    >>> water['pore.a'] = 1
    >>> water['pore.b'] = 2
    >>> water['pore.c'] = 3
    >>> water['pore.x'] = sp.random.random(water.Np)
    >>> a, b, c, x = _syp.symbols('a,b,c,x')
    >>> y = a*x**b + c
    >>> arg_map = {'a':'pore.a', 'b':'pore.b', 'c':'pore.c', 'x':'pore.x'}
    >>> water.add_model(propname='pore.general',
    ...                 model=gst.general_symbolic,
    ...                 eqn=y, arg_map=arg_map,
    ...                 regen_mode='normal')
    >>> assert 'pore.general.rate' in water.props()
    >>> assert 'pore.general.S1' in water.props()
    >>> assert 'pore.general.S1' in water.props()
    '''
    # First make sure all the symbols have been allocated dict items
    for arg in _syp.postorder_traversal(eqn):
        if _syp.srepr(arg)[:6] == 'Symbol':
            key = _syp.srepr(arg)[7:].strip('(').strip(')').strip("'")
            if key not in arg_map.keys():
                raise Exception('argument mapping incomplete, missing '+key)
    if 'x' not in arg_map.keys():
        raise Exception('argument mapping must contain "x" for the ' +
                        'independent variable')
    # Get the data
    data = {}
    args = {}
    for key in arg_map.keys():
        data[key] = target[arg_map[key]]
        # Callable functions
        args[key] = _syp.symbols(key)
    r, s1, s2 = _build_func(eqn, **args)
    r_val = r(*data.values())
    s1_val = s1(*data.values())
    s2_val = s2(*data.values())
    values = {'S1': s1_val, 'S2': s2_val, 'rate': r_val}
    return values
