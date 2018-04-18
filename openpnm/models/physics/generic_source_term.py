r"""
===============================================================================
Submodule -- generic_source_term
===============================================================================

"""

import scipy as _sp
import sympy as syp
import collections


def standard_kinetics(target, quantity, prefactor, exponent):
    r"""

    """
    X = target[quantity]
    A = target[prefactor]
    b = target[exponent]

    r = A*(X**b)
    S1 = A*b*(X**(b - 1))
    S2 = A*(1 - b)*(X**b)
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def linear(target, A1='', A2='', x='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x  +  A_{2}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 , A2 : string
        The property name of the coefficients in the source term model.
        With A2 set to zero this equation takes on the familiar for of r=kx.
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----
    Because this source term is linear in concentration (x) is it not necessary
    to iterate during the solver step.  Thus, when using the
    ``set_source_term`` method for an algorithm, it is recommended to set the
    ``maxiter``
    argument to 0.  This will save 1 unncessary solution of the system, since
    the solution would coverge after the first pass anyway.

    """
    if x is '':
        X = _sp.ones(target.Np) * _sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = target[x]
            except KeyError:
                raise Exception(target.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    a = {}
    source_params = [A1, A2]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A is '':
            a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = target[A]
                except KeyError:
                    raise Exception(target.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    r = a['1'] * X + a['2']
    S1 = a['1']
    S2 = a['2']
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def power_law(target, A1='', A2='', A3='', x='',
              return_rate=True, **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x^{A_{2}}  +  A_{3}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A3 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    if x is '':
        X = _sp.ones(target.Np) * _sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = target[x]
            except KeyError:
                raise Exception(target.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    a = {}
    source_params = [A1, A2, A3]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A is '':
            a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = target[A]
                except KeyError:
                    raise Exception(target.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    r=a['1'] * X ** a['2'] + a['3']
    S1 = a['1'] * a['2'] * X ** (a['2'] - 1)
    S2 = a['1'] * X ** a['2'] * (1 - a['2']) + a['3']
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def exponential(target, A1='', A2='', A3='', A4='', A5='', A6='',
                x='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A6 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    if x is '':
        X = _sp.ones(target.Np) * _sp.nan
    else:
        if type(x) == str:
            x = 'pore.'+x.split('.')[-1]
            try:
                X = target[x]
            except KeyError:
                raise Exception(target.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    a = {}
    source_params = [A1, A2, A3, A4, A5, A6]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A is '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = target[A]
                except KeyError:
                    raise Exception(target.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    r = a['1'] * a['2'] ** (a['3'] * X ** a['4'] + a['5']) + a['6']
    S1 = a['1'] * a['3'] * a['4'] * \
        X ** (a['4'] - 1) * _sp.log(a['2']) * \
        a['2'] ** (a['3'] * X ** a['4'] + a['5'])
    S2 = a['1'] * a['2'] ** (a['3'] * X ** a['4'] + a['5']) * \
        (1 - a['3'] * a['4'] * _sp.log(a['2']) * X ** a['4']) + a['6']
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def natural_exponential(target, A1='', A2='', A3='', A4='', A5='',
                        x='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A5 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    if x is '':
        X = _sp.ones(target.Np)*_sp.nan
    else:
        if type(x) == str:
            x = 'pore.'+x.split('.')[-1]
            try:
                X = target[x]
            except KeyError:
                raise Exception(target.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    a = {}
    source_params = [A1, A2, A3, A4, A5]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A is '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = target[A]
                except KeyError:
                    raise Exception(target.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    r = a['1'] * _sp.exp(a['2'] * X ** a['3'] + a['4']) + a['5']
    S1 = a['1'] * a['2'] * \
        a['3'] * X ** (a['3'] - 1) * \
        _sp.exp(a['2'] * X ** a['3'] + a['4'])
    S2 = a['1'] * (1 - a['2'] * a['3'] * X ** a['3']) * \
        _sp.exp(a['2'] * X ** a['3'] + a['4']) + a['5']
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def logarithm(target, A1='', A2='', A3='', A4='', A5='', A6='',
              x='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1}   Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A6 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    if x is '':
        X = _sp.ones(target.Np)*_sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = target[x]
            except KeyError:
                raise Exception(target.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    a = {}
    source_params = [A1, A2, A3, A4, A5, A6]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A is '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = target[A]
                except KeyError:
                    raise Exception(target.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    r = (a['1'] * _sp.log(a['3'] * X ** a['4'] + a['5']) /
               _sp.log(a['2']) + a['6'])
    S1 = a['1'] * a['3'] * a['4'] * \
        X ** (a['4'] - 1) / \
        (_sp.log(a['2']) * (a['3'] * X ** a['4'] + a['5']))
    S2 = a['1'] * _sp.log(a['3'] * X ** a['4'] + a['5']) / \
        _sp.log(a['2']) + a['6'] - a['1'] * a['3'] * \
        a['4'] * X ** a['4'] / \
        (_sp.log(a['2']) * (a['3'] * X ** a['4'] + a['5']))
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values


def natural_logarithm(target, A1='', A2='', A3='', A4='', A5='',
                      x='', return_rate=True, **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A5 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    if x is '':
        X = _sp.ones(target.Np)*_sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = target[x]
            except KeyError:
                raise Exception(target.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    a = {}
    source_params = [A1, A2, A3, A4, A5]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A is '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = target[A]
                except KeyError:
                    raise Exception(target.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    r = a['1'] * _sp.log(a['2'] * X ** a['3'] + a['4']) + a['5']
    S1 = a['1'] * a['2'] * a['3'] * \
        X ** (a['3'] - 1) / \
        (a['2'] * X ** a['3'] + a['4'])
    S2 = a['1'] * _sp.log(a['2'] * X ** a['3'] + a['4']) + \
        a['5'] - a['1'] * a['2'] * a['3'] * \
        X ** a['3'] / (a['2'] * X ** a['3'] + a['4'])
    values = {'pore.S1': S1, 'pore.S2': S2, 'pore.rate': r}
    return values

# Symbols used in all symbolic functions except general_symbolic which gets its
# own : a-f are coefficients and x is the independent variable
a, b, c, d, e, f, x = syp.symbols('a,b,c,d,e,f,x')


def build_func(eq, args=None):
    r'''
    Take a symbolic equation and return the lambdified version plus the
    Linearization of form S1 * x + S2
    '''
    eq_prime = eq.diff(x)
    s1 = eq_prime
    s2 = eq - eq_prime*x
    EQ = syp.lambdify(args, eq, 'numpy')
    S1 = syp.lambdify(args, s1, 'numpy')
    S2 = syp.lambdify(args, s2, 'numpy')
    return EQ, S1, S2


def linear_sym(target, A1='', A2='', X='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x  +  A_{2}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 , A2 : string
        The property name of the coefficients in the source term model.
        With A2 set to zero this equation takes on the familiar for of r=kx.
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----
    Because this source term is linear in concentration (x) is it not necessary
    to iterate during the solver step.  Thus, when using the
    ``set_source_term`` method for an algorithm, it is recommended to set the
    ``maxiter``
    argument to 0.  This will save 1 unncessary solution of the system, since
    the solution would coverge after the first pass anyway.

    """
    A = target['pore.' + A1.split('.')[-1]]
    B = target['pore.' + A2.split('.')[-1]]
    X = target['pore.' + X.split('.')[-1]]
    # Equation
    y = a*x + b
    # Callable functions
    r, s1, s2 = build_func(y, (a, b, x))
    # Values
    r_val=r(A, B, X)
    s1_val=s1(A, B, X)
    s2_val=s2(A, B, X)
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values


def power_law_sym(target, A1='', A2='', A3='', X='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x^{A_{2}}  +  A_{3}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A3 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    A = target['pore.' + A1.split('.')[-1]]
    B = target['pore.' + A2.split('.')[-1]]
    C = target['pore.' + A3.split('.')[-1]]
    X = target['pore.' + X.split('.')[-1]]
    # Equation
    y =  a*x**b + c
    # Callable functions
    r, s1, s2 = build_func(y, (a, b, c, x))
    # Values
    r_val=r(A, B, C, X)
    s1_val=s1(A, B, C, X)
    s2_val=s2(A, B, C, X)
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values


def exponential_sym(target, A1='', A2='', A3='',
                    A4='', A5='', A6='', X='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A6 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    A = target['pore.' + A1.split('.')[-1]]
    B = target['pore.' + A2.split('.')[-1]]
    C = target['pore.' + A3.split('.')[-1]]
    D = target['pore.' + A4.split('.')[-1]]
    E = target['pore.' + A5.split('.')[-1]]
    F = target['pore.' + A6.split('.')[-1]]
    X = target['pore.' + X.split('.')[-1]]
    # Equation
    y =  a*b**(c*x**d + e)+f
    # Callable functions
    r, s1, s2 = build_func(y, (a, b, c, d, e, f, x))
    # Values
    r_val=r(A, B, C, D, E, F, X)
    s1_val=s1(A, B, C, D, E, F, X)
    s2_val=s2(A, B, C, D, E, F, X)
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values


def natural_exponential_sym(target, A1='', A2='', A3='',
                            A4='', A5='', X='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A5 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    A = target['pore.' + A1.split('.')[-1]]
    B = target['pore.' + A2.split('.')[-1]]
    C = target['pore.' + A3.split('.')[-1]]
    D = target['pore.' + A4.split('.')[-1]]
    E = target['pore.' + A5.split('.')[-1]]
    X = target['pore.' + X.split('.')[-1]]
    # Equation
    y =  a*syp.exp(b*x**c + d)+e
    # Callable functions
    r, s1, s2 = build_func(y, (a, b, c, d, e, x))
    # Values
    r_val=r(A, B, C, D, E, X)
    s1_val=s1(A, B, C, D, E, X)
    s2_val=s2(A, B, C, D, E, X)
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values


def logarithm_sym(target, A1='', A2='', A3='',
                  A4='', A5='', A6='', X='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1}   Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A6 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    A = target['pore.' + A1.split('.')[-1]]
    B = target['pore.' + A2.split('.')[-1]]
    C = target['pore.' + A3.split('.')[-1]]
    D = target['pore.' + A4.split('.')[-1]]
    E = target['pore.' + A5.split('.')[-1]]
    F = target['pore.' + A6.split('.')[-1]]
    X = target['pore.' + X.split('.')[-1]]
    # Equation
    y =  a*syp.log((c*x**d + e), b)+f
    # Callable functions
    r, s1, s2 = build_func(y, (a, b, c, d, e, f, x))
    # Values
    r_val=r(A, B, C, D, E, F, X)
    s1_val=s1(A, B, C, D, E, F, X)
    s2_val=s2(A, B, C, D, E, F, X)
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values


def natural_logarithm_sym(target, A1='', A2='', A3='',
                          A4='', A5='', X='', **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}
    If return_rate is True, it returns the value of source term for the
    provided x in each pore.
    If return_rate is False, it calculates the slope and intercept for the
    following linear form :
        .. math::
            r = S_{1}   x  +  S_{2}

    Parameters
    ----------
    A1 -> A5 : string
        The property name of the coefficients in the source term model
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----

    """
    A = target['pore.' + A1.split('.')[-1]]
    B = target['pore.' + A2.split('.')[-1]]
    C = target['pore.' + A3.split('.')[-1]]
    D = target['pore.' + A4.split('.')[-1]]
    E = target['pore.' + A5.split('.')[-1]]
    X = target['pore.' + X.split('.')[-1]]
    # Equation
    y =  a*syp.ln(b*x**c + d)+e
    # Callable functions
    r, s1, s2 = build_func(y, (a, b, c, d, e, x))
    # Values
    r_val=r(A, B, C, D, E, X)
    s1_val=s1(A, B, C, D, E, X)
    s2_val=s2(A, B, C, D, E, X)
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values


def general_symbolic(target, eqn=None, arg_map=None, **kwargs):
    r'''
    A general function to interpret a sympy equation and evaluate the linear
    components of the source term.

    Parameters
    ----------
    target : OpenPNM object
    eqn: sympy symbolic expression for the source terms
        e.g. y = a*x**b + c
    arg_map : Dict mapping the symbols in the expression to OpenPNM data
        on the target. Must contain 'x' which is the independent variable.
        e.g.
        arg_map={'a':'pore.a', 'b':'pore.b', 'c':'pore.c', 'x':'pore.x'}

    Example
    ----------
    >>> import openpnm as op
    >>> from openmpnm.models.physics import generic_source_term as gst
    >>> import scipy as sp
    >>> import sympy as syp
    >>> pn = op.network.Cubic(shape=[5, 5, 5], spacing=0.0001)
    >>> water = op.phases.Water(network=pn)
    >>> water['pore.a'] = 1
    >>> water['pore.b'] = 2
    >>> water['pore.c'] = 3
    >>> water['pore.x'] = sp.random.random(water.Np)
    >>> a,b,c,x = syp.symbols('a,b,c,x')
    >>> y = a*x**b + c
    >>> arg_map={'a':'pore.a', 'b':'pore.b', 'c':'pore.c', 'x':'pore.x'}
    >>> water.add_model(propname='pore.general',
    ...                 model=gst.general_symbolic,
    ...                 eqn=y, arg_map=arg_map,
    ...                 regen_mode='normal')
    >>> assert 'pore.general.rate' in water.props()
    >>> assert 'pore.general.S1' in water.props()
    >>> assert 'pore.general.S1' in water.props()
    '''
    # First make sure all the symbols have been allocated dict items
    for arg in syp.postorder_traversal(eqn):
        if syp.srepr(arg)[:6] == 'Symbol':
            key = syp.srepr(arg)[7:].strip('(').strip(')').strip("'")
            if key not in arg_map.keys():
                raise Exception('argument mapping incomplete, missing '+key)
    if 'x' not in arg_map.keys():
        raise Exception('argument mapping must contain "x" for the '+
                        'independent variable')
    # Get the data
    data = collections.OrderedDict()
    for key in arg_map.keys():
        data[key] = target[arg_map[key]]
        # Callable functions
    symbols = tuple(arg_map.keys())
    r, s1, s2 = build_func(eqn, symbols)
    r_val=r(*data.values())
    s1_val=s1(*data.values())
    s2_val=s2(*data.values())
    values = {'pore.S1': s1_val, 'pore.S2': s2_val, 'pore.rate': r_val}
    return values
