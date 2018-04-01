import scipy as _sp


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


def linear(target, A1='', A2='', quantity=''):
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

    quantity : string
        The property name for the main quantity being produced or consumed.

    Returns
    -------
    A dictionary containing the values for *S1*, *S2*, and *rate*:

        d = {'pore.S1': vals_1, 'pore.S2': vals_2, 'pore.rate': vals_3}

    Starting in V2, OpenPNM allows the assigning of such ``dicts`` to a Core
    object ``dict``, but it breaks the received ``dict`` into separate arrays
    and prepends the assigned propname.  For instance:

        obj['pore.source_term']  = d

    is equivalent to:

        obj['pore.source_term_S1'] = d['pore.S1']
        obj['pore.source_term_S2'] = d['pore.S2']
        obj['pore.source_term_rate'] = d['pore.rate']

    """
    X = target[quantity]
    r = target[A1] * X + target[A2]
    S1 = target[A1]
    S2 = target[A2]
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def power_law(target, A1='', A2='', A3='', quantity=''):
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

    quantity : string
        The property name for the main quantity being produced or consumed.

    Returns
    -------
    A dictionary containing the values for *S1*, *S2*, and *rate*:

        d = {'pore.S1': vals_1, 'pore.S2': vals_2, 'pore.rate': vals_3}

    Starting in V2, OpenPNM allows the assigning of such ``dicts`` to a Core
    object ``dict``, but it breaks the received ``dict`` into separate arrays
    and prepends the assigned propname.  For instance:

        obj['pore.source_term']  = d

    is equivalent to:

        obj['pore.source_term_S1'] = d['pore.S1']
        obj['pore.source_term_S2'] = d['pore.S2']
        obj['pore.source_term_rate'] = d['pore.rate']

    """
    X = target[quantity]

    r = target[A1] * X ** target[A2] + target[A3]
    S1 = target[A1] * target[A2] * X ** (target[A2] - 1)
    S2 = target[A1] * X ** target[A2] * (1 - target[A2]) + target[A3]
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def exponential(target, A1='', A2='', A3='', A4='', A5='', A6='',
                quantity='', return_rate=True):
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

    quantity : string
        The property name for the main quantity being produced or consumed.

    Returns
    -------
    A dictionary containing the values for *S1*, *S2*, and *rate*:

        d = {'pore.S1': vals_1, 'pore.S2': vals_2, 'pore.rate': vals_3}

    Starting in V2, OpenPNM allows the assigning of such ``dicts`` to a Core
    object ``dict``, but it breaks the received ``dict`` into separate arrays
    and prepends the assigned propname.  For instance:

        obj['pore.source_term']  = d

    is equivalent to:

        obj['pore.source_term_S1'] = d['pore.S1']
        obj['pore.source_term_S2'] = d['pore.S2']
        obj['pore.source_term_rate'] = d['pore.rate']

    """
    X = target[quantity]
    a = {'1': A1, '2': A2, '3': A3, '4': A4, '5': A5, '6': A6}
    for item in a.keys():
        if a[item] == '':
            a[item] = 0
        else:
            a[item] = target[a[item]]

    r = a['1'] * a['2'] ** (a['3'] * X ** a['4'] + a['5']) + a['6']
    S1 = a['1'] * a['3'] * a['4'] * X ** (a['4'] - 1) * _sp.log(a['2']) *\
        a['2'] ** (a['3'] * X ** a['4'] + a['5'])
    S2 = a['1'] * a['2'] ** (a['3'] * X ** a['4'] + a['5']) * \
        (1 - a['3'] * a['4'] * _sp.log(a['2']) * X ** a['4']) + a['6']
    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def natural_exponential(target, A1='', A2='', A3='', A4='', A5='',
                        quantity=''):
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

    quantity : string
        The property name for the main quantity being produced or consumed.

    Returns
    -------
    A dictionary containing the values for *S1*, *S2*, and *rate*:

        d = {'pore.S1': vals_1, 'pore.S2': vals_2, 'pore.rate': vals_3}

    Starting in V2, OpenPNM allows the assigning of such ``dicts`` to a Core
    object ``dict``, but it breaks the received ``dict`` into separate arrays
    and prepends the assigned propname.  For instance:

        obj['pore.source_term']  = d

    is equivalent to:

        obj['pore.source_term_S1'] = d['pore.S1']
        obj['pore.source_term_S2'] = d['pore.S2']
        obj['pore.source_term_rate'] = d['pore.rate']

    """
    X = target[quantity]
    a = {'1': A1, '2': A2, '3': A3, '4': A4, '5': A5}
    for item in a.keys():
        if a[item] == '':
            a[item] = 0
        else:
            a[item] = target[a[item]]

    r = a['1'] * _sp.exp(a['2'] * X ** a['3'] + a['4']) + a['5']

    S1 = a['1'] * a['2'] * a['3'] * X ** (a['3'] - 1) * \
        _sp.exp(a['2'] * X ** a['3'] + a['4'])

    S2 = a['1'] * (1 - a['2'] * a['3'] * X ** a['3']) * \
        _sp.exp(a['2'] * X ** a['3'] + a['4']) + a['5']

    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def logarithm(target, A1='', A2='', A3='', A4='', A5='', A6='', quantity=''):
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

    quantity : string
        The property name for the main quantity being produced or consumed.

    Returns
    -------
    A dictionary containing the values for *S1*, *S2*, and *rate*:

        d = {'pore.S1': vals_1, 'pore.S2': vals_2, 'pore.rate': vals_3}

    Starting in V2, OpenPNM allows the assigning of such ``dicts`` to a Core
    object ``dict``, but it breaks the received ``dict`` into separate arrays
    and prepends the assigned propname.  For instance:

        obj['pore.source_term']  = d

    is equivalent to:

        obj['pore.source_term_S1'] = d['pore.S1']
        obj['pore.source_term_S2'] = d['pore.S2']
        obj['pore.source_term_rate'] = d['pore.rate']

    """
    X = target[quantity]
    a = {'1': A1, '2': A2, '3': A3, '4': A4, '5': A5, '6': A6}
    for item in a.keys():
        if a[item] == '':
            a[item] = 0
        else:
            a[item] = target[a[item]]

    r = a['1'] * _sp.log(a['3'] * X ** a['4'] + a['5']) / _sp.log(a['2']) + \
        a['6']

    S1 = a['1'] * a['3'] * a['4'] * X ** (a['4'] - 1) / \
        (_sp.log(a['2']) * (a['3'] * X ** a['4'] + a['5']))

    S2 = a['1'] * _sp.log(a['3'] * X ** a['4'] + a['5']) / _sp.log(a['2']) + \
        a['6'] - a['1'] * a['3'] * a['4'] * X ** a['4'] / \
        (_sp.log(a['2']) * (a['3'] * X ** a['4'] + a['5']))

    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values


def natural_logarithm(target, A1='', A2='', A3='', A4='', A5='', quantity=''):
    r"""
    For the following source term:

        .. math::
            r =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}

    This model calculates the rate of generation or consumption based on the
    current value of ``quantity``, as well as the slope and interecpt of the
    above equation at th current value of ``quantity`` giving the following
    linear form:

        .. math::
            r = S_{1}   x  +  S_{2}

    which can be used in linear solvers to produce a better guess of
    ``quantity``, which can

    Parameters
    ----------
    A1 -> A5 : string
        The property name of the coefficients in the source term model

    quantity : string
        The property name for the main quantity being produced or consumed.

    Returns
    -------
    A dictionary containing the values for *S1*, *S2*, and *rate*:

        d = {'pore.S1': vals_1, 'pore.S2': vals_2, 'pore.rate': vals_3}

    Starting in V2, OpenPNM allows the assigning of such ``dicts`` to a Core
    object ``dict``, but it breaks the received ``dict`` into separate arrays
    and prepends the assigned propname.  For instance:

        obj['pore.source_term']  = d

    is equivalent to:

        obj['pore.source_term_S1'] = d['pore.S1']
        obj['pore.source_term_S2'] = d['pore.S2']
        obj['pore.source_term_rate'] = d['pore.rate']

    """
    X = target[quantity]
    a = {'1': A1, '2': A2, '3': A3, '4': A4, '5': A5}
    for item in a.keys():
        if a[item] == '':
            a[item] = 0
        else:
            a[item] = target[a[item]]

    r = a['1'] * _sp.log(a['2'] * X ** a['3'] + a['4']) + a['5']

    S1 = a['1'] * a['2'] * a['3'] * X ** (a['3'] - 1) / \
        (a['2'] * X ** a['3'] + a['4'])

    S2 = a['1'] * _sp.log(a['2'] * X ** a['3'] + a['4']) + \
        a['5'] - a['1'] * a['2'] * a['3'] * X ** a['3'] / \
        (a['2'] * X ** a['3'] + a['4'])

    values = {'S1': S1, 'S2': S2, 'rate': r}
    return values
