r"""
===============================================================================
Submodule -- generic_source_term
===============================================================================

"""

import scipy as _sp


def linear(physics, phase, A1='', A2='', x='', return_rate=True, **kwargs):
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
    if x == '':
        X = _sp.ones(physics.Np) * _sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = physics[x]
            except KeyError:
                raise Exception(physics.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    length_X = _sp.size(X)
    if length_X != physics.Np:
        if length_X == 1:
            X = X * _sp.ones(physics.Np)
        elif length_X >= phase.Np:
            X = X[physics.map_pores()]
        else:
            raise Exception('Wrong size for the numerical array of x!')

    a = {}
    source_params = [A1, A2]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A == '':
            a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = physics[A]
                except KeyError:
                    raise Exception(physics.name + '/' + phase.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    if return_rate:
        return(a['1'] * X ** a['2'])
    else:
        S1 = a['1']
        S2 = a['2']
        return(_sp.vstack((S1, S2)).T)


def power_law(physics, phase, A1='', A2='', A3='', x='',
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
    if x == '':
        X = _sp.ones(physics.Np) * _sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = physics[x]
            except KeyError:
                raise Exception(physics.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    length_X = _sp.size(X)
    if length_X != physics.Np:
        if length_X == 1:
            X = X * _sp.ones(physics.Np)
        elif length_X >= phase.Np:
            X = X[physics.map_pores()]
        else:
            raise Exception('Wrong size for the numerical array of x!')

    a = {}
    source_params = [A1, A2, A3]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A == '':
            a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = physics[A]
                except KeyError:
                    raise Exception(physics.name + '/' + phase.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    if return_rate:
        return(a['1'] * X ** a['2'] + a['3'])
    else:
        S1 = a['1'] * a['2'] * X ** (a['2'] - 1)
        S2 = a['1'] * X ** a['2'] * (1 - a['2']) + a['3']
        return(_sp.vstack((S1, S2)).T)


def exponential(physics, phase, A1='', A2='', A3='', A4='', A5='', A6='',
                x='', return_rate=True, **kwargs):
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
    if x == '':
        X = _sp.ones(physics.Np) * _sp.nan
    else:
        if type(x) == str:
            x = 'pore.'+x.split('.')[-1]
            try:
                X = physics[x]
            except KeyError:
                raise Exception(physics.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    length_X = _sp.size(X)
    if length_X != physics.Np:
        if length_X == 1:
            X = X * _sp.ones(physics.Np)
        elif length_X >= phase.Np:
            X = X[physics.map_pores()]
        else:
            raise Exception('Wrong size for the numerical array of x!')

    a = {}
    source_params = [A1, A2, A3, A4, A5, A6]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A == '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = physics[A]
                except KeyError:
                    raise Exception(physics.name + '/' + phase.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    if return_rate:
        return a['1'] * a['2'] ** (a['3'] * X ** a['4'] + a['5']) + a['6']
    else:
        S1 = a['1'] * a['3'] * a['4'] * \
            X ** (a['4'] - 1) * _sp.log(a['2']) * \
            a['2'] ** (a['3'] * X ** a['4'] + a['5'])
        S2 = a['1'] * a['2'] ** (a['3'] * X ** a['4'] + a['5']) * \
            (1 - a['3'] * a['4'] * _sp.log(a['2']) * X ** a['4']) + a['6']
        return(_sp.vstack((S1, S2)).T)


def natural_exponential(physics, phase, A1='', A2='', A3='', A4='', A5='',
                        x='', return_rate=True, **kwargs):
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
    if x == '':
        X = _sp.ones(physics.Np)*_sp.nan
    else:
        if type(x) == str:
            x = 'pore.'+x.split('.')[-1]
            try:
                X = physics[x]
            except KeyError:
                raise Exception(physics.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    length_X = _sp.size(X)
    if length_X != physics.Np:
        if length_X == 1:
            X = X * _sp.ones(physics.Np)
        elif length_X >= phase.Np:
            X = X[physics.map_pores()]
        else:
            raise Exception('Wrong size for the numerical array of x!')

    a = {}
    source_params = [A1, A2, A3, A4, A5]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A == '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = physics[A]
                except KeyError:
                    raise Exception(physics.name + '/' + phase.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    if return_rate:
        return(a['1'] * _sp.exp(a['2'] * X ** a['3'] + a['4']) + a['5'])
    else:
        S1 = a['1'] * a['2'] * \
            a['3'] * X ** (a['3'] - 1) * \
            _sp.exp(a['2'] * X ** a['3'] + a['4'])
        S2 = a['1'] * (1 - a['2'] * a['3'] * X ** a['3']) * \
            _sp.exp(a['2'] * X ** a['3'] + a['4']) + a['5']
        return(_sp.vstack((S1, S2)).T)


def logarithm(physics, phase, A1='', A2='', A3='', A4='', A5='', A6='',
              x='', return_rate=True, **kwargs):
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
    if x == '':
        X = _sp.ones(physics.Np)*_sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = physics[x]
            except KeyError:
                raise Exception(physics.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    length_X = _sp.size(X)
    if length_X != physics.Np:
        if length_X == 1:
            X = X * _sp.ones(physics.Np)
        elif length_X >= phase.Np:
            X = X[physics.map_pores()]
        else:
            raise Exception('Wrong size for the numerical array of x!')

    a = {}
    source_params = [A1, A2, A3, A4, A5, A6]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A == '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = physics[A]
                except KeyError:
                    raise Exception(physics.name + '/' + phase.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    if return_rate:
        return(a['1'] * _sp.log(a['3'] * X ** a['4'] + a['5']) /
               _sp.log(a['2']) + a['6'])
    else:
        S1 = a['1'] * a['3'] * a['4'] * \
            X ** (a['4'] - 1) / \
            (_sp.log(a['2']) * (a['3'] * X ** a['4'] + a['5']))
        S2 = a['1'] * _sp.log(a['3'] * X ** a['4'] + a['5']) / \
            _sp.log(a['2']) + a['6'] - a['1'] * a['3'] * \
            a['4'] * X ** a['4'] / \
            (_sp.log(a['2']) * (a['3'] * X ** a['4'] + a['5']))
        return(_sp.vstack((S1, S2)).T)


def natural_logarithm(physics, phase, A1='', A2='', A3='', A4='', A5='',
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
    if x == '':
        X = _sp.ones(physics.Np)*_sp.nan
    else:
        if type(x) == str:
            x = 'pore.' + x.split('.')[-1]
            try:
                X = physics[x]
            except KeyError:
                raise Exception(physics.name +
                                ' does not have the pore property :' + x + '!')
        else:
            X = _sp.array(x)

    length_X = _sp.size(X)

    if length_X != physics.Np:
        if length_X == 1:
            X = X * _sp.ones(physics.Np)
        elif length_X >= phase.Np:
            X = X[physics.map_pores()]
        else:
            raise Exception('Wrong size for the numerical array of x!')

    a = {}
    source_params = [A1, A2, A3, A4, A5]
    for ind in _sp.arange(_sp.size(source_params)):
        A = source_params[ind]
        if A == '':
            if ind == 0:
                a[str(ind+1)] = 1
            else:
                a[str(ind+1)] = 0
        else:
            if type(A) == str:
                A = 'pore.' + A.split('.')[-1]
                try:
                    a[str(ind+1)] = physics[A]
                except KeyError:
                    raise Exception(physics.name + '/' + phase.name +
                                    ' does not have the pore property :' +
                                    A + '!')
            else:
                raise Exception('source_term parameters can only be string '
                                'type!')

    if return_rate:
        return(a['1'] * _sp.log(a['2'] * X ** a['3'] + a['4']) + a['5'])
    else:
        S1 = a['1'] * a['2'] * a['3'] * \
            X ** (a['3'] - 1) / \
            (a['2'] * X ** a['3'] + a['4'])
        S2 = a['1'] * _sp.log(a['2'] * X ** a['3'] + a['4']) + \
            a['5'] - a['1'] * a['2'] * a['3'] * \
            X ** a['3'] / (a['2'] * X ** a['3'] + a['4'])
        return(_sp.vstack((S1, S2)).T)
