r"""
===============================================================================
Submodule -- generic_source_term
===============================================================================

"""

import scipy as _sp

def linear(physics,
           phase,
           A1,
           A2,
           x=None,
           **kwargs):
    r"""
    For the following source term:
        .. math::
            r = S_{1}   x  +  S_{2} 
    It returns the slope and intercept for the corresponding pores                 

    Parameters
    ----------
    
    Notes
    -----
    
    """  
    S1 = A1*_sp.ones(phase.Np)
    S2 = A2*_sp.ones(phase.Np)
    r = _sp.vstack((S1,S2)).T
    return r[physics.map_pores()]

def power_law(physics,
            phase,
            A1,
            A2,
            A3,
            x=None,        
            **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x^{A_{2}}  +  A_{3}
    It calculates the slope and intercept for the following linear form:                 
        .. math::            
            r = S_{1}   x  +  S_{2} 
    
    Parameters
    ----------
    
    Notes
    -----
    
    """
    if x is None:   
        X = _sp.zeros(phase.Np)
        length = _sp.size(X)
    else:
        length = _sp.size(x)
        if length==phase.Np or length>phase.Np:    X = x
        elif length<phase.Np:
            if _sp.size(x)==physics.Np:
                X = _sp.zeros(phase.Np)
                X[physics.map_pores()] = x
            else:   raise Exception('Wrong size for the guess!')     

    if _sp.size(A1) != length:    A1 = A1*_sp.ones_like(X)   
    if _sp.size(A2) != length:    A2 = A2*_sp.ones_like(X)   
    if _sp.size(A3) != length:    A3 = A3*_sp.ones_like(X)   
    S1 = A1*A2*X**(A2-1)
    S2 = A1*X**A2*(1-A2)+A3
    r = _sp.vstack((S1,S2)).T    
    return r[physics.map_pores()]

def exponential(physics,
                phase,
                A1,
                A2,
                A3,
                A4,
                A5,
                A6,
                x=None,        
                **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6} 
    It calculates the slope and intercept for the following linear form:                 
        .. math::            
            r = S_{1}   x  +  S_{2} 
    
    Parameters
    ----------
    
    Notes
    -----
    
    """       
    
    if x is None:   
        X = _sp.zeros(phase.Np)
        length = _sp.size(X)
    else:
        length = _sp.size(x)
        if length==phase.Np or length>phase.Np:    X = x
        elif length<phase.Np:
            if _sp.size(x)==physics.Np:
                X = _sp.zeros(phase.Np)
                X[physics.map_pores()] = x
            else:   raise Exception('Wrong size for the guess!')     

    if _sp.size(A1) != length:    A1 = A1*_sp.ones_like(X)   
    if _sp.size(A2) != length:    A2 = A2*_sp.ones_like(X)   
    if _sp.size(A3) != length:    A3 = A3*_sp.ones_like(X)   
    if _sp.size(A4) != length:    A4 = A4*_sp.ones_like(X)   
    if _sp.size(A5) != length:    A5 = A5*_sp.ones_like(X)   
    if _sp.size(A6) != length:    A6 = A6*_sp.ones_like(X)   

    S1 = A1*A3*A4*_sp.log(A2)*A2**(A3*X**A4+A5)*X**(A4-1)
    S2 = A1*A2**(A3*X**A4+A5)*(1-A3*A4*_sp.log(A2)*X**A4)+A6
    r = _sp.vstack((S1,S2)).T    
    return r[physics.map_pores()]


def natural_exponential(physics,
                        phase,
                        A1,
                        A2,
                        A3,
                        A4,
                        A5,
                        x=None,        
                        **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5} 
    It calculates the slope and intercept for the following linear form:                 
        .. math::            
            r = S_{1}   x  +  S_{2} 
    
    Parameters
    ----------
    
    Notes
    -----
    
    """       
    
    if x is None:   
        X = _sp.zeros(phase.Np)
        length = _sp.size(X)
    else:
        length = _sp.size(x)
        if length==phase.Np or length>phase.Np:    X = x
        elif length<phase.Np:
            if _sp.size(x)==physics.Np:
                X = _sp.zeros(phase.Np)
                X[physics.map_pores()] = x
            else:   raise Exception('Wrong size for the guess!')     

    if _sp.size(A1) != length:    A1 = A1*_sp.ones_like(X)   
    if _sp.size(A2) != length:    A2 = A2*_sp.ones_like(X)   
    if _sp.size(A3) != length:    A3 = A3*_sp.ones_like(X)   
    if _sp.size(A4) != length:    A4 = A4*_sp.ones_like(X)   
    if _sp.size(A5) != length:    A5 = A5*_sp.ones_like(X)   

    S1 = A1*A2*A3*X**(A3-1)*_sp.exp(A2*X**A3+A4)
    S2 = A1*_sp.exp(A2*X**A3+A4)*(1-A2*A3*X**A3)+A5
    r = _sp.vstack((S1,S2)).T    
    return r[physics.map_pores()]
    
    
def logarithm(physics,
              phase,
              A1,
              A2,
              A3,
              A4,
              A5,
              A6,
              x=None,        
              **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1}   Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}  
    It calculates the slope and intercept for the following linear form:                 
        .. math::            
            r = S_{1}   x  +  S_{2} 
    
    Parameters
    ----------
    
    Notes
    -----
    
    """       
    
    if x is None:   
        X = _sp.zeros(phase.Np)
        length = _sp.size(X)
    else:
        length = _sp.size(x)
        if length==phase.Np or length>phase.Np:    X = x
        elif length<phase.Np:
            if _sp.size(x)==physics.Np:
                X = _sp.zeros(phase.Np)
                X[physics.map_pores()] = x
            else:   raise Exception('Wrong size for the guess!')     

    if _sp.size(A1) != length:    A1 = A1*_sp.ones_like(X)   
    if _sp.size(A2) != length:    A2 = A2*_sp.ones_like(X)   
    if _sp.size(A3) != length:    A3 = A3*_sp.ones_like(X)   
    if _sp.size(A4) != length:    A4 = A4*_sp.ones_like(X)   
    if _sp.size(A5) != length:    A5 = A5*_sp.ones_like(X)   
    if _sp.size(A6) != length:    A6 = A6*_sp.ones_like(X)   

    S1 = A1*A3*A4*X**(A4-1)/(_sp.log(A2)*(A3*X**A4+A5))
    S2 = A1*_sp.log(A3*X**A4+A5)/_sp.log(A2)+A6-A1*A3*A4*X**A4/(_sp.log(A2)*(A3*X**A4+A5))
    r = _sp.vstack((S1,S2)).T    
    return r[physics.map_pores()]
    
def natural_logarithm(physics,
                      phase,
                      A1,
                      A2,
                      A3,
                      A4,
                      A5,
                      x=None,        
                      **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}    
    It calculates the slope and intercept for the following linear form:                 
        .. math::            
            r = S_{1}   x  +  S_{2} 
    
    Parameters
    ----------
    
    Notes
    -----
    
    """       
    
    if x is None:   
        X = _sp.zeros(phase.Np)
        length = _sp.size(X)
    else:
        length = _sp.size(x)
        if length==phase.Np or length>phase.Np:    X = x
        elif length<phase.Np:
            if _sp.size(x)==physics.Np:
                X = _sp.zeros(phase.Np)
                X[physics.map_pores()] = x
            else:   raise Exception('Wrong size for the guess!')     

    if _sp.size(A1) != length:    A1 = A1*_sp.ones_like(X)   
    if _sp.size(A2) != length:    A2 = A2*_sp.ones_like(X)   
    if _sp.size(A3) != length:    A3 = A3*_sp.ones_like(X)   
    if _sp.size(A4) != length:    A4 = A4*_sp.ones_like(X)   
    if _sp.size(A5) != length:    A5 = A5*_sp.ones_like(X)   

    S1 = A1*A2*A3*X**(A3-1)/(A2*X**A3+A4)
    S2 = A1*_sp.log(A2*X**A3+A4)+A5-A1*A2*A3*X**A3/(A2*X**A3+A4)
    r = _sp.vstack((S1,S2)).T    
    return r[physics.map_pores()]