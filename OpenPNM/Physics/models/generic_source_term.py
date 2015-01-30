r"""
===============================================================================
Submodule -- generic_source_term
===============================================================================

"""

import scipy as _sp

def linear(physics,
           phase,
           A1=0,
           A2=0,
           x=None,
           **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x  +  A_{2} 
    It returns the slope and intercept for the corresponding pores                 

    Parameters
    ----------
    A1, A2 : float or array
        These are the kinetic constants to be applied.  With A2 set to zero
        this equation takes on the familiar for of r=kx, where k is A1.  
        
    Notes
    -----
    It is possible to change the values of these constants at a later point.
    The values are stored in the Physics object's models dictionary under the
    chosen ``propname``.  The individual constants can be acesses as:
    phys.models[propname][A1].    
    
    Because this source term is linear in concentration (x) is it not necessary
    to iterate during the solver step.  Thus, when using the 
    ``set_source_term`` method it is recommended to set the ``maxiter`` 
    argument to 0.  This will save 1 unncessary solution of the system, since
    the solution would coverge after the first pass anyway.  
    
    """
    length_A1 = _sp.size(A1)
    length_A2 = _sp.size(A2)

    if length_A1==physics.Np:
        S1 = A1
    elif length_A1==1:
        S1 = A1*_sp.ones(physics.Np)        
    elif length_A1>=phase.Np:
        S1 = A1[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A1!') 

    if length_A2==physics.Np:
        S2 = A2
    elif length_A2==1:
        S2 = A2*_sp.ones(physics.Np)        
    elif length_A2>=phase.Np:
        S2 = A2[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A2!') 
    
    r = _sp.vstack((S1,S2)).T
    return r

def power_law(physics,
            phase,
            A1=0,
            A2=0,
            A3=0,
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
    A1 -> A3 : float or array
        The numerical values of the source term model's constant coefficients
        
    Notes
    -----
    It is possible to change the values of these constants at a later point.
    The values are stored in the Physics object's models dictionary under the
    chosen ``propname``.  The individual constants can be acesses as:
    phys.models[propname][A1].
    
    """
    if x is None:   
        X = _sp.zeros(physics.Np)
        length_X = _sp.size(x)
    else:
        length_X = _sp.size(x)
        if length_X==physics.Np:
            X = x
        elif length_X==1:
            X = x*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = x[physics.map_pores()] 
        else:   raise Exception('Wrong size for the guess value!')     

    length_A1 = _sp.size(A1)
    length_A2 = _sp.size(A2)
    length_A3 = _sp.size(A3)

    if length_A1==physics.Np:
        a1 = A1
    elif length_A1==1:
        a1 = A1*_sp.ones(physics.Np)        
    elif length_A1>=phase.Np:
        a1 = A1[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A1!') 

    if length_A2==physics.Np:
        a2 = A2
    elif length_A2==1:
        a2 = A2*_sp.ones(physics.Np)        
    elif length_A2>=phase.Np:
        a2 = A2[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A2!') 

    if length_A3==physics.Np:
        a3 = A3
    elif length_A3==1:
        a3 = A3*_sp.ones(physics.Np)        
    elif length_A3>=phase.Np:
        a3 = A3[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A3!') 
 
    S1 = a1*a2*X**(a2-1)
    S2 = a1*X**a2*(1-a2)+a3
    r = _sp.vstack((S1,S2)).T    
    return r

def exponential(physics,
                phase,
                A1=0,
                A2=0,
                A3=0,
                A4=0,
                A5=0,
                A6=0,
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
    A1 -> A6 : float or array
        The numerical values of the source term model's constant coefficients
    
    Notes
    -----
    It is possible to change the values of these constants at a later point.
    The values are stored in the Physics object's models dictionary under the
    chosen ``propname``.  The individual constants can be acesses as:
    phys.models[propname][A1].
    
    """        
    
    if x is None:   
        X = _sp.zeros(physics.Np)
        length_X = _sp.size(x)
    else:
        length_X = _sp.size(x)
        if length_X==physics.Np:
            X = x
        elif length_X==1:
            X = x*_sp.ones(physics.Np)  
        elif length_X>=phase.Np:
            X = x[physics.map_pores()] 
        else:   raise Exception('Wrong size for the guess value!')     

    length_A1 = _sp.size(A1)
    length_A2 = _sp.size(A2)
    length_A3 = _sp.size(A3)
    length_A4 = _sp.size(A4)
    length_A5 = _sp.size(A5)
    length_A6 = _sp.size(A6)

    if length_A1==physics.Np:
        a1 = A1
    elif length_A1==1:
        a1 = A1*_sp.ones(physics.Np)        
    elif length_A1>=phase.Np:
        a1 = A1[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A1!') 

    if length_A2==physics.Np:
        a2 = A2
    elif length_A2==1:
        a2 = A2*_sp.ones(physics.Np)        
    elif length_A2>=phase.Np:
        a2 = A2[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A2!') 

    if length_A3==physics.Np:
        a3 = A3
    elif length_A3==1:
        a3 = A3*_sp.ones(physics.Np)        
    elif length_A3>=phase.Np:
        a3 = A3[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A3!')  

    if length_A4==physics.Np:
        a4 = A4
    elif length_A4==1:
        a4 = A4*_sp.ones(physics.Np)        
    elif length_A4>=phase.Np:
        a4 = A4[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A4!') 
        
        
    if length_A5==physics.Np:
        a5 = A5
    elif length_A5==1:
        a5 = A5*_sp.ones(physics.Np)        
    elif length_A5>=phase.Np:
        a5 = A5[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A5!')         

    if length_A6==physics.Np:
        a6 = A6
    elif length_A6==1:
        a6 = A6*_sp.ones(physics.Np)        
    elif length_A6>=phase.Np:
        a6 = A6[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A6!') 

    S1 = a1*a3*a4*_sp.log(a2)*a2**(a3*X**a4+a5)*X**(a4-1)
    S2 = a1*a2**(a3*X**a4+a5)*(1-a3*a4*_sp.log(a2)*X**a4)+a6
    r = _sp.vstack((S1,S2)).T    
    return r


def natural_exponential(physics,
                        phase,
                        A1=0,
                        A2=0,
                        A3=0,
                        A4=0,
                        A5=0,
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
    A1 -> A5 : float or array
        The numerical values of the source term model's constant coefficients
    
    Notes
    -----
    It is possible to change the values of these constants at a later point.
    The values are stored in the Physics object's models dictionary under the
    chosen ``propname``.  The individual constants can be acesses as:
    phys.models[propname][A1].
    
    """        
    
    if x is None:   
        X = _sp.zeros(physics.Np)
        length_X = _sp.size(x)
    else:
        length_X = _sp.size(x)
        if length_X==physics.Np:
            X = x
        elif length_X==1:
            X = x*_sp.ones(physics.Np)  
        elif length_X>=phase.Np:
            X = x[physics.map_pores()] 
        else:   raise Exception('Wrong size for the guess value!')     

    length_A1 = _sp.size(A1)
    length_A2 = _sp.size(A2)
    length_A3 = _sp.size(A3)
    length_A4 = _sp.size(A4)
    length_A5 = _sp.size(A5)

    if length_A1==physics.Np:
        a1 = A1
    elif length_A1==1:
        a1 = A1*_sp.ones(physics.Np)        
    elif length_A1>=phase.Np:
        a1 = A1[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A1!') 

    if length_A2==physics.Np:
        a2 = A2
    elif length_A2==1:
        a2 = A2*_sp.ones(physics.Np)        
    elif length_A2>=phase.Np:
        a2 = A2[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A2!') 

    if length_A3==physics.Np:
        a3 = A3
    elif length_A3==1:
        a3 = A3*_sp.ones(physics.Np)        
    elif length_A3>=phase.Np:
        a3 = A3[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A3!')  

    if length_A4==physics.Np:
        a4 = A4
    elif length_A4==1:
        a4 = A4*_sp.ones(physics.Np)        
    elif length_A4>=phase.Np:
        a4 = A4[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A4!') 
        
        
    if length_A5==physics.Np:
        a5 = A5
    elif length_A5==1:
        a5 = A5*_sp.ones(physics.Np)        
    elif length_A5>=phase.Np:
        a5 = A5[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A5!')         


    S1 = a1*a2*a3*X**(a3-1)*_sp.exp(a2*X**a3+a4)
    S2 = a1*_sp.exp(a2*X**a3+a4)*(1-a2*a3*X**a3)+a5
    r = _sp.vstack((S1,S2)).T    
    return r
    
    
def logarithm(physics,
              phase,
              A1=0,
              A2=0,
              A3=0,
              A4=0,
              A5=0,
              A6=0,
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
    A1 -> A6 : float or array
        The numerical values of the source term model's constant coefficients
    
    Notes
    -----
    It is possible to change the values of these constants at a later point.
    The values are stored in the Physics object's models dictionary under the
    chosen ``propname``.  The individual constants can be acesses as:
    phys.models[propname][A1].
    
    """     
    
    if x is None:   
        X = _sp.zeros(physics.Np)
        length_X = _sp.size(x)
    else:
        length_X = _sp.size(x)
        if length_X==physics.Np:
            X = x
        elif length_X==1:
            X = x*_sp.ones(physics.Np)  
        elif length_X>=phase.Np:
            X = x[physics.map_pores()] 
        else:   raise Exception('Wrong size for the guess value!')     

    length_A1 = _sp.size(A1)
    length_A2 = _sp.size(A2)
    length_A3 = _sp.size(A3)
    length_A4 = _sp.size(A4)
    length_A5 = _sp.size(A5)
    length_A6 = _sp.size(A6)

    if length_A1==physics.Np:
        a1 = A1
    elif length_A1==1:
        a1 = A1*_sp.ones(physics.Np)        
    elif length_A1>=phase.Np:
        a1 = A1[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A1!') 

    if length_A2==physics.Np:
        a2 = A2
    elif length_A2==1:
        a2 = A2*_sp.ones(physics.Np)        
    elif length_A2>=phase.Np:
        a2 = A2[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A2!') 

    if length_A3==physics.Np:
        a3 = A3
    elif length_A3==1:
        a3 = A3*_sp.ones(physics.Np)        
    elif length_A3>=phase.Np:
        a3 = A3[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A3!')  

    if length_A4==physics.Np:
        a4 = A4
    elif length_A4==1:
        a4 = A4*_sp.ones(physics.Np)        
    elif length_A4>=phase.Np:
        a4 = A4[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A4!') 
        
        
    if length_A5==physics.Np:
        a5 = A5
    elif length_A5==1:
        a5 = A5*_sp.ones(physics.Np)        
    elif length_A5>=phase.Np:
        a5 = A5[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A5!')         

    if length_A6==physics.Np:
        a6 = A6
    elif length_A6==1:
        a6 = A6*_sp.ones(physics.Np)        
    elif length_A6>=phase.Np:
        a6 = A6[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A6!') 

    S1 = a1*a3*a4*X**(a4-1)/(_sp.log(a2)*(a3*X**a4+a5))
    S2 = a1*_sp.log(a3*X**a4+a5)/_sp.log(a2)+a6-a1*a3*a4*X**a4/(_sp.log(a2)*(a3*X**a4+a5))
    r = _sp.vstack((S1,S2)).T    
    return r
    
def natural_logarithm(physics,
                      phase,
                      A1=0,
                      A2=0,
                      A3=0,
                      A4=0,
                      A5=0,
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
    A1 -> A5 : float or array
        The numerical values of the source term model's constant coefficients
    
    Notes
    -----
    It is possible to change the values of these constants at a later point.
    The values are stored in the Physics object's models dictionary under the
    chosen ``propname``.  The individual constants can be acesses as:
    phys.models[propname][A1].
    
    """         
    
    if x is None:   
        X = _sp.zeros(physics.Np)
        length_X = _sp.size(x)
    else:
        length_X = _sp.size(x)
        if length_X==physics.Np:
            X = x
        elif length_X==1:
            X = x*_sp.ones(physics.Np)  
        elif length_X>=phase.Np:
            X = x[physics.map_pores()] 
        else:   raise Exception('Wrong size for the guess value!')     

    length_A1 = _sp.size(A1)
    length_A2 = _sp.size(A2)
    length_A3 = _sp.size(A3)
    length_A4 = _sp.size(A4)
    length_A5 = _sp.size(A5)

    if length_A1==physics.Np:
        a1 = A1
    elif length_A1==1:
        a1 = A1*_sp.ones(physics.Np)        
    elif length_A1>=phase.Np:
        a1 = A1[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A1!') 

    if length_A2==physics.Np:
        a2 = A2
    elif length_A2==1:
        a2 = A2*_sp.ones(physics.Np)        
    elif length_A2>=phase.Np:
        a2 = A2[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A2!') 

    if length_A3==physics.Np:
        a3 = A3
    elif length_A3==1:
        a3 = A3*_sp.ones(physics.Np)        
    elif length_A3>=phase.Np:
        a3 = A3[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A3!')  

    if length_A4==physics.Np:
        a4 = A4
    elif length_A4==1:
        a4 = A4*_sp.ones(physics.Np)        
    elif length_A4>=phase.Np:
        a4 = A4[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A4!') 
        
        
    if length_A5==physics.Np:
        a5 = A5
    elif length_A5==1:
        a5 = A5*_sp.ones(physics.Np)        
    elif length_A5>=phase.Np:
        a5 = A5[physics.map_pores()]   
    else:
        raise Exception('Wrong size for the parameter A5!')         


    S1 = a1*a2*a3*X**(a3-1)/(a2*X**a3+a4)
    S2 = a1*_sp.log(a2*X**a3+a4)+a5-a1*a2*a3*X**a3/(a2*X**a3+a4)
    r = _sp.vstack((S1,S2)).T    
    return r