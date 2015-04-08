r"""
===============================================================================
Submodule -- generic_source_term
===============================================================================

"""

import scipy as _sp

def linear(physics,
           phase,
           A1='',
           A2='',
           x='',
           return_rate=True,
           **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x  +  A_{2} 
    If return_rate is True, it returns the value of source term for the provided x in each pore. 
    If return_rate is False, it calculates the slope and intercept for the following linear form :                 
        .. math::            
            r = S_{1}   x  +  S_{2}               

    Parameters
    ----------
    A1 , A2 : string
        The property name of the coefficients in the source term model. With A2 set to zero
        this equation takes on the familiar for of r=kx.  
    x : string or float/int or array/list
        The property name or numerical value or array for the main quantity
    Notes
    -----    
    Because this source term is linear in concentration (x) is it not necessary
    to iterate during the solver step.  Thus, when using the 
    ``set_source_term`` method for an algorithm, it is recommended to set the ``maxiter`` 
    argument to 0.  This will save 1 unncessary solution of the system, since
    the solution would coverge after the first pass anyway.  
    
    """
    if x=='':   
        X = _sp.ones(physics.Np)*_sp.nan
    else:
        if type(x)==str:
            x = 'pore.'+x.split('.')[-1]
            try:    X = physics[x]
            except: raise Exception(physics.name+' does not have the pore property :'+x+'!')
        else:
            X = _sp.array(x)      

    length_X = _sp.size(X)
    if length_X!=physics.Np:
        if length_X==1:
            X = X*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = X[physics.map_pores()] 
        else:   raise Exception('Wrong size for the numerical array of x!')     

    if A1 == '':    a1 = 0
    else:
        if type(A1)==str:
            A1 = 'pore.'+A1.split('.')[-1]
            try:    a1 = physics[A1]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A1+'!')
        else:  raise Exception('A1 can be only string!')

    if A2 == '':    a2 = 0
    else:
        if type(A2)==str:
            A2 = 'pore.'+A2.split('.')[-1]
            try:    a2 = physics[A2]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A2+'!')
        else:  raise Exception('A2 can be only string!')

    if return_rate:
        return(a1*X**a2)
    else:
        S1 = a1
        S2 = a2
        return(_sp.vstack((S1,S2)).T)   

def power_law(physics,
            phase,
            A1='',
            A2='',
            A3='',
            x='',
            return_rate=True,
            **kwargs):
    r"""
    For the following source term:
        .. math::
            r = A_{1}   x^{A_{2}}  +  A_{3}
    If return_rate is True, it returns the value of source term for the provided x in each pore. 
    If return_rate is False, it calculates the slope and intercept for the following linear form :                 
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
    if x=='':   
        X = _sp.ones(physics.Np)*_sp.nan
#        length_X = _sp.size(x)
    else:
        if type(x)==str:
            x = 'pore.'+x.split('.')[-1]
            try:    X = physics[x]
            except: raise Exception(physics.name+' does not have the pore property :'+x+'!')
        else:
            X = _sp.array(x)      

    length_X = _sp.size(X)
    if length_X!=physics.Np:
        if length_X==1:
            X = X*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = X[physics.map_pores()] 
        else:   raise Exception('Wrong size for the numerical array of x!')     

    if A1 == '':    a1 = 0
    else:
        if type(A1)==str:
            A1 = 'pore.'+A1.split('.')[-1]
            try:    a1 = physics[A1]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A1+'!')
        else:  raise Exception('A1 can be only string!')

    if A2 == '':    a2 = 0
    else:
        if type(A2)==str:
            A2 = 'pore.'+A2.split('.')[-1]
            try:    a2 = physics[A2]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A2+'!')
        else:  raise Exception('A2 can be only string!')

    if A3 == '':    a3 = 0
    else:
        if type(A3)==str:
            A3 = 'pore.'+A3.split('.')[-1]
            try:    a3 = physics[A3]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A3+'!')
        else:  raise Exception('A3 can be only string!')

    if return_rate:
        return(a1*X**a2 +a3)
    else:
        S1 = a1*a2*X**(a2-1)
        S2 = a1*X**a2*(1-a2)+a3
        return(_sp.vstack((S1,S2)).T)    


def exponential(physics,
            phase,
            A1='',
            A2='',
            A3='',
            A4='',
            A5='',
            A6='',
            x='',
            return_rate=True,
            **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1} A_{2}^{( A_{3} x^{ A_{4} } + A_{5})} + A_{6} 
    If return_rate is True, it returns the value of source term for the provided x in each pore. 
    If return_rate is False, it calculates the slope and intercept for the following linear form :                 
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
    if x=='':   
        X = _sp.ones(physics.Np)*_sp.nan
#        length_X = _sp.size(x)
    else:
        if type(x)==str:
            x = 'pore.'+x.split('.')[-1]
            try:    X = physics[x]
            except: raise Exception(physics.name+' does not have the pore property :'+x+'!')
        else:
            X = _sp.array(x)      

    length_X = _sp.size(X)
    if length_X!=physics.Np:
        if length_X==1:
            X = X*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = X[physics.map_pores()] 
        else:   raise Exception('Wrong size for the numerical array of x!')     

    if A1 == '':    a1 = 1
    else:
        if type(A1)==str:
            A1 = 'pore.'+A1.split('.')[-1]
            try:    a1 = physics[A1]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A1+'!')
        else:  raise Exception('A1 can be only string!')

    if A2 == '':    a2 = 0
    else:
        if type(A2)==str:
            A2 = 'pore.'+A2.split('.')[-1]
            try:    a2 = physics[A2]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A2+'!')
        else:  raise Exception('A2 can be only string!')

    if A3 == '':    a3 = 0
    else:
        if type(A3)==str:
            A3 = 'pore.'+A3.split('.')[-1]
            try:    a3 = physics[A3]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A3+'!')
        else:  raise Exception('A3 can be only string!')


    if A4 == '':    a4 = 0
    else:
        if type(A4)==str:
            A4 = 'pore.'+A4.split('.')[-1]
            try:    a4 = physics[A4]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A4+'!')
        else:  raise Exception('A4 can be only string!')

    if A5 == '':    a5 = 0
    else:
        if type(A5)==str:
            A5 = 'pore.'+A5.split('.')[-1]
            try:    a5 = physics[A5]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A5+'!')
        else:  raise Exception('A5 can be only string!')

    if A6 == '':    a6 = 0
    else:
        if type(A6)==str:
            A6 = 'pore.'+A6.split('.')[-1]
            try:    a6 = physics[A6]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A6+'!')
        else:  raise Exception('A6 can be only string!')

    if return_rate:
        return(a1*a2**(a3*X**a4 +a5)+a6)
    else:
        S1 = a1*a3*a4*_sp.log(a2)*a2**(a3*X**a4+a5)*X**(a4-1)
        S2 = a1*a2**(a3*X**a4+a5)*(1-a3*a4*_sp.log(a2)*X**a4)+a6
        return(_sp.vstack((S1,S2)).T) 


def natural_exponential(physics,
                        phase,
                        A1='',
                        A2='',
                        A3='',
                        A4='',
                        A5='',
                        x='',
                        return_rate=True,
                        **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1} exp( A_{2}  x^{ A_{3} } + A_{4} )+ A_{5} 
    If return_rate is True, it returns the value of source term for the provided x in each pore. 
    If return_rate is False, it calculates the slope and intercept for the following linear form :                 
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
    if x=='':   
        X = _sp.ones(physics.Np)*_sp.nan
#        length_X = _sp.size(x)
    else:
        if type(x)==str:
            x = 'pore.'+x.split('.')[-1]
            try:    X = physics[x]
            except: raise Exception(physics.name+' does not have the pore property :'+x+'!')
        else:
            X = _sp.array(x)      

    length_X = _sp.size(X)
    if length_X!=physics.Np:
        if length_X==1:
            X = X*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = X[physics.map_pores()] 
        else:   raise Exception('Wrong size for the numerical array of x!')     

    if A1 == '':    a1 = 1
    else:
        if type(A1)==str:
            A1 = 'pore.'+A1.split('.')[-1]
            try:    a1 = physics[A1]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A1+'!')
        else:  raise Exception('A1 can be only string!')

    if A2 == '':    a2 = 0
    else:
        if type(A2)==str:
            A2 = 'pore.'+A2.split('.')[-1]
            try:    a2 = physics[A2]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A2+'!')
        else:  raise Exception('A2 can be only string!')

    if A3 == '':    a3 = 0
    else:
        if type(A3)==str:
            A3 = 'pore.'+A3.split('.')[-1]
            try:    a3 = physics[A3]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A3+'!')
        else:  raise Exception('A3 can be only string!')


    if A4 == '':    a4 = 0
    else:
        if type(A4)==str:
            A4 = 'pore.'+A4.split('.')[-1]
            try:    a4 = physics[A4]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A4+'!')
        else:  raise Exception('A4 can be only string!')

    if A5 == '':    a5 = 0
    else:
        if type(A5)==str:
            A5 = 'pore.'+A5.split('.')[-1]
            try:    a5 = physics[A5]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A5+'!')
        else:  raise Exception('A5 can be only string!')

    if return_rate:
        return(a1*_sp.exp(a2*X**a3+a4) +a5)
    else:
        S1 = a1*a2*a3*X**(a3-1)*_sp.exp(a2*X**a3+a4)
        S2 = a1*_sp.exp(a2*X**a3+a4)*(1-a2*a3*X**a3)+a5
        return(_sp.vstack((S1,S2)).T) 

def logarithm(physics,
            phase,
            A1='',
            A2='',
            A3='',
            A4='',
            A5='',
            A6='',
            x='',
            return_rate=True,
            **kwargs):
    r"""
    For the following source term:
        .. math::
            r =  A_{1}   Log_{ A_{2} }( A_{3} x^{ A_{4} }+ A_{5})+ A_{6}  
    If return_rate is True, it returns the value of source term for the provided x in each pore. 
    If return_rate is False, it calculates the slope and intercept for the following linear form :                 
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
    if x=='':   
        X = _sp.ones(physics.Np)*_sp.nan
#        length_X = _sp.size(x)
    else:
        if type(x)==str:
            x = 'pore.'+x.split('.')[-1]
            try:    X = physics[x]
            except: raise Exception(physics.name+' does not have the pore property :'+x+'!')
        else:
            X = _sp.array(x)      

    length_X = _sp.size(X)
    if length_X!=physics.Np:
        if length_X==1:
            X = X*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = X[physics.map_pores()] 
        else:   raise Exception('Wrong size for the numerical array of x!')     

    if A1 == '':    a1 = 1
    else:
        if type(A1)==str:
            A1 = 'pore.'+A1.split('.')[-1]
            try:    a1 = physics[A1]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A1+'!')
        else:  raise Exception('A1 can be only string!')

    if A2 == '':    a2 = 0
    else:
        if type(A2)==str:
            A2 = 'pore.'+A2.split('.')[-1]
            try:    a2 = physics[A2]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A2+'!')
        else:  raise Exception('A2 can be only string!')

    if A3 == '':    a3 = 0
    else:
        if type(A3)==str:
            A3 = 'pore.'+A3.split('.')[-1]
            try:    a3 = physics[A3]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A3+'!')
        else:  raise Exception('A3 can be only string!')


    if A4 == '':    a4 = 0
    else:
        if type(A4)==str:
            A4 = 'pore.'+A4.split('.')[-1]
            try:    a4 = physics[A4]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A4+'!')
        else:  raise Exception('A4 can be only string!')

    if A5 == '':    a5 = 0
    else:
        if type(A5)==str:
            A5 = 'pore.'+A5.split('.')[-1]
            try:    a5 = physics[A5]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A5+'!')
        else:  raise Exception('A5 can be only string!')

    if A6 == '':    a6 = 0
    else:
        if type(A6)==str:
            A6 = 'pore.'+A6.split('.')[-1]
            try:    a6 = physics[A6]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A6+'!')
        else:  raise Exception('A6 can be only string!')

    if return_rate:
        return(a1*_sp.log(a3*X**a4 +a5)/_sp.log(a2)+a6)
    else:
        S1 = a1*a3*a4*X**(a4-1)/(_sp.log(a2)*(a3*X**a4+a5))
        S2 = a1*_sp.log(a3*X**a4+a5)/_sp.log(a2)+a6-a1*a3*a4*X**a4/(_sp.log(a2)*(a3*X**a4+a5))
        return(_sp.vstack((S1,S2)).T) 
        

def natural_logarithm(physics,
                      phase,
                      A1='',
                      A2='',
                      A3='',
                      A4='',
                      A5='',
                      x='',
                      return_rate=True,
                      **kwargs):
    r"""
    For the following source term:
        .. math::
            r =   A_{1}  Ln( A_{2} x^{ A_{3} }+ A_{4})+ A_{5}    
    If return_rate is True, it returns the value of source term for the provided x in each pore. 
    If return_rate is False, it calculates the slope and intercept for the following linear form :                 
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
    if x=='':   
        X = _sp.ones(physics.Np)*_sp.nan
#        length_X = _sp.size(x)
    else:
        if type(x)==str:
            x = 'pore.'+x.split('.')[-1]
            try:    X = physics[x]
            except: raise Exception(physics.name+' does not have the pore property :'+x+'!')
        else:
            X = _sp.array(x)      

    length_X = _sp.size(X)
    if length_X!=physics.Np:
        if length_X==1:
            X = X*_sp.ones(physics.Np) 
        elif length_X>=phase.Np:
            X = X[physics.map_pores()] 
        else:   raise Exception('Wrong size for the numerical array of x!')     

    if A1 == '':    a1 = 1
    else:
        if type(A1)==str:
            A1 = 'pore.'+A1.split('.')[-1]
            try:    a1 = physics[A1]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A1+'!')
        else:  raise Exception('A1 can be only string!')

    if A2 == '':    a2 = 0
    else:
        if type(A2)==str:
            A2 = 'pore.'+A2.split('.')[-1]
            try:    a2 = physics[A2]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A2+'!')
        else:  raise Exception('A2 can be only string!')

    if A3 == '':    a3 = 0
    else:
        if type(A3)==str:
            A3 = 'pore.'+A3.split('.')[-1]
            try:    a3 = physics[A3]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A3+'!')
        else:  raise Exception('A3 can be only string!')


    if A4 == '':    a4 = 0
    else:
        if type(A4)==str:
            A4 = 'pore.'+A4.split('.')[-1]
            try:    a4 = physics[A4]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A4+'!')
        else:  raise Exception('A4 can be only string!')

    if A5 == '':    a5 = 0
    else:
        if type(A5)==str:
            A5 = 'pore.'+A5.split('.')[-1]
            try:    a5 = physics[A5]
            except: raise Exception(physics.name+'/'+phase.name+' does not have the pore property :'+A5+'!')
        else:  raise Exception('A5 can be only string!')


    if return_rate:
        return(a1*_sp.log(a2*X**a3 +a4)+a5)
    else:
        S1 = a1*a2*a3*X**(a3-1)/(a2*X**a3+a4)
        S2 = a1*_sp.log(a2*X**a3+a4)+a5-a1*a2*a3*X**a3/(a2*X**a3+a4)
        return(_sp.vstack((S1,S2)).T)         