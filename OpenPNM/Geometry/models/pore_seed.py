r"""
===============================================================================
pore_seeds -- Methods for generating fields of values for use as seeds in
statistical pore size distributions
===============================================================================

"""
import scipy as _sp

def perlin_noise(geometry,freq=1,octaves=4,mode='classic',**kwargs):
    r'''
    Generate pore seed values using the Perlin noise algorithm.  This approach
    imparts some spatial clumpiness to the pore seeds.
    
    Parameters
    ----------
    freq, octaves : int
        Parameters that control the qualities of the noise.  Lower frequency
        results in more smaller clumps.  Higher octaves gives a more textured
        noise.
    mode : {'classic','simplex'}
        The algorithm to use when generating the noise.
        
        * 'classic' : (default) The original algorithm developed by Perlin 
        * 'simplex' : A newer algorithm that is supposedly faster and results
        in a more natural texture.
        
    Returns
    -------
    A list of pore seeds values between 0:1 that are spatially correlated (i.e.
    similar values are clumped together)
        
    Notes
    -----
    - This method uses image analysis type tools, so only works on Cubic 
    networks
    - This method requires the 'noise' module is installed
    
    '''
    from noise import pnoise3, snoise3
    import scipy.stats as spst
    
    net = geometry._net
    if mode == 'classic':
        model = pnoise3
    elif mode == 'simplex':
        model = snoise3
    freq = freq * octaves
    #The following will only work on Cubic networks
    x = net._shape[0]
    y = net._shape[1]
    z = net._shape[2]
    temp = _sp.ndarray((x,y,z))
    for k in range(z):
        for j in range(y):
            for i in range(x):
                temp[i,j,k] = model(i / freq, j / freq, k / freq, octaves) + 0.5
    #Assuming that the noise is normally distributed, find seeds of that dist
    temp = _sp.reshape(temp,(temp.size,))
    x_mean = _sp.mean(temp)
    x_sigma = _sp.sqrt(1/(temp.size-1)*_sp.sum((temp - x_mean)**2))
    fn1 = spst.norm(loc=x_mean,scale=x_sigma)
    values = fn1.cdf(temp)
    values = values[geometry['pore.map']]
    return values.flatten()

def distance_from_inclusion(geometry,p,**kwargs):
    r'''
    Genrate spatially correlated pore seeds by calculating distance from random 
    locations (inclusions) in the domain
    
    Parameters
    ----------
    p : float
        The fraction of pores in the domain that are set as 'seeds' for the 
        distance calculation
    
    Returns
    -------
    A list of distance values (in voxels) between each pore and it nearest
    seed pore.  A list of voxel distances is returned rather than normalized 
    seeds between 0:1 so that the user can manipulate the map as desired, by 
    applying desired thresholds and/or scaling to get 0:1 seeds.
    
    Notes
    -----
    This method uses image analysis type tools, so only works on Cubic networks
    '''
    import scipy.ndimage as _spim
    net = geometry._net
    #The following will only work on Cubic networks
    x = net._shape[0]
    y = net._shape[1]
    z = net._shape[2]
    img = _sp.rand(x,y,z)>p
    #Pad image by tiling
    a = _sp.tile(img,[3,3,3])
    b = a[x:-x,y:-y,z:-z]
    #Perform distance transform
    img = _spim.distance_transform_bf(b)
    #Convert back to pore-list
    values = img.flatten()
    values = values[geometry['pore.map']]
    return values
    
    