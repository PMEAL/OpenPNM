r"""
===============================================================================
pore_seeds -- 
===============================================================================

"""
import scipy as _sp

def perlin_noise(divisions,freq=1,octaves=4,mode='classic',**kwargs):
    r'''
    Generate pore seed values using the Perlin noise algorithm.  This approach
    imparts some spatial clumpiness to the pore seeds.
    
    Parameters
    ----------
    divisions : list or tuple
        The x, y and z dimensions of the cubic network
    freq, octaves : int
        Parameters that control the qualities of the noise.  Lower frequency
        results in more smaller clumps.  Higher octaves gives a more textured
        noise.
    mode : {'classic','simplex'}
        The algorithm to use when generating the noise.
        
        * 'classic' : (default) The original algorithm developed by Perlin 
        * 'simplex' : A newer algorithm that is supposedly faster and results
        in a more natural texture.
    '''
    from noise import pnoise3, snoise3
    import scipy.stats as spst
    if mode == 'classic':
        model = pnoise3
    elif mode == 'simplex':
        model = snoise3
    freq = freq * octaves
    x = divisions[0]
    y = divisions[1]
    z = divisions[2]
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
    return values

def distance_from_solid(geometry,array,**kwargs):
    r'''
    '''
    import scipy.ndimage as _spim
    #Pad image by reflecting
#    array = _sp.tile(array,(3,3,3))
    #Perform distance transform
    img = _spim.distance_transform_bf(array)
#    img = img[25:50,25:50,25:50]
    #Convert back to pore-list
    seeds = _sp.reshape(img,(geometry.Np,),order='F')
    return seeds
    
    
    
















