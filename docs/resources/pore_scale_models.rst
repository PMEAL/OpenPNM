.. _pore_scale_models:

###############################################################################
Spatially correlated pore seeds using Perlin noise
###############################################################################
Perlin noise [ref?] is an algorithm for creating spatially correlated random
noise.  

.. code-block:: python

    def perlin_noise(geometry, freq=1, octaves=4, mode='classic', **kwargs):
        r"""
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

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.Cubic(shape=[50, 50, 50])
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.Ps,throats=pn.Ts)
        >>> geom.add_model(propname='pore.seed',
        ...                model=OpenPNM.Geometry.models.pore_seed.perlin_noise)
        >>> im = pn.asarray(geom['pore.seed'])

        Visualizing the end result can be done with:

        .. code-block:: python

            matplotlib.pyplot.imshow(im[:, 25, :],interpolation='none')

        """
        from noise import pnoise3, snoise3
        import scipy.stats as spst

        net = geometry._net
        if mode == 'classic':
            model = pnoise3
        elif mode == 'simplex':
            model = snoise3
        freq = freq * octaves
        # The following will only work on Cubic networks
        x = net._shape[0]
        y = net._shape[1]
        z = net._shape[2]
        temp = _sp.ndarray((x, y, z))
        for k in range(z):
            for j in range(y):
                for i in range(x):
                    temp[i, j, k] = model(i / freq, j / freq, k / freq, octaves) + 0.5
        # Assuming that the noise is normally distributed, find seeds of that dist
        temp = _sp.reshape(temp, (temp.size,))
        x_mean = _sp.mean(temp)
        x_sigma = _sp.sqrt(1/(temp.size-1)*_sp.sum((temp - x_mean)**2))
        fn1 = spst.norm(loc=x_mean, scale=x_sigma)
        values = fn1.cdf(temp)
        values = values[geometry.map_pores(target=net, pores=geometry.Ps)]
        return values.flatten()