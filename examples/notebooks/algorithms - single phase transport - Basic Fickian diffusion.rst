Fickian Diffusion
=================

One of the main applications of ``OpenPNM`` is simulating transport
phenomena such as Fickian diffusion, advection diffusion, reactive
transport, etc. In this example, we will learn how to perform Fickian
diffusion on a ``Cubic`` network. The algorithm works fine with every
other network type, but for now we want to keep it simple.

.. code:: ipython3

    import numpy as np
    import openpnm as op
    %matplotlib inline
    np.random.seed(10)
    ws = op.Workspace()
    ws.settings["loglevel"] = 40
    np.set_printoptions(precision=5)

Generating network
------------------

First, we need to generate a ``Cubic`` network. For now, we stick to a
2d network, but you might as well try it in 3d!

.. code:: ipython3

    net = op.network.Cubic(shape=[1, 10, 10], spacing=1e-5)

Adding geometry
---------------

Next, we need to add a geometry to the generated network. A geometry
contains information about size of the pores/throats in a network.
``OpenPNM`` has tons of prebuilt geometries that represent the
microstructure of different materials such as Toray090 carbon papers,
sand stone, electrospun fibers, etc. For now, we stick to a sample
geometry called ``StickAndBall`` that assigns random values to
pore/throat diameters.

.. code:: ipython3

    geom = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

Adding phase
------------

Next, we need to add a phase to our simulation. A phase object(s)
contain(s) thermophysical information about the working fluid(s) in the
simulation. ``OpenPNM`` has tons of prebuilt phases as well! For this
simulation, we use air as our working fluid.

.. code:: ipython3

    air = op.phases.Air(network=net)

Adding physics
--------------

Finally, we need to add a physics. A physics object contains information
about the working fluid in the simulation that depend on the geometry of
the network. A good example is diffusive conductance, which not only
depends on the thermophysical properties of the working fluid, but also
depends on the geometry of pores/throats.

.. code:: ipython3

    phys_air = op.physics.Standard(network=net, phase=air, geometry=geom)

Performing Fickian diffusion
============================

Now that everything’s set up, it’s time to perform our Fickian diffusion
simulation. For this purpose, we need to add the ``FickianDiffusion``
algorithm to our simulation. Here’s how we do it:

.. code:: ipython3

    fd = op.algorithms.FickianDiffusion(network=net, phase=air)

Note that ``network`` and ``phase`` are required parameters for pretty
much every algorithm we add, since we need to specify on which network
and for which phase do we want to run the algorithm.

Adding boundary conditions
--------------------------

Next, we need to add some boundary conditions to the simulation. By
default, ``OpenPNM`` assumes zero flux for the boundary pores.

.. code:: ipython3

    inlet  = net.pores('left') 
    outlet = net.pores('right')
    fd.set_value_BC(pores=inlet, values=1.0)
    fd.set_value_BC(pores=outlet, values=0.0)

``set_value_BC`` applies the so-called “Dirichlet” boundary condition to
the specified pores. Note that unless you want to apply a single value
to all of the specified pores (like we just did), you must pass a list
(or ``ndarray``) as the ``values`` parameter.

Running the algorithm
---------------------

Now, it’s time to run the algorithm. This is done by calling the ``run``
method attached to the algorithm object.

.. code:: ipython3

    fd.run()

Post processing
===============

When an algorithm is successfully run, the results are attached to the
same object. To access the results, you need to know the quantity for
which the algorithm was solving. For instance, ``FickianDiffusion``
solves for the quantity ``pore.concentration``, which is somewhat
intuitive. However, if you ever forget it, or wanted to manually check
the quantity, you can take a look at the algorithm ``settings``:

.. code:: ipython3

    print(fd.settings)


.. parsed-literal::

    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    key                                 value
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    cache_A                             True
    cache_b                             True
    conductance                         throat.diffusive_conductance
    phase                               phase_01
    quantity                            pore.concentration
    solver_atol                         None
    solver_family                       scipy
    solver_max_iter                     5000
    solver_preconditioner               jacobi
    solver_rtol                         None
    solver_tol                          1e-08
    solver_type                         spsolve
    prefix                              alg
    nlin_max_iter                       5000
    relaxation_quantity                 1.0
    relaxation_source                   1.0
    sources                             []
    variable_props                      []
    ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
    

Now that we know the quantity for which ``FickianDiffusion`` was solved,
let’s take a look at the results:

.. code:: ipython3

    c = fd['pore.concentration']
    print(c)


.. parsed-literal::

    [1.      1.      1.      1.      1.      1.      1.      1.      1.
     1.      0.91337 0.9027  0.90122 0.912   0.91149 0.89046 0.87493 0.87743
     0.85069 0.84281 0.80997 0.80214 0.77949 0.82358 0.84042 0.82003 0.81628
     0.77634 0.80587 0.7963  0.71402 0.69612 0.70543 0.70796 0.70568 0.68883
     0.68759 0.66212 0.69278 0.73835 0.63602 0.63711 0.59846 0.6043  0.58074
     0.56718 0.57254 0.59279 0.55155 0.57938 0.50945 0.48409 0.44291 0.47565
     0.51677 0.51963 0.51152 0.49198 0.45323 0.42504 0.39419 0.35525 0.35866
     0.33035 0.37632 0.43706 0.37632 0.34388 0.32492 0.30221 0.24301 0.25058
     0.25166 0.25424 0.22324 0.25012 0.21119 0.2167  0.20701 0.20723 0.0852
     0.12852 0.1457  0.12419 0.11529 0.11315 0.107   0.11574 0.12471 0.10392
     0.      0.      0.      0.      0.      0.      0.      0.      0.
     0.     ]
    

Heatmap
-------

Let’s visualize the results. Since the network is 2d, we can simply
reshape the results in form of a 2d array similar to the shape of the
network and plot the heatmap of it using ``matplotlib``.

.. code:: ipython3

    print('Network shape:', net._shape)
    c2d = c.reshape((net._shape))


.. parsed-literal::

    Network shape: (1, 10, 10)
    

.. code:: ipython3

    #NBVAL_IGNORE_OUTPUT
    import matplotlib.pyplot as plt
    plt.imshow(c2d[0,:,:]);
    plt.title('Concentration (mol/m$^3$)')
    plt.colorbar()




.. parsed-literal::

    <matplotlib.colorbar.Colorbar at 0x29bf969e860>




.. image:: algorithms%20-%20single%20phase%20transport%20-%20Basic%20Fickian%20diffusion_files%5Calgorithms%20-%20single%20phase%20transport%20-%20Basic%20Fickian%20diffusion_25_1.png


Calculating heat flux
---------------------

You might as well be interested in calculating the mass flux from a
boundary! This is easily done in ``OpenPNM`` via calling the ``rate``
method attached to the algorithm. Let’s see how it works:

.. code:: ipython3

    rate_inlet = fd.rate(pores=inlet)[0]
    print(f'Mass flow rate from inlet: {rate_inlet:.5e} mol/s')


.. parsed-literal::

    Mass flow rate from inlet: 7.16784e-12 mol/s
    
