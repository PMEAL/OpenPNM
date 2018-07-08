.. _transport_guide:

================================================================================
Transport Simulations
================================================================================

.. contents:: Page Contents
    :depth: 3

The Transport algorithms are broken into 3 separate classes, each inheriting from the one above it.  The top level class is :ref:`generic_transport_api` which implements the steady-state transport solver. Next is :ref:`reactive_transport_api` which adds the ability to apply source terms.  Finally, the bottom class is :ref:`transient_reactive_transport_api` which manages the transient solution of the transport problems.  Each of these generic classes is then subclassed for specific phyical problems such as *StokesFlow* or *TransientFickianDiffusion*.

Consider the following example:

.. code-block:: python

  >>> import openpnm as op
  >>> pn = op.network.Cubic(shape=[5, 5, 5], spacing=1e-4)
  >>> phase = op.phases.GenericPhase(network=pn)
  >>> phase['throat.conductance'] = 1.0  # Use dummy values

--------------------------------------------------------------------------------
Steady-State Transport
--------------------------------------------------------------------------------

The steady-state transport class contains all of the methods for setting up and solving the linear system of equations describing the material balance around each pore to find the unknown quantity *x*:

.. math::

  x = A^{-1} b

The ``GenericTransport`` class handles the following steps:

* Building the coefficient matrix *A*
* Building the right-hand side matrix *b*
* Specifying boundary values and locations
* Applying the boundary conditions by updating *A* and *b*
* Calling the specified solver to find values of the unknown quantity *x*

Continuing with the example from above:

.. code-block:: python

  >>> alg = op.algorithms.GenericTransport(network=pn)
  >>> alg.setup(phase=phase, quantity='pore.quantity', conductance='throat.conductance')

................................................................................
Building the Coefficient and Right-Hand Size Matrices
................................................................................

The material balance around a given pore *i* assuming steady-state and no reaction results in the following:

.. math::

    \sum_{j=1}^n {g_{i,j} \cdot (x_i - x_j)} = 0

where *j* is the index of the neighboring pores, of which there are a total of *n*, *g* is the conductance between pores *i* and *j*, and *x* is the unknown quantity being solved for.  Expanding this equation into matrix form, for i = 1 and j = [2, 3, 4] yields:

.. math::

  -( g_{1,2} + g_{1,3} + g_{1,4} ) \cdot x_1 + g_{1,2} \cdot x_2 + g_{1,3} \cdot x_3 + g_{1,4} \cdot x_4 = 0

The *A* matrix is composed of a set of linear equation analogous to this one for each each pore in the network, so a network with 1000 pores will have an A matrix of 1000 by 1000.  This could very quickly grow unreasonable in memory requirements, so OpenPNM always uses `sparse matrices <https://en.wikipedia.org/wiki/Sparse_matrix>`_ from the `scipy.sparse module <https://docs.scipy.org/doc/scipy/reference/sparse.html>`_.  The sparse *A* matrix is constructed using the ``_build_A`` method of the ``GenericTransport`` class automatically and is accessible via the ``A`` attribute.  The *A* matrix can be inspected as follows:

.. code-block:: python

  >>> alg.A.shape, alg.A.format, alg.A.nnz
  ((125, 125), 'coo', 725)

For steady-state transport with no source or sink terms the right-hand side matrix, *b* is all 0's (*b* is not stored as a sparse matrix since is just single column).

................................................................................
Specifying and Applying Boundary Conditions
................................................................................

OpenPNM supports the two most common type of boundary conditions: constant values (also known as Dirichlet) and constant rate (similar to Neumann).  The constant value condition is equivalent to writing:

.. math::

  x_i = b_i

This is applied to the *A* and *b* matrices by removing the balance equation on row *i* and replacing it with this.  The constant rate boundary condition is even easier to implement since it is simply a matter of setting the right-hand side of each balance equation to some non-zero value.

The operation of altering *A* and *b* accordingly is performed by the ``_apply_BCs`` method, which occurs automatically when the user calls ``run``. In terms of the *A* and *b* matrices, constant rate BCs require setting :math:`b_i` to a value of :math:`n_i`, and *A* is untouched, while constant value BCs require setting :math:`b_i` to 1 while replacing all elments in row *i* with 0, and setting the diagonal to :math:`x_i`.

Boundary conditions are specified by the user with the ``set_value_BC`` and/or ``set_rate_BC`` methods.  Each of these methods requires the numerical value of the BC as well as which pores to apply them (BCs can only be applied in pores at this time).  These methods store the values they recieve on the Algorithm object under 'pore.value_BC' and 'pore.rate_BC', with the given values at the corresponding locations, and NaNs elsewhere.

.. code-block:: python

  >>> alg.set_value_BC(pores=pn.pores('left'), values=1)
  >>> alg.set_rate_BC(pores=pn.pores('right'), values=1)
  >>> print(alg)
  ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
  openpnm.algorithms.GenericTransport : alg_01
  ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
  #     Properties                                    Valid Values
  ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
  1     pore.bc_rate                                     25 / 125
  2     pore.bc_value                                    25 / 125
  ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
  #     Labels                                        Assigned Locations
  ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――
  1     pore.all                                      125
  2     throat.all                                    300
  ――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――――

................................................................................
Choosing a Solver and Running the Simulation
................................................................................

The final step is solving the system of equations is simply calling the ``run`` method.  Before doing so, however, you can specify which solver to use in ``settings['solver']``.  The default is ``'spsolve'`` which in turn uses the default Scipy sparse solver.  On a vanilla install of Scipy, this will likely be SuperLU, which is very stable but slow.  If the ``scikit-umfpack`` package has been installed, then Scipy will automatically use this by default, which is much faster.  It is also possible to specify any of the `iterative solvers offered by Scipy <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems>`_.  For instance, to use conjugate gradient, use ``settings['solver'] = 'cg'``.  Iterative solvers a much faster and can handle larger systems, but they are not always stable, so much be used with care.  Most OpenPNM problems are well handled by the conjugate gradient solver.

The solution is produced by calling ``run`` method, which actually calls a few other methods behind the scenes.  For the sake of illustration, let's call these explicitly:

.. code-block:: python

  >>> alg._apply_BCs()
  >>> x = alg._solve(A=alg.A, b=alg.b)

The ``_solve`` method computes *x* and returns it.  The ``run`` method, which essentially just calls the above 2 methods, captures the received value of *x* and stores it on the Algorithm under ``'pore.quantity'``.  The name of the quantity is specified in the ``settings['quantity']`` and is given sensible names by default for the various subclasses (e.g.for  StokesFlow it is 'pore.pressure').

--------------------------------------------------------------------------------
Reactive Transport
--------------------------------------------------------------------------------

The ReactiveTransport class inherits directly from GenericTransport, so inherits all of the mechanism described above, plus the ability to include non-linear source and sink terms.  The balance equation around a pore *i* in the presence of a source/sink term is:

.. math::

    \sum_{j=1}^n {g_{i,j} \cdot (x_i - x_j)} + R(x_i) = 0

where :math:`R(x_i)` is some function of the quantity being solved for, and can be non-linear.  A common example is the standard 2nd order kinetics: :math:`R = A \cdot x^2`.

To have access to the source/sink machinery we must use an instance of ReactiveTransport:

.. code-block:: python

  >>> alg = op.algorithms.ReactiveTransport(network=pn)
  >>> alg.setup(phase=phase, quantity='pore.quantity', conductance='throat.conductance')
  >>> alg.set_value_BC(pores=pn.pores('left'), values=1)

................................................................................
Setting Up and Specifying Source Terms
................................................................................

Specifying a source/sink term requires first defining the form of the equation and its constants.  This is done on the Physics or Phase objects as a pore-scale model:

.. code-block:: python

  >>> phase.add_model(propname='pore.rxn',
  ...                 model=op.models.physics.generic_source_term.power_law,
  ...                 A1='pore.A', A2='pore.n', A3='pore.A3',
  ...                 X='pore.quantity')
  >>> phase['pore.A'] = -1.0
  >>> phase['pore.n'] = 2
  >>> phase['pore.A3'] = 0

Now the Algorithm can be told where to look for the source term, and where to apply it:

.. code-block:: python

  >>> alg.set_source(propname='pore.rxn', pores=63)  # A single pore near the middle

The process of setting a source/sink term does two things.  It places 'pore.rxn' in ``alg.settings['sources']`` and it creates a label called 'pore.rxn' indicating which pores it applied to.

.. note:: All Transport Classes Have Reactions

  If no source terms are specified to the algorithm then no attempt is made by the algorithm to add source/sink terms to the matrices, and no iterations are performed.  In other words, when no source/sink terms are specified (using ``set_source``) the ReactiveTransport class behaves *exactly* the same at the GenericTransport.  Hence, all transport Algorithms are subclasses of ReactiveTransport.

................................................................................
Adjusting the Coefficent and Right-Hand Side Matrix
................................................................................

If the source/sink term were linear (e.g. :math:`R_i = k \cdot x_i)`, then simply adding :math:`-k_i` to the diagonal of the *ith* row of *A* would be sufficient.  However, to handle the general case of non-linear source terms, OpenPNM uses a method based on Newton's method adapted from `Numerical Heat Transfer and Fluid Flow by Patankar <https://www.crcpress.com/Numerical-Heat-Transfer-and-Fluid-Flow/Patankar/p/book/9780891165224>`_.  This involves linearizing the source term about the current value of *x* such that :math:`R_i = S_1 \cdot x_i + S_2`, which means that :math:`S_1` is added the diagonal of *A* and :math:`-S_2` is added to *b*.

The ReactiveTransport class has a ``_apply_sources`` method which check ``alg.settings['sources']`` to see which source/sink terms have been added (if any).  For each sourc/sink term if finds, it regenerates that model in the associated Physics and/or Phase objects using the current value 'pore.quantity' (which defaults to 0).  Next the :math:`S_1` and :math:`S_2` terms are fetched from each Physics and/or Phase and applied to the corresponding pores by adding :math:`S_(1,i)` to the diagonal of row *i* of *A*, and :math:`_S_(2,i)` to row *i* of *b*.

................................................................................
Iteratively Solving for Non-Linear Source/Sink Terms
................................................................................

The ReactiveTransport has a ``run`` method that will apply all the necessary steps for the user, including ``_apply_sources``, and most importantly calling ``_solve`` of the GenericTransport class repeatedly until convergance is acheived on the quantity *x*.  In this case, convergence means that the value of *x* used when regenerating the source/sink terms to determine :math:`S_1` and :math:`S_2`, is sufficiently close the value of *x* returned by ``_solve``.

To run the simulation using explicit steps:

.. code-block:: python

  >>> alg._build_A()
  >>> alg._build_b()
  >>> alg._apply_BCs()
  >>> alg['pore.quantity'] = 0  # Make initial guess of quantity
  >>> alg._apply_sources()
  >>> x = alg._run_reactive(x=None)

--------------------------------------------------------------------------------
Transient and Reactive Transport
--------------------------------------------------------------------------------

Transient algorithms inherit from ReactiveTransport, so possess all the machinery described above, plus some extra methods for setting up and performing the simulation transiently.  The only additional method is ``set_IC`` for setting the initial conditions, plus there are a number of extra ``settings`` required, specifically, the start time, end time, and time step (``'t_start', 't_stop', 't_step'``).

OpenPNM offers two transient transport solvers: The implicit scheme is used by default, but the Crank-Nicolson scheme can be used by setting ``self.settings['t_scheme'] = 'cranknicolson'``.
