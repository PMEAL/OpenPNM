.. _getting_started:


***********************************
Getting started with sphinx and rst
***********************************

.. _rst-guide:

RST tutorials
=============

The text markup is based on rst. There are quite a few good exmaples
and tutorials on the web:

.. list-table:: rst tutorials.
   :widths: 10 90
   :header-rows: 1

   * - Website
     - Description
   * - `basic <http://people.ee.ethz.ch/~creller/web/tricks/reST.html>`__
     - A short tutorial covering the basics. Well organized
   * - `math <http://sphinx-doc.org/ext/math.html>`__
     - A short document on how to include mathematical formulas
   * - `inria <http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html>`__
     - A good introduction / example at INRIA

Docscrape module
----------------

Known sections of the docsrape module:

.. warning:: I have not figured this out yet. It only allows specific sections in the class documentation.

The doctest module
==================

.. testcode::

   1+1        # this will give no output!
   print(2+2)  # this will give output

.. testoutput::

   4



.. _ipython-highlighting:

ipython sessions
================

Michael Droettboom contributed a sphinx extension which does `pygments
<http://pygments.org>`_ syntax highlighting on `ipython
<http://ipython.scipy.org>`_ sessions.  Just use ipython as the
language in the ``sourcecode`` directive::

    .. sourcecode:: ipython

        In [69]: lines = plot([1,2,3])

        In [70]: setp(lines)
          alpha: float
          animated: [True | False]
          antialiased or aa: [True | False]
          ...snip


and you will get the syntax highlighted output below.

.. sourcecode:: ipython

    In [69]: lines = plot([1,2,3])

    In [70]: setp(lines)
      alpha: float
      animated: [True | False]
      antialiased or aa: [True | False]
      ...snip

This support is included in this template, but will also be included
in a future version of Pygments by default.

.. _using-math:

Using math
==========

In sphinx you can include inline math :math:`x\leftarrow y\ x\forall
y\ x-y` or display math

.. math::

  W^{3\beta}_{\delta_1 \rho_1 \sigma_2} = U^{3\beta}_{\delta_1 \rho_1} + \frac{1}{8 \pi 2} \int^{\alpha_2}_{\alpha_2} d \alpha^\prime_2 \left[\frac{ U^{2\beta}_{\delta_1 \rho_1} - \alpha^\prime_2U^{1\beta}_{\rho_1 \sigma_2} }{U^{0\beta}_{\rho_1 \sigma_2}}\right]

To include math in your document, just use the math directive; here is
a simpler equation::

    .. math::

      W^{3\beta}_{\delta_1 \rho_1 \sigma_2} \approx U^{3\beta}_{\delta_1 \rho_1}

which is rendered as

.. math::

   W^{3\beta}_{\delta_1 \rho_1 \sigma_2} \approx U^{3\beta}_{\delta_1 \rho_1}

This documentation framework includes a Sphinx extension,
:file:`sphinxext/mathmpl.py`, that uses matplotlib to render math
equations when generating HTML, and LaTeX itself when generating a
PDF.  This can be useful on systems that have matplotlib, but not
LaTeX, installed.  To use it, add ``mathmpl`` to the list of
extensions in :file:`conf.py`.

:math:`\mbox{\LaTeX}`

Current SVN versions of Sphinx now include built-in support for math.
There are two flavors:

  - pngmath: uses dvipng to render the equation

  - jsmath: renders the math in the browser using Javascript

To use these extensions instead, add ``sphinx.ext.pngmath`` or
``sphinx.ext.jsmath`` to the list of extensions in :file:`conf.py`.

All three of these options for math are designed to behave in the same
way.

See the matplotlib `mathtext guide
<http://matplotlib.sourceforge.net/users/mathtext.html>`_ for lots
more information on writing mathematical expressions in matplotlib.

You can also inline code for plots directly, and the code will be
executed at documentation build time and the figure inserted into your
docs; the following code::

   .. plot::

      import matplotlib.pyplot as plt
      import numpy as np
      x = np.random.randn(1000)
      plt.hist( x, 20)
      plt.grid()
      plt.title(r'Normal: $\mu=%.2f, \sigma=%.2f$'%(x.mean(), x.std()))
      plt.show()

produces this output:

.. plot::

    import matplotlib.pyplot as plt
    import numpy as np
    x = np.random.randn(1000)
    plt.hist( x, 20)
    plt.grid()
    plt.title(r'Normal: $\mu=%.2f, \sigma=%.2f$'%(x.mean(), x.std()))
    plt.show()


See the matplotlib `pyplot tutorial
<http://matplotlib.sourceforge.net/users/pyplot_tutorial.html>`_ and
the `gallery <http://matplotlib.sourceforge.net/gallery.html>`_ for
lots of examples of matplotlib plots.


Inheritance diagrams
====================

Inheritance diagrams can be inserted directly into the document by
providing a list of class or module names to the
``inheritance-diagram`` directive.

For example::

  .. inheritance-diagram:: codecs

produces:

.. inheritance-diagram:: codecs


.. _extensions-literal:

