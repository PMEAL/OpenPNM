.. highlight:: python
   :linenothreshold: 3


Developer Guide
===============

.. contents:: Table of Contents
   :local:
   :backlinks: top

This section purports to document the *SfePy* internals. It is mainly useful
for those who wish to contribute to the development of  *SfePy* and understand
the inner workings of the code.

OpenPNM Directory Structure
---------------------------

Here we list and describe the directories that are in the main sfepy
directory.

.. list-table:: Top directory structure.
   :widths: 10 90
   :header-rows: 1

   * - name
     - description
   * - `build/`
     - directory created by the build process (generated)
   * - `docs/`
     - source files of this documentation
   * - `examples/`
     - example problem description files
   * - `pnmdata/`
     - pore network data in various formats shared by the examples
   * - `output/`
     - default output directory for storing results of the examples
   * - `output-tests/`
     - output directory for tests
   * - `script/`
     - various small scripts (simple mesh generators, mesh format
       convertors etc.)
   * - `OpenPNM/`
     - the source code
   * - `tests/`
     - the tests run by `runTests.py`
   * - `tmp/`
     - directory for temporary files (generated)

New users/developers 
should explore the `examples/` directory. For developers, the principal
directory is `sfepy/`, which has the following contents:

.. list-table:: `OpenPNM/` directory structure.
   :widths: 10 80 10
   :header-rows: 1

   * - name
     - description
     - field-specific
   * - `misc`
     - common utilities and classes used by most of the othe rmodules
     -
   * - `TOP`
     - Storage and manipulations of network topoologies and data stored on them.
     -
   * - `Generators`
     - Generators for pore networks. (Random cubic, image based, Voronoi). Should also contain
       a mapper of the pore network back on the segmented image.
     -
   * - `Algorithms`
     - Module containing all algorithmic classes for networks.
     -
   * - `IO`
     - Input output routines
     -
   * - `Visualization`
     - Mayavi-based post-processing modules (`postproc.py`)
     -   
   * - `interactive/`
     - setup of IPython-based shell `isfepy`
     -
   * - `linalg/`
     - linear algebra functions not covered by NumPy and SciPy
     -

.. _python_packaging:

Python Packaging
----------------

OpenPNM shoudl be a self consistent python package, similar to

Packaging tutorials
^^^^^^^^^^^^^^^^^^^


.. list-table:: Packaging tutorials.
   :widths: 10 90
   :header-rows: 1
   
   * - Website
     - Description
   * - `classes <http://docs.python.org/2/tutorial/classes.html>`__
     - A basic intro into python classes
   * - `basic <http://guide.python-distribute.org/creation.html>`__
     - A basic intro into package construction
   * - `example <http://www.blog.pythonlibrary.org/2012/07/08/python-201-creating-modules-and-packages/>`__
     - A short and basic example
   * - `distutils <http://docs.python.org/2/distutils/setupscript.html>`__
     - Distutils and setup scripts examples


.. _coding_style:

Coding style
------------

All the code in SfePy should try to adhere to python style guidelines, see
`PEP-0008`.

There are some additional recommendations:

- Prefer whole words to abbreviations in public APIs - there is completion
  after all. If some abbreviation is needed (*really* too long name), try to
  make it as comprehensible as possible. Also check the code for similar
  names - try to name things consistently with the existing code. Examples:

  - yes: ``equation``, ``transform_variables()``, ``filename``
  - rather not: ``eq``, ``transvar()``, ``fname``

- Functions have usually form ``<action>_<subject>()`` e.g.: ``save_data()``,
  ``transform_variables()``, do not use ``data_save()``,
  ``variable_transform()`` etc.
- Variables like ``V``, ``c``, ``A``, ``b``, ``x`` should be tolerated only
  locally when expressing mathematical ideas.

Really minor recommendations:

- Avoid single letter names, if you can:

  - not even for loop variables - use e.g. ir, ic, ... instead of i, j for rows
    and columns
  - not even in generators, as they "leak" (this is fixed in Python 3.x)

These are recommendations only, we will not refuse code just on the ground that
it uses slightly different formatting, as long as it follows the PEP.

Note: some old parts of the code might not follow the PEP, yet. We fix them
progressively as we update the code.



Docstring standard
^^^^^^^^^^^^^^^^^^

We use `sphinx` with the `numpydoc` extension to generate this
documentation. Refer to the sphinx site for the possible markup constructs.

Basically (with a little tweak), we try to follow the NumPy/SciPy docstring
standard as described in `NumPy documentation guide`. See also the complete
`docstring example`. It is exaggerated a bit to show all the
possibilities. Use your common sense here - the docstring should be sufficient
for a new user to use the documented object. A good way to remember the format
is to type::

    In [1]: import numpy as nm
    In [2]: nm.sin?

in `ipython`. The little tweak mentioned above is the starting newline::

    def function(arg1, arg2):
        """
	This is a function.

        Parameters
        ----------
        arg1 : array
            The coordinates of ...
        arg2 : int
            The dimension ...

        Returns
        -------
        out : array
           The resulting array of shape ....
        """

It seems visually better than::

    def function(arg1, arg2):
        """This is a function.

        Parameters
        ----------
        arg1 : array
            The coordinates of ...
        arg2 : int
            The dimension ...

        Returns
        -------
        out : array
           The resulting array of shape ....
        """

When using :math:`\mbox{\LaTeX}` in a docstring, use a raw string::

    def function():
        r"""
	This is a function with :math:`\mbox{\LaTeX}` math:
        :math:`\frac{1}{\pi}`.
	"""

to prevent Python from interpreting and consuming the backslashes in common
escape sequences like '\\n', '\\f' etc.

.. _how_to_regenerate_documentation:

How to Regenerate Documentation
-------------------------------

The following steps summarize how to regenerate this documentation.

#. Install `sphinx` and `numpydoc`. Do not forget to set the path to numpydoc
   in site_cfg.py if it is not installed in a standard location for Python
   packages on your platform. A recent :math:`\mbox{\LaTeX}` distribution is
   required, too, for example `TeX Live`. Depending on your OS/platform, it
   can be in the form of one or several packages.

#. Edit the rst files in `doc/` directory using your favorite text editor - the
   ReST format is really simple, so nothing fancy is needed. Follow the
   existing files in `doc/`; for reference also check `reStructuredText
   Primer`, `Sphinx Markup Constructs` and `docutils reStructuredText`.

   - When adding a new Python module, add a corresponding documentation file
     into `doc/src/sfepy/<path>`, where `<path>` should reflect the location of
     the module in `sfepy/`.

   - Figures belong to `doc/images`; subdirectories can be used.

#. (Re)generate the documentation (assuming GNU make is installed)::

    cd doc
    make html

#. View it (substitute your favorite browser)::

    firefox _build/html/index.html

