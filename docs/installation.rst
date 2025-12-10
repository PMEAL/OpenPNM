.. _installation:

Installation
############

``uv``
------

The recommended way to install Python and desired packages (like OpenPNM) is using `uv <https://docs.astral.sh/uv>`_. ``uv`` has rapidly become the most popular package manager for Python. You'll be able to find lots of blog articles extolling the virtues of ``uv``, but its main feature is speed.

1. First you need to install ``uv``. The official installation instructions are provided on the `astral website <https://docs.astral.sh/uv/getting-started/installation>`_. This will install ``uv`` on your system so it will be available from the terminal or command line.
2. Next you navigate to the folder where the project files and data will be or already are stored, like scripts and tomograms.
3. Then you create a virtual environment in that directory using ``uv venv --python 3.12``. This adds a ``.venv`` folder, where ``uv`` will store all its information. Note that OpenPNM only supports python version 3.11, 3.12 and 3.13. This is largely because of dependencies which have firm version requirements.
4. Finally, you run `uv pip install openpnm`. The first time you do this ``uv`` will download and compile a few things, which may take some time, but it will store all of this so subsequent usage will be much faster.

When using ``uv`` we recommend `VSCode <https://code.visualstudio.com/>`_ as the IDE to write and edit scripts. VSCode will automatically find the ``venv`` inside the current folder, which is very handy. Spyder is a good IDE for scientific programming, but it does not (yet) have automatic support for finding ``venvs`` so is cumbersome to use when switching between projects.

``conda``
---------

OpenPNM is also available via the `conda <https://docs.conda.io/en/latest/>`_ package manager. If you prefer to use ``conda``, you can install OpenPNM by running the following command in your terminal or command prompt:

.. code-block:: bash

   $ conda install conda-forge::openpnm

Installing the ``dev`` version
##############################

If you are a OpenPNM contributor or want to get the newest updates as they roll in, you need to clone the OpenPNM repository from Github and install it locally. It's not as difficult as it sounds, just follow these steps:

Open up the terminal/cmd and use ``cd`` to navigate to the directory where you want to store the OpenPNM code. Clone the repo to your disk using:

.. code-block:: bash

   $ git clone https://github.com/PMEAL/OpenPNM

Since your terminal is currently in the OpenPNM directory, you can now install OpenPNM and all its dependencies with:

.. code-block:: bash

   $ uv pip install -e .

Voila! You can now use the latest features available on the ``dev`` branch. To keep your "local" OpenPNM installation up to date, you should regularly pull the latest changes:

.. code-block:: bash

   $ git pull

.. warning::
   For the development version of OpenPNM to work, you need to first remove
   any versions that you may have previously installed using ``uv pip uninstall openpnm``
