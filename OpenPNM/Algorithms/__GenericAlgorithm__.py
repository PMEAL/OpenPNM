"""
module __GenericAlgorithm__: Base class to construct pore networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms/__init__.py' file.

"""
import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
from time import clock
import heapq
import itertools

class GenericAlgorithm(OpenPNM.Utilities.Tools):
    r"""
    GenericAlgorithm - Base class to execute algorithms

    Parameters
    ----------
    network : Descendent of OpenPNM.Network.GenericNetwork
        A valid network for this algorithm
    name : name of this algorithm



    Examples
    --------
    >>> print('nothing yet')

    .. note::
    n/a

    """

    def __init__(self,name,network=None,**kwords):
        r"""
        Initialize
        """
        super(GenericAlgorithm,self).__init__(**kwords)
        self._logger.debug("Construct class")
        self.name = name
        self._net = network
        self.set_pore_info(label='all',locations=self._net.get_pore_info(label='all')) #This is necessary for the methods from 'tools' to work.  They must know network size.
        self.set_throat_info(label='all',locations=self._net.get_throat_info(label='all'))  

    def run(self,**params):
        r"""
        Main run command for the algorithm
        """
#        self._logger.info("Execute run(): Basic version")
        self._setup(**params)
        self._do_outer_iteration_stage()

    def _setup(self, **params):
        r"""
        Perform applicable preliminary checks and calculations required for algorithm

        Notes
        -----
        This method is not implemented in the GenericAlgorithm, and must be subclassed in each algorithm as necessary
        """
#        self._logger.error("_setup: not implemented")

    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm
        (e.g. time or parametric study)
        """
#        self._logger.info("One Outer Iteration: Basic version")
        self._do_inner_iteration_stage()

    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
#        self._logger.info("Outer Iteration Stage: Basic version")
        self._do_one_outer_iteration()

    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
#        self._logger.warning("One Inner Iteration: Implement me")

    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
#        self._logger.info("Inner Iteration Stage: Basic version")
        self._do_one_inner_iteration()
        
    def update(self,**kwargs):
        print('not implemented')

    def set_result(self,**kwargs):
        
        self.update(**kwargs)

if __name__ =="__main__":
    print('    ************Testing Generic Algorithm**************')
    pn = OpenPNM.Geometry.Cubic().generate()
    test = GenericAlgorithm(loggername="TestGenericAlg")
    test.run(pn)
