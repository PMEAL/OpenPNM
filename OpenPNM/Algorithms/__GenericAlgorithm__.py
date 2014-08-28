"""
module __GenericAlgorithm__: Base class to build custom algorithms
==================================================================

.. warning:: The classes of this module should be loaded through the 'Algorithms/__init__.py' file.

"""
import OpenPNM
import sys

class GenericAlgorithm(OpenPNM.Base.Core):
    r"""
    GenericAlgorithm - Base class to execute algorithms

    Parameters
    ----------
    network : OpenPNM Network Object
        The network object to which this algorithm will apply.  
        
    name : string, optional
        Name of this algorithm


    Notes
    -----
    If no network is supplied an empty algorithm object is returned.  This is
    useful for loading in a saved algorithm from memory.

    """

    def __init__(self,network=None,name=None,**kwords):
        r"""
        Initialize
        """
        super(GenericAlgorithm,self).__init__(**kwords)
        self._logger.debug("Construct class")
        
        if network == None:
            self._net = OpenPNM.Network.GenericNetwork()
        else:
            self._net = network
        self.name = name

        # Initialize label 'all' in the object's own info dictionaries
        self['pore.all'] = self._net['pore.all']
        self['throat.all'] = self._net['throat.all']

    def run(self,**params):
        r"""
        Main run command for the algorithm
        """
        self._logger.debug(sys._getframe().f_code.co_name)
        self._do_outer_iteration_stage()

    def _do_outer_iteration_stage(self):
        r"""
        Executes the outer iteration stage
        """
        self._logger.debug(sys._getframe().f_code.co_name)
        self._do_one_outer_iteration()
        
    def _do_one_outer_iteration(self):
        r"""
        One iteration of an outer iteration loop for an algorithm
        (e.g. time or parametric study)
        """
        self._logger.debug(sys._getframe().f_code.co_name)
        self._do_inner_iteration_stage()

    def _do_inner_iteration_stage(self):
        r"""
        Executes the inner iteration stage
        """
        self._logger.debug(sys._getframe().f_code.co_name)
        self._do_one_inner_iteration()
        
    def _do_one_inner_iteration(self):
        r"""
        Executes one inner iteration
        """
        self._logger.debug(sys._getframe().f_code.co_name)
        
    def return_results(self,**kwargs):
        self._logger.debug(sys._getframe().f_code.co_name)
        
    def update(self,**kwargs):
        self.return_results(**kwargs)



if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Algorithms.GenericAlgorithm(network=pn,loglevel=10)
    test.run()
