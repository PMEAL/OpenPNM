"""
===============================================================================
module __GenericAlgorithm__: Base class to build custom algorithms
==================================================================

This generic class contains the recommended methods for subclassed algorithms.
It inherits from Core, so is Python Dict with the OpenPNM data control methods.

"""
import OpenPNM
import sys
import scipy as sp

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

    def update_results(self,**kwargs):
        self._logger.debug(sys._getframe().f_code.co_name)

    def set_boundary_conditions(self,component=None,bctype='',bcvalue=[],pores=[],throats=[],mode='merge'):
        r'''
        Apply boundary conditions to specified pores or throats

        Parameters
        ----------
        bctype : string
            Specifies the type or the name of boundary condition to apply.  The types can be one one of the followings:

            - 'Dirichlet' : Specify the quantity in each location
            - 'Neumann' : Specify the flow rate into each location
            - 'Neumann_group' : Specify the net flow rate into a group of pores/throats
        component : OpenPNM Phase object
            The Phase object to which this BC applies
        bcvalue : array_like
            The boundary value to apply, such as concentration or rate
        pores : array_like
            The pores where the boundary conditions should be applied
        throats : array_like
            The throats where the boundary conditions should be applied
        mode : string, optional
            Controls how the conditions are applied.  Options are:

            - 'merge': Inserts the specified values, leaving existing values elsewhere
            - 'overwrite': Inserts specified values, clearing all other values
            - 'remove': Removes boundary conditions from specified locations

        Notes
        -----
        1. It is not possible to have multiple boundary conditions for a specified location in just one algorithm. 
        So when new condition is going to be applied to a specific location, any existing one
        should be removed or overwritten..
        '''
        try: self._existing_BC
        except: self._existing_BC = []
        if component==None:
            try:
                component = self._phase
            except:
                raise Exception('For using set_boundary_conditions method, a component/phase should be specified.')
        self._logger.debug('BC applies to the component: '+component.name)

        #If mode is 'remove', also bypass checks
        if mode == 'remove':
            if pores == [] and throats == []:
                if bctype=='':
                    raise Exception('No bctype/pore/throat is specified')
                else:
                    for item in self.labels():
                        if bctype == (item.split('.')[-1]).replace(self._phase.name+'_',"") :
                            element = item.split('.')[0]
                            try:
                                del self[element+'.'+component.name+'_bcval_'+bctype]
                            except: pass
                            del self[item]
                    self._logger.debug('Removing '+bctype+' from all locations in '+self.name)
                    self._existing_BC.remove(bctype)
            else:
                if pores!=[]:
                    if bctype!='':
                        self['pore.'+component.name+'_bcval_'+bctype][pores] = sp.nan
                        self['pore.'+component.name+'_'+bctype][pores] = False
                        self._logger.debug('Removing '+bctype+' from the specified pores in '+self.name)
                    else:   raise Exception('Cannot remove BC from the pores unless bctype is specified')

                if throats!=[]:
                    if bctype!='':
                        self['throat.'+component.name+'_bcval_'+bctype][throats] = sp.nan
                        self['throat.'+component.name+'_'+bctype][throats] = False
                        self._logger.debug('Removing '+bctype+' from the specified throats in '+self.name)
                    else:   raise Exception('Cannot remove BC from the throats unless bctype is specified')

            return
        #Validate bctype
        if bctype == '':
            raise Exception('bctype must be specified')
        #Validate pores/throats
        if pores == [] and throats == []:
            raise Exception('pores/throats must be specified')
        elif pores!=[] and throats!=[]:
            raise Exception('BC for pores and throats must be specified independently')
        elif  pores!=[]:
            element ='pore'
            loc = sp.array(pores,ndmin=1)
            all_length = self.num_pores()
        if throats!=[]:
            element ='throat'
            loc = sp.array(throats,ndmin=1)
            all_length = self.num_throats()
        #Validate bcvalue
        if bcvalue != []:
            #Check bcvalues are compatible with bctypes
            if bctype == 'Neumann_group':  #Only scalars are acceptable
                if sp.size(bcvalue) != 1:
                    raise Exception('When specifying Neumann_group, bcval should be a scalar')
            else: #Only scalars or Np/Nt-long are acceptable
                if sp.size(bcvalue) == 1:
                    bcvalue = sp.ones(sp.shape(loc))*bcvalue
                elif sp.size(bcvalue) != sp.size(loc):
                    raise Exception('The pore/throat list and bcvalue list are different lengths')

        #Confirm that prop and label arrays exist
        if element+'.'+component.name+'_bcval_'+bctype not in self.props():
            self[element+'.'+component.name+'_bcval_'+bctype] = sp.ones((all_length,),dtype=float)*sp.nan
        if element+'.'+component.name+'_'+bctype not in self.labels():
            self[element+'.'+component.name+'_'+bctype] = sp.zeros((all_length,),dtype=bool)

        #Check all BC from specified locations, prior to setting new ones
        for item in self.labels():
            bcname = (item.split('.')[-1]).replace(self._phase.name+'_',"")
            if bcname in self._existing_BC  and item.split('.')[0]==element:
                if mode=='merge':
                    if not (sp.isnan(self[element+'.'+component.name+'_bcval_'+bcname][loc]).all() and sp.sum(self[element+'.'+component.name+'_'+bcname][loc])==0):
                        raise Exception('Because of the existing BCs, the method cannot apply new BC with the merge mode to the specified pore/throat.')
        #Set boundary conditions based on supplied mode
        if mode == 'merge':
            if bcvalue != []:   self[element+'.'+component.name+'_bcval_'+bctype][loc] = bcvalue
            self[element+'.'+component.name+'_'+bctype][loc] = True
            if bctype not in self._existing_BC: self._existing_BC.append(bctype)
        elif mode == 'overwrite':
            self[element+'.'+component.name+'_bcval_'+bctype] = sp.ones((all_length,),dtype=float)*sp.nan
            if bcvalue != []:   self[element+'.'+component.name+'_bcval_'+bctype][loc] = bcvalue
            self[element+'.'+component.name+'_'+bctype] = sp.zeros((all_length,),dtype=bool)
            self[element+'.'+component.name+'_'+bctype][loc] = True
            if bctype not in self._existing_BC: self._existing_BC.append(bctype)

if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    test = OpenPNM.Algorithms.GenericAlgorithm(network=pn,loglevel=10)
    test.run()
