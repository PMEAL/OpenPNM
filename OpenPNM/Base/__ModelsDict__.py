'''
###############################################################################
ModelsDict:  Abstract Class for Containing Models
###############################################################################
'''
import scipy as sp
from collections import OrderedDict
from OpenPNM.Base import logging, Controller
logger = logging.getLogger()

class GenericModel(dict):
    r"""
    Accepts a model from the OpenPNM model library, as well as all required 
    and optional argumnents, then wraps it in a custom dictionary with 
    various methods for working with the models.
    
    Notes
    -----
    This class in ONLY used inside the ModelsDict class, so it is a
    nested class.  This may change in the future if we need it elsewhere.

    """
    def __init__(self,**kwargs):
        self.update(**kwargs)

    def __call__(self):
        return self['model'](**self)
        
    def __str__(self):
        header = '-'*60
        print(header)
        print(self['model'].__module__+'.'+self['model'].__name__)
        print(header)
        print("{a:<20s} {b}".format(a='Argument Name',b='Value'))
        print(header)
        for item in self.keys():
            if item not in ['model','network','geometry','phase','physics','propname']:
                print("{a:<20s} {b}".format(a=item, b=self[item]))
        print(header)
        return ' '
        
    def regenerate(self):
        r'''
        '''
        master = self._find_master()
        #Determine object type, and assign associated objects
        self_type = [item.__name__ for item in master.__class__.__mro__]
        kwargs = {}
        if 'GenericGeometry' in self_type:
            kwargs['network'] = master._net
            kwargs['geometry'] = master
        elif 'GenericPhase' in self_type:
            kwargs['network'] = master._net
            kwargs['phase'] = master
        elif 'GenericPhysics' in self_type:
            kwargs['network'] = master._net
            kwargs['phase'] = master._phases[0]
            kwargs['physics'] = master
        else:
            kwargs['network'] = master
        kwargs.update(self)
        return self['model'](**kwargs)
        
    def _find_master(self):
        sim = Controller()
        master = []
        for item in sim.keys():
            if sim[item].models is not None:
                for model in sim[item].models.keys():
                    if sim[item].models[model] is self:
                        master.append(sim[item])
        if len(master) > 1:
            raise Exception('More than one master found! This model dictionary has been associated with multiple objects. To use the same dictionary multiple times use the copy method.')
        return master[0]

class ModelsDict(OrderedDict):
    r"""
    Accepts a model from the OpenPNM model library, as well as all required and
    optional argumnents, then wraps it in a custom dictionary with various 
    methods for working with the models.

    """
    
    def __setitem__(self,propname,model):
        temp = GenericModel(propname=propname,model=None)
        temp.update(**model)
        super(ModelsDict,self).__setitem__(propname,temp)
        
    def __str__(self):
        header = '-'*60
        print(header)
        print("{n:<5s} {a:<30s} {b:<20s}".format(n='#', a='Property Name', b='Regeneration Mode'))
        print(header)
        count = 0
        for item in self.keys():
            print("{n:<5d} {a:<30s} {b:<20s}".format(n=count, a=item, b=self[item]['regen_mode']))
            count += 1
        print(header)
        return ' '
        
    def keys(self):
        return list(super(ModelsDict,self).keys())
            
    def regenerate(self, props='', mode='inclusive'):
        r'''
        This updates properties using any models on the object that were
        assigned using ``add_model``

        Parameters
        ----------
        props : string or list of strings
            The names of the properties that should be updated, defaults to 'all'
        mode : string
            This controls which props are regenerated and how.  Options are:

            * 'inclusive': (default) This regenerates all given properties
            * 'exclude': This generates all given properties EXCEPT the given ones

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.pores(),throats=pn.throats())
        >>> geom['pore.diameter'] = 1
        >>> import OpenPNM.Geometry.models as gm  # Import Geometry model library
        >>> f = gm.pore_area.cubic
        >>> geom.add_model(propname='pore.area',model=f)  # Add model to Geometry object
        >>> geom['pore.area'][0]  # Look at area value in pore 0
        1
        >>> geom['pore.diameter'] = 2
        >>> geom.models.regenerate()  # Regenerate all models
        >>> geom['pore.area'][0]  # Look at pore area calculated with new diameter
        4

        '''
        master = self._find_master()
        if props == '':  # If empty, assume all models are to be regenerated
            props = list(self.keys())
            for item in props:  # Remove models if they are meant to be regenerated 'on_demand' only
                if self[item]['regen_mode'] == 'on_demand':
                    props.remove(item)
        elif type(props) == str:
            props = [props]
        if mode == 'exclude':
            temp = list(self.keys())
            for item in props:
                temp.remove(item)
            props = temp
        for item in self.keys():
            if self[item]['regen_mode'] == 'constant':
                props.remove(item)
        logger.info('Models are being recalculated in the following order: ')
        count = 0
        for item in props:
            if item in list(self.keys()):
                master[item] = self[item].regenerate()
                logger.info(str(count)+' : '+item)
                count += 1
            else:
                logger.warning('Requested proptery is not a dynamic model: '+item)
            
    def add(self,propname,model,regen_mode='normal',**kwargs):
        r'''
        Add specified property estimation model to the object.

        Parameters
        ----------
        propname : string
            The name of the property to use as dictionary key, such as
            'pore.diameter' or 'throat.length'

        model : function
            The property estimation function to use

        regen_mode : string
            Controls when and if the property is regenerated. Options are:

            * 'normal' : The property is stored as static data and is only regenerated when the object's ``regenerate`` is called

            * 'constant' : The property is calculated once when this method is first run, but always maintains the same value

            * 'deferred' : The model is stored on the object but not run until ``regenerate`` is called

            * 'on_demand' : The model is stored on the object but not run, AND will only run if specifically requested in ``regenerate``

        Notes
        -----
        This method is inherited by all net/geom/phys/phase objects.  It takes
        the received model and stores it on the object under private dictionary
        called _models.  This dict is an 'OrderedDict', so that the models can
        be run in the same order they are added.

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn)
        >>> import OpenPNM.Geometry.models as gm
        >>> f = gm.pore_misc.random  # Get model from Geometry library
        >>> geom.add_model(propname='pore.seed',model=f)
        >>> print(geom.models)  # Look in private dict to verify model was added
        ['pore.seed']

        '''
        master = self._find_master()
        if master == None:
            logger.warning('ModelsDict has no master, changing regen_mode to deferred')
            regen_mode = 'deferred'
        #Build dictionary containing given kwargs plus other required info
        f = {'model':model,'regen_mode':regen_mode}
        f.update(**kwargs)
        # Add model to ModelsDict
        self[propname] = f
        # Now generate data as necessary
        if regen_mode in ['normal','constant']:
            master[propname] = self[propname].regenerate()
        if regen_mode in ['deferred','on_demand']:
            pass

    def reorder(self,new_order):
        r'''
        Reorders the models on the object to change the order in which they
        are regenerated, where item 0 is calculated first.

        Parameters
        ----------
        new_order : dict
            A dictionary containing the model name(s) as the key, and the
            location(s) in the new order as the value

        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.TestGeometry(network=pn,pores=pn.Ps,throats=pn.Ts)
        >>> print(list(geom.models))
        ['pore.seed', 'throat.seed', 'throat.length']
        >>> geom.models.reorder({'pore.seed':1,'throat.length':0})
        ['throat.length', 'pore.seed', 'throat.seed']

        '''
        #Generate numbered list of current models
        order = [item for item in list(self.keys())]
        #Remove supplied models from list
        for item in new_order:
            order.remove(item)
        #Add models back to list in new order
        inv_dict = {v: k for k, v in new_order.items()}
        for item in inv_dict:
            order.insert(item,inv_dict[item])
        #Now rebuild models OrderedDict in new order
        for item in order:
            self.move_to_end(item)
        
    def _find_master(self):
        sim = Controller()
        master = []
        for item in sim.keys():
            if sim[item].models is self:
                master.append(sim[item])
        if len(master) > 1:
            raise Exception('More than one master found! This model dictionary has been associated with multiple objects. To use the same dictionary multiple times use the copy method.')
        return master[0]

if __name__ == '__main__':
    pn = OpenPNM.Network.TestNet()
    a = ModelsDict()



