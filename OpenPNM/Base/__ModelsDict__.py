"""
###############################################################################
ModelsDict:  Abstract Class for Containing Models
###############################################################################
"""
import inspect
from collections import OrderedDict
from OpenPNM.Base import logging, Controller
logger = logging.getLogger()


class ModelWrapper(dict):
    r"""
    Accepts a model from the OpenPNM model library, as well as all required
    and optional argumnents, then wraps it in a custom dictionary with
    various methods for working with the models.

    """

    COMPONENTS = ['model', 'network', 'geometry', 'phase', 'physics', 'propname']

    def __init__(self, **kwargs):
        self.update(**kwargs)

    def __call__(self):
        return self['model'](**self)

    def __str__(self):
        if self['model'] is None:
            return 'No model specified.'
        horizontal_rule = '-' * 60
        lines = [horizontal_rule]
        lines.append(self['model'].__module__ + '.' + self['model'].__name__)
        lines.append(horizontal_rule)
        lines.append('{0:<20s} {1}'.format('Argument Name', 'Value / (Default)'))
        lines.append(horizontal_rule)
        # Scan default argument names and values of model
        defs = {}
        if self['model'].__defaults__ is not None:
            vals = list(inspect.getargspec(self['model']).defaults)
            keys = inspect.getargspec(self['model']).args[-len(vals):]
            # Put defaults into the dict
            defs.update(zip(keys, vals))
        keys = list(self.keys())
        keys.sort()
        for item in keys:
            if item not in self.COMPONENTS:
                if item not in defs.keys():
                    defs[item] = '---'
                lines.append('{0:<20s} {1} / ({2})'.format(item,
                                                           self[item],
                                                           defs[item]))
        lines.append(horizontal_rule)
        return '\n'.join(lines)

    def regenerate(self):
        r"""
        Regenerate the model
        """
        master = self._find_master()
        # Determine object type, and assign associated objects
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
        ctrl = Controller()
        master = []
        for item in ctrl.keys():
            if ctrl[item].models is not None:
                for model in ctrl[item].models.keys():
                    if ctrl[item].models[model] is self:
                        master.append(ctrl[item])
        if len(master) > 1:
            raise Exception('More than one master found! This model dictionary '
                            'has been associated with multiple objects. To use '
                            'the same dictionary multiple times use the copy'
                            'method.')
        return master[0]


class GenericModel(ModelWrapper):
    r"""
    This class was deprecated, and replaced by ModelWrapper.  Unfortunately,
    this broke the ability to load files that were saved prior to the
    deprecation.  This is placed here as a 'wrapper' to the new ModelsWraper.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)


class ModelsDict(OrderedDict):
    r"""
    This custom dictionary stores the models that are associated with each
    OpenPNM object.  This is an ordered dict with a few additional methods.
    This ModelsDict class can be created as a standalone object, then
    associated with an OpenPNM object, and ModelsDicts from one object can
    be copied and attached to another.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> Ps = pn.pores(labels='top',mode='not')
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, throats=pn.Ts)

    Add a model to the object's ModelsDict:

    >>> f = OpenPNM.Geometry.models.pore_seed.random
    >>> geom.models.add(propname='pore.seed', model=f, seed=0)

    It is possible to use the ModelsDict from one object with another object:

    >>> Ps = pn.pores('top')
    >>> boun = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps)
    >>> boun.models  # The boun object has no models in its Models dict
    ModelsDict()
    >>> mod = geom.models.copy()  # Create a copy of geom's models
    >>> boun.models = mod  # Use the same set of models on boun as geom

    Because ``copy`` was used above, the ModelsDict on the two objects are
    equal but different:

    >>> boun.models == geom.models
    True
    >>> boun.models is geom.models
    False
    """

    def __setitem__(self, propname, model):
        temp = ModelWrapper(propname=propname, model=None)
        temp.update(**model)
        super().__setitem__(propname, temp)

    def __str__(self):
        horizontal_rule = '-' * 60
        lines = [horizontal_rule]
        lines.append('{0:<5s} {1:<30s} {2}'.format('#',
                                                   'Property Name',
                                                   'Regeneration Mode'))
        lines.append(horizontal_rule)
        for i, item in enumerate(self.keys()):
            str = '{0:<5d} {1:<30s} {2:<20s}'
            lines.append(str.format(i + 1, item, self[item]['regen_mode']))
        lines.append(horizontal_rule)
        return '\n'.join(lines)

    def keys(self):
        return list(super().keys())

    def regenerate(self, props='', mode='inclusive'):
        r"""
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
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
        ...                                         pores=pn.pores(),
        ...                                         throats=pn.throats())
        >>> geom['pore.diameter'] = 1
        >>> import OpenPNM.Geometry.models as gm  # Import Geometry model library
        >>> f = gm.pore_area.cubic
        >>> geom.add_model(propname='pore.area', model=f)
        >>> geom['pore.area'][0]  # Look at area value in pore 0
        1
        >>> geom['pore.diameter'] = 2
        >>> geom.models.regenerate()  # Regenerate all models
        >>> geom['pore.area'][0]  # Look at pore area calculated with new diameter
        4

        """

        master = self._find_master()
        # If empty, assume all models are to be regenerated
        if props == '':
            props = list(self.keys())
            # Remove models if they are meant to be regenerated 'on_demand' only
            for item in props:
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
                logger.info(str(count) + ' : ' + item)
                count += 1
            else:
                logger.warning('Requested proptery is not a dynamic model: ' + item)

    def add(self, propname, model, regen_mode='normal', **kwargs):
        r"""
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

            * 'normal' : The property is stored as static data and is only
                         regenerated when the object's ``regenerate`` is called

            * 'constant' : The property is calculated once when this method is first
                           run, but always maintains the same value

            * 'deferred' : The model is stored on the object but not run until
                           ``regenerate`` is called

            * 'on_demand' : The model is stored on the object but not run, AND will
                            only run if specifically requested in ``regenerate``

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
        >>> geom.models.keys()  # Look in dict to verify model was added
        ['pore.seed']
        >>> print(geom.models['pore.seed'])  # Look at arguments for model
        ------------------------------------------------------------
        OpenPNM.Geometry.models.pore_misc.random
        ------------------------------------------------------------
        Argument Name        Value / (Default)
        ------------------------------------------------------------
        num_range            [0, 1] / ([0, 1])
        regen_mode           normal / (---)
        seed                 None / (None)
        ------------------------------------------------------------
        """
        master = self._find_master()
        if master is None:
            logger.warning('ModelsDict has no master, changing regen_mode \
                            to deferred')
            regen_mode = 'deferred'
        # Build dictionary containing default model values, plus other required info
        f = {'model': model, 'regen_mode': regen_mode}
        # Scan default argument names and values of model
        if model.__defaults__ is not None:
            vals = list(inspect.getargspec(model).defaults)
            keys = inspect.getargspec(model).args[-len(vals):]
            # Put defaults into the dict
            f.update(zip(keys, vals))
        # Update dictionary with supplied arguments, overwriting defaults
        f.update(**kwargs)
        # Add model to ModelsDict
        self[propname] = f
        # Now generate data as necessary
        if regen_mode in ['normal', 'constant']:
            master[propname] = self[propname].regenerate()
        if regen_mode in ['deferred', 'on_demand']:
            pass

    def remove(self, propname):
        r"""
        Removes selected model from the dictionary, as well as removing its
        associated data from the master Core object.

        Parameters
        ----------
        propname : string
            The name of the model to remove
        """
        master = self._find_master()
        temp = master.pop(propname, None)
        del self[propname]

    def reorder(self, new_order):
        r"""
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
        >>> geom = OpenPNM.Geometry.TestGeometry(network=pn,
        ...                                      pores=pn.Ps,
        ...                                      throats=pn.Ts)
        >>> geom.models.keys()
        ['pore.seed', 'throat.seed', 'throat.length']
        >>> geom.models.reorder({'pore.seed': 1, 'throat.length': 0})
        >>> geom.models.keys()
        ['throat.length', 'pore.seed', 'throat.seed']

        """
        # Generate numbered list of current models
        order = [item for item in list(self.keys())]
        # Remove supplied models from list
        for item in new_order:
            order.remove(item)
        # Add models back to list in new order
        inv_dict = {v: k for k, v in new_order.items()}
        for item in inv_dict:
            order.insert(item, inv_dict[item])
        # Now rebuild models OrderedDict in new order
        for item in order:
            self.move_to_end(item)

    def _find_master(self):
        ctrl = Controller()
        master = []
        for item in ctrl.keys():
            if ctrl[item].models is self:
                master.append(ctrl[item])
        if len(master) > 1:
            raise Exception('More than one master found! This model dictionary '
                            'has been associated with multiple objects. To use the '
                            'same dictionary multiple times use the copy method.')
        elif len(master) == 0:
            raise Exception('ModelsDict has no master.')

        return master[0]
