from collections import namedtuple
import numpy as np
from openpnm.utils import PrintableDict, PrintableList

__all__ = ['ParamMixin', 'LabelMixin', 'LegacyMixin']


class ParamMixin:
    """Brief explanation of ParamMixin"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._params = PrintableDict()
        self._params._key = "Parameters"
        self._params._value = "Values"

    def __getitem__(self, key):
        # If the key is a just a numerical value, the kick it directly back
        # This allows one to do either value='pore.blah' or value=1.0
        if isinstance(key, (int, float, bool, complex)):
            return key

        if key.startswith('param'):
            try:
                vals = self._params[key]
            except KeyError:
                vals = self.network._params[key]
        else:
            vals = super().__getitem__(key)
        return vals

    def __setitem__(self, key, value):
        if key.startswith('param'):
            self._params[key] = value
        else:
            super().__setitem__(key, value)

    def __str__(self):
        s = super().__str__()
        s = s.rpartition('\n')[0]
        s = s + '\n' + self._params.__str__()
        return s

    def params(self):
        r"""
        Return parameter names and values in a dictionary
        """
        return self._params


class LegacyMixin:
    """Brief explanation of LegacyMixin"""

    def tomask(self, *args, **kwargs):
        """Brief explanation of tomask"""
        return self.to_mask(*args, **kwargs)

    def toindices(self, *args, **kwargs):
        """Brief explanation of tomask"""
        return self.to_indices(*args, **kwargs)

    def _map(self, ids, element, filtered):
        ids = np.array(ids, dtype=np.int64)
        locations = self._get_indices(element=element)
        self_in_ids = np.isin(ids, self[element+'._id'], assume_unique=True)
        ids_in_self = np.isin(self[element+'._id'], ids, assume_unique=True)
        mask = np.zeros(shape=ids.shape, dtype=bool)
        mask[self_in_ids] = True
        ind = np.ones_like(mask, dtype=np.int64) * -1
        ind[self_in_ids] = locations[ids_in_self]
        if filtered:
            return ind[mask]
        t = namedtuple('index_map', ('indices', 'mask'))
        return t(ind, mask)

    def map_pores(self, pores, origin, filtered=True):
        r"""
        Given a list of pore on a target object, finds indices of those pores
        on the calling object

        Parameters
        ----------
        pores : array_like
            The indices of the pores on the object specifiedin ``origin``

        origin : Base
            The object corresponding to the indices given in ``pores``

        filtered : bool (default is ``True``)
            If ``True`` then a ndarray of indices is returned with missing
            indices removed, otherwise a named-tuple containing both the
            ``indices`` and a boolean ``mask`` with ``False`` indicating
            which locations were not found.

        Returns
        -------
        Pore indices on the calling object corresponding to the same pores
        on the ``origin`` object.  Can be an array or a tuple containing an
        array and a mask, depending on the value of ``filtered``.

        See Also
        --------
        pores
        map_throats

        """
        ids = origin['pore._id'][pores]
        return self._map(element='pore', ids=ids, filtered=filtered)

    def map_throats(self, throats, origin, filtered=True):
        r"""
        Given a list of throats on a target object, finds indices of
        those throats on the calling object

        Parameters
        ----------
        throats : array_like
            The indices of the throats on the object specified in ``origin``
        origin : Base
            The object corresponding to the indices given in ``throats``
        filtered : bool, default is ``True``
            If ``True`` then a ndarray of indices is returned with missing
            indices removed, otherwise a named-tuple containing both the
            ``indices`` and a boolean ``mask`` with ``False`` indicating
            which locations were not found.

        Returns
        -------
        ndarray
            Throat indices on the calling object corresponding to the same
            throats on the ``origin`` object.  Can be an array or a tuple
            containing an array and a mask, depending on the value of
            ``filtered``.

        See Also
        --------
        throats
        map_pores

        """
        ids = origin['throat._id'][throats]
        return self._map(element='throat', ids=ids, filtered=filtered)


class LabelMixin:
    """Brief explanation of LabelMixin"""

    def _get_labels(self, element, locations, mode):
        r"""
        This is the actual label getter method, but it should not be called
        directly.  Use ``labels`` instead.
        """
        # Parse inputs
        locations = self._parse_indices(locations)
        element = self._parse_element(element=element)
        # Collect list of all pore OR throat labels
        labels = self.keys(mode='labels', element=element)
        labels.sort()
        labels = np.array(labels)  # Convert to ndarray for following checks
        # Make an 2D array with locations in rows and labels in cols
        arr = np.vstack([self[item][locations] for item in labels]).T
        num_hits = np.sum(arr, axis=0)  # Number of locations with each label
        if mode in ['or', 'union', 'any']:
            temp = labels[num_hits > 0]
        elif mode in ['and', 'intersection']:
            temp = labels[num_hits == locations.size]
        elif mode in ['xor', 'exclusive_or']:
            temp = labels[num_hits == 1]
        elif mode in ['nor', 'not', 'none']:
            temp = labels[num_hits == 0]
        elif mode in ['nand']:
            temp = labels[num_hits == (locations.size - 1)]
        elif mode in ['xnor', 'nxor']:
            temp = labels[num_hits > 1]
        else:
            raise Exception('Unrecognized mode:'+str(mode))
        return PrintableList(temp)

    def labels(self, pores=[], throats=[], element=None, mode='union'):
        r"""
        Returns a list of labels present on the object

        Additionally, this function can return labels applied to a specified
        set of pores or throats

        Parameters
        ----------
        element : str
            Controls whether pore or throat labels are returned.  If empty then
            both are returned (default).

        pores (or throats) : array_like
            The pores (or throats) whose labels are sought.  If left empty a
            list containing all pore and throat labels is returned.

        mode : str, optional
            Controls how the query should be performed.  Only applicable
            when ``pores`` or ``throats`` are specified:

            ==============  =====================================================
            mode            meaning
            ==============  =====================================================
            'or'            Returns the labels that are assigned to *any* of the
                            given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the given
                            locations
            'xnor'          Labels that are present on *more than one* of the given
                            locations.
            ==============  =====================================================

        Returns
        -------
        A list containing the labels on the object.  If ``pores`` or
        ``throats`` are given, the results are filtered according to the
        specified ``mode``.

        See Also
        --------
        props
        keys

        Notes
        -----
        Technically, *'nand'* and *'xnor'* should also return pores with *none*
        of the labels but these are not included.  This makes the returned list
        more useful.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.labels(pores=[11, 12])
        ['pore.all', 'pore.internal', 'pore.left', 'pore.surface']
        """
        # Short-circuit query when no pores or throats are given
        if (np.size(pores) == 0) and (np.size(throats) == 0):
            labels = PrintableList(self.keys(element=element, mode='labels'))
        elif (np.size(pores) > 0) and (np.size(throats) > 0):
            raise Exception('Cannot perform label query on pores and '
                            + 'throats simultaneously')
        elif np.size(pores) > 0:
            labels = self._get_labels(element='pore', locations=pores,
                                      mode=mode)
        elif np.size(throats) > 0:
            labels = self._get_labels(element='throat', locations=throats,
                                      mode=mode)
        return labels

    def set_label(self, label, pores=None, throats=None, mode='add'):
        r"""
        Creates or updates a label array

        Parameters
        ----------
        label : str
            The label to apply to the specified locations
        pores : array_like
            A list of pore indices or a boolean mask of where given label
            should be added or removed (see ``mode``)
        throats : array_like
            A list of throat indices or a boolean mask of where given label
            should be added or removed (see ``mode``)
        mode : str
            Controls how the labels are handled.  Options are:

            =========== ======================================================
            mode        description
            =========== ======================================================
            'add'       (default) Adds the given label to the specified
                        locations while keeping existing labels

            'overwrite' Removes existing label from all locations before
                        adding the label in the specified locations

            'remove'    Removes the given label from the specified locations
                        leaving the remainder intact

            'purge'     Removes the specified label from the object completely.
                        This ignores the ``pores`` and ``throats`` arguments.

            'clear'     Sets all the labels to ``False`` but does not remove
                        the label array
            =========== ======================================================

        """
        self._parse_mode(mode=mode,
                         allowed=['add', 'overwrite', 'remove', 'purge',
                                  'clear'])

        if label.split('.')[0] in ['pore', 'throat']:
            label = label.split('.', 1)[1]

        if (pores is not None) and (throats is not None):
            self.set_label(label=label, pores=pores, mode=mode)
            self.set_label(label=label, throats=throats, mode=mode)
            return
        elif pores is not None:
            locs = self._parse_indices(pores)
            element = 'pore'
        elif throats is not None:
            locs = self._parse_indices(throats)
            element = 'throat'

        if mode == 'add':
            if element + '.' + label not in self.keys():
                self[element + '.' + label] = False
            self[element + '.' + label][locs] = True
        if mode == 'overwrite':
            self[element + '.' + label] = False
            self[element + '.' + label][locs] = True
        if mode == 'remove':
            self[element + '.' + label][locs] = False
        if mode == 'clear':
            self['pore' + '.' + label] = False
            self['throat' + '.' + label] = False
        if mode == 'purge':
            _ = self.pop('pore.' + label, None)
            _ = self.pop('throat.' + label, None)

    def _get_indices(self, element, labels='all', mode='or'):
        r"""
        This is the actual method for getting indices, but should not be called
        directly.  Use ``pores`` or ``throats`` instead.
        """
        # Parse and validate all input values.
        element = self._parse_element(element, single=True)
        labels = self._parse_labels(labels=labels, element=element)
        if element+'.all' not in self.keys():
            raise Exception('Cannot proceed without {}.all'.format(element))

        # Begin computing label array
        if mode in ['or', 'any', 'union']:
            union = np.zeros_like(self[element+'.all'], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                union = union + self[element+'.'+item.split('.')[-1]]
            ind = union
        elif mode in ['and', 'all', 'intersection']:
            intersect = np.ones_like(self[element+'.all'], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                intersect = intersect*self[element+'.'+item.split('.')[-1]]
            ind = intersect
        elif mode in ['xor', 'exclusive_or']:
            xor = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                xor = xor + np.int8(info)
            ind = (xor == 1)
        elif mode in ['nor', 'not', 'none']:
            nor = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                nor = nor + np.int8(info)
            ind = (nor == 0)
        elif mode in ['nand']:
            nand = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                nand = nand + np.int8(info)
            ind = (nand < len(labels)) * (nand > 0)
        elif mode in ['xnor', 'nxor']:
            xnor = np.zeros_like(self[element+'.all'], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.')[-1]]
                xnor = xnor + np.int8(info)
            ind = (xnor > 1)
        else:
            raise Exception('Unsupported mode: '+mode)
        # Extract indices from boolean mask
        ind = np.where(ind)[0]
        ind = ind.astype(dtype=int)
        return ind

    def pores(self, labels='all', mode='or', asmask=False, to_global=False):
        r"""
        Returns pore indicies where given labels exist, according to the logic
        specified by the ``mode`` argument.

        Parameters
        ----------
        labels : str or list[str]
            The label(s) whose pores locations are requested.  This argument
            also accepts '*' for wildcard searches.
        mode : str
            Specifies how the query should be performed.  The options are:

            ==============  =====================================================
            mode            meaning
            ==============  =====================================================
            'or'            Pores with *one or more* of the given labels are
                            returned. Also accepts 'union' and 'any'.
            'and'           Pores with all of the given labels are returned.
                            Also accepts 'intersection' and 'all'.
            'xor'           Pores with *only one* of the given labels are returned.
                            Also accepts 'exclusive_or'.
            'nor'           Pores with *none* of the given labels are returned.
                            Also accepts 'not' and 'none'.
            'nand'          Pores with *not all* of the given labels are
                            returned.
            'xnor'          Pores with *more than one* of the given labels are
                            returned.
            ==============  =====================================================

        asmask : bool
            If ``True`` then a boolean array of length Np is returned with
            ``True`` values indicating the pores that satisfy the query.
        to_global : bool
            If ``True``, the returned indices will be indexed relative to the
            full domain.  This only has an effect when the calling object
            is a Subdomain.

        Returns
        -------
        A Numpy array containing pore indices filtered by the logic specified
        in ``mode``.

        See Also
        --------
        throats

        Notes
        -----
        Technically, *nand* and *xnor* should also return pores with *none* of
        the labels but these are not included.  This makes the returned list
        more useful.

        To perform more complex or compound queries, you can opt to receive
        the result a a boolean mask (``asmask=True``), then manipulate the
        arrays manually.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> Ps = pn.pores(labels=['top', 'back'], mode='union')
        >>> Ps[:5]  # Look at first 5 pore indices
        array([ 4,  9, 14, 19, 20])
        >>> pn.pores(labels=['top', 'back'], mode='xnor')
        array([ 24,  49,  74,  99, 124])
        """
        ind = self._get_indices(element='pore', labels=labels, mode=mode)
        if to_global and hasattr(self, 'to_global'):
            ind = self.to_global(pores=ind)
            if asmask:
                ind = self._domain.to_mask(pores=ind)
        else:
            if asmask:
                ind = self.to_mask(pores=ind)
        return ind

    def throats(self, labels='all', mode='or', asmask=False, to_global=False):
        r"""
        Returns throat locations where given labels exist, according to the
        logic specified by the ``mode`` argument.

        Parameters
        ----------
        labels : str or list[str]
            The throat label(s) whose locations are requested.  If omitted,
            'all' throat inidices are returned.  This argument also accepts
            '*' for wildcard searches.
        mode : str
            Specifies how the query should be performed. The options are:

            ==============  =====================================================
            mode            meaning
            ==============  =====================================================
            'or'            Throats with *one or more* of the given labels are
                            returned. Also accepts 'union' and 'any'.
            'and'           Throats with all of the given labels are returned.
                            Also accepts 'intersection' and 'all'.
            'xor'           Throats with *only one* of the given labels are returned.
                            Also accepts 'exclusive_or'.
            'nor'           Throats with *none* of the given labels are returned.
                            Also accepts 'not' and 'none'.
            'nand'          Throats with *not all* of the given labels are
                            returned.
            'xnor'          Throats with *more than one* of the given labels are
                            returned.
            ==============  =====================================================

        asmask : bool
            If ``True`` then a boolean array of length Nt is returned with
            ``True`` values indicating the throats that satisfy the query.
        to_global : bool
            If ``True``, the returned indices will be indexed relative to the
            full domain.  This only has an effect when the calling object
            is a Subdomain.

        Returns
        -------
        A Numpy array containing throat indices filtered by the logic specified
        in ``mode``.

        See Also
        --------
        pores

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[3, 3, 3])
        >>> Ts = pn.throats()
        >>> Ts[0:5]  # Look at first 5 throat indices
        array([0, 1, 2, 3, 4])

        """
        ind = self._get_indices(element='throat', labels=labels, mode=mode)
        if to_global and hasattr(self, 'to_global'):
            ind = self.to_global(throats=ind)
            if asmask:
                ind = self._domain.to_mask(throats=ind)
        else:
            if asmask:
                ind = self.to_mask(throats=ind)
        return ind

    def filter_by_label(self, pores=[], throats=[], labels=None, mode='or'):
        r"""
        Returns which of the supplied pores (or throats) has the specified
        label(s)

        Parameters
        ----------
        pores, or throats : array_like
            List of pores or throats to be filtered
        labels : list of strings
            The labels to apply as a filter
        mode : str
            Controls how the filter is applied. The default value is
            'or'. Options include:

            ===========  =====================================================
            mode         meaning
            ===========  =====================================================
            'or'         Returns a list of the given locations where *any* of
                         the given labels exist. Also accepts 'union' and 'any'.
            'and'        Only locations where *all* the given labels are
                         found. Also accepts 'intersection' and 'all'.
            'xor'        Only locations where exactly *one* of the given
                         labels are found.
            'nor'        Only locations where *none* of the given labels are
                         found. Also accepts 'none' and 'not'
            'nand'       Only locations with *some but not all* of the given
                         labels are returned
            'xnor'       Only locations with *more than one* of the given
                         labels are returned
            ===========  =====================================================

        Returns
        -------
        A list of pores (or throats) that have been filtered according the
        given criteria.  The returned list is a subset of the received list of
        pores (or throats).

        See Also
        --------
        pores
        throats

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.filter_by_label(pores=[0, 1, 25, 32], labels='left')
        array([0, 1])
        >>> Ps = pn.pores(['top', 'bottom', 'back'], mode='or')
        >>> pn.filter_by_label(pores=Ps, labels=['top', 'back'],
        ...                    mode='and')
        array([ 24,  49,  74,  99, 124])
        """
        # Convert inputs to locations and element
        if (np.size(throats) > 0) and (np.size(pores) > 0):
            raise Exception('Can only filter either pores OR labels')
        if np.size(pores) > 0:
            element = 'pore'
            locations = self._parse_indices(pores)
        elif np.size(throats) > 0:
            element = 'throat'
            locations = self._parse_indices(throats)
        else:
            return(np.array([], dtype=int))
        labels = self._parse_labels(labels=labels, element=element)
        labels = [element+'.'+item.split('.')[-1] for item in labels]
        all_locs = self._get_indices(element=element, labels=labels, mode=mode)
        mask = self._tomask(indices=all_locs, element=element)
        ind = mask[locations]
        return locations[ind]

    def num_pores(self, labels='all', mode='or'):
        r"""
        Returns the number of pores of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The pore labels that should be included in the count.
            If not supplied, all pores are counted.
        labels : list of strings
            Label of pores to be returned
        mode : str, optional
            Specifies how the count should be performed. The options are:

            ==============  =====================================================
            mode            meaning
            ==============  =====================================================
            'or'            Pores with *one or more* of the given labels are
                            counted. Also accepts 'union' and 'any'.
            'and'           Pores with all of the given labels are returned.
                            Also accepts 'intersection' and 'all'.
            'xor'           Pores with *only one* of the given labels are
                            counted. Also accepts 'exclusive_or'.
            'nor'           Pores with *none* of the given labels are counted.
                            Also accepts 'not' and 'none'.
            'nand'          Pores with *not all* of the given labels are
                            counted.
            'xnor'          Pores with *more than one* of the given labels are
                            counted.
            ==============  =====================================================


        Returns
        -------
        Np : int
            Number of pores with the specified labels

        See Also
        --------
        num_throats
        count

        Notes
        -----
        Technically, *'nand'* and *'xnor'* should also count pores with *none*
        of the labels, however, to make the count more useful these are not
        included.

        Examples
        --------
        >>> import openpnm as op
        >>> pn = op.network.Cubic(shape=[5, 5, 5])
        >>> pn.num_pores()
        125
        >>> pn.num_pores(labels=['top'])
        25
        >>> pn.num_pores(labels=['top', 'front'], mode='or')
        45
        >>> pn.num_pores(labels=['top', 'front'], mode='xnor')
        5

        """
        # Count number of pores of specified type
        Ps = self._get_indices(labels=labels, mode=mode, element='pore')
        Np = np.shape(Ps)[0]
        return Np

    def num_throats(self, labels='all', mode='union'):
        r"""
        Return the number of throats of the specified labels

        Parameters
        ----------
        labels : list of strings, optional
            The throat labels that should be included in the count.
            If not supplied, all throats are counted.
        mode : str, optional
            Specifies how the count should be performed.  The options are:

            ==============  =====================================================
            mode            meaning
            ==============  =====================================================
            'or'            Throats with *one or more* of the given labels are
                            counted. Also accepts 'union' and 'any'.
            'and'           Throats with all of the given labels are returned.
                            Also accepts 'intersection' and 'all'.
            'xor'           Throats with *only one* of the given labels are
                            counted. Also accepts 'exclusive_or'.
            'nor'           Throats with *none* of the given labels are counted.
                            Also accepts 'not' and 'none'.
            'nand'          Throats with *not all* of the given labels are
                            counted.
            'xnor'          Throats with *more than one* of the given labels are
                            counted.
            ==============  =====================================================

        Returns
        -------
        Nt : int
            Number of throats with the specified labels

        See Also
        --------
        num_pores
        count

        Notes
        -----
        Technically, *'nand'* and *'xnor'* should also count throats with
        *none* of the labels, however, to make the count more useful these are
        not included.

        """
        # Count number of pores of specified type
        Ts = self._get_indices(labels=labels, mode=mode, element='throat')
        Nt = np.shape(Ts)[0]
        return Nt

    def props(self, *args, **kwargs):
        # Overload props on base to remove labels
        props = set(super().props(*args, **kwargs))
        labels = set(self.labels())
        props = props.difference(labels)
        return PrintableList(props)

    def __str__(self):
        s = super().__str__()
        # s = s.rpartition('\n')[0]
        horizontal_rule = '―' * 78
        lines = []
        lines.append(s)
        lines.append("{0:<5s} {1:<45s} {2:<10s}".format('#',
                                                        'Labels',
                                                        'Assigned Locations'))
        lines.append(horizontal_rule)
        labels = self.labels()
        labels.sort()
        fmt = "{0:<5d} {1:<45s} {2:<10d}"
        for i, item in enumerate(labels):
            prop = item
            if len(prop) > 35:
                prop = prop[0:32] + '...'
            if '._' not in prop:
                lines.append(fmt.format(i + 1, prop, np.sum(self[item])))
        lines.append(horizontal_rule)
        return '\n'.join(lines)
