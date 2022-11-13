import numpy as np
from openpnm.utils import PrintableList


__all__ = [
    'ParserMixin',
    'LabelMixin',
]


class ParserMixin:

    def _parse_indices(self, indices):
        r"""
        This private method accepts a list of pores or throats and returns a
        properly structured Numpy array of indices.

        Parameters
        ----------
        indices : int or array_like
            This argument can accept numerous different data types including
            boolean masks, integers and arrays.

        Returns
        -------
        A Numpy array of indices.

        Notes
        -----
        This method should only be called by the method that is actually using
        the locations, to avoid calling it multiple times.

        """
        if indices is None:
            indices = np.array([], ndmin=1, dtype=int)
        locs = np.array(indices, ndmin=1)
        # If boolean array, convert to indices
        if locs.dtype == bool:
            if np.size(locs) == self.Np:
                locs = self.Ps[locs]
            elif np.size(locs) == self.Nt:
                locs = self.Ts[locs]
            else:
                raise Exception('Mask of locations must be either '
                                + 'Np nor Nt long')
        locs = locs.astype(dtype=int)
        return locs

    def _parse_element(self, element, single=False):
        r"""
        This private method is used to parse the keyword \'element\' in many
        of the above methods.

        Parameters
        ----------
        element : str or List[str]
            The element argument to check.  If is None is recieved, then a list
            containing both \'pore\' and \'throat\' is returned.
        single : bool (default is False)
            When set to True only a single element is allowed and it will also
            return a string containing the element.

        Returns
        -------
        When ``single`` is ``False`` (default) a list containing the element(s)
        is returned.  When ``single`` is ``True`` a bare string containing the
        element is returned.

        """
        if element is None:
            element = ['pore', 'throat']
        # Convert element to a list for subsequent processing
        if isinstance(element, str):
            element = [element]
        # Convert 'pore.prop' and 'throat.prop' into just 'pore' and 'throat'
        element = [item.split('.', 1)[0] for item in element]
        # Make sure all are lowercase
        element = [item.lower() for item in element]
        # Deal with an plurals
        element = [item.rsplit('s', maxsplit=1)[0] for item in element]
        for item in element:
            if item not in ['pore', 'throat']:
                raise Exception('All keys must start with either pore or throat')
        # Remove duplicates if any
        _ = [element.remove(L) for L in element if element.count(L) > 1]
        if single:
            if len(element) > 1:
                raise Exception('Both elements recieved when single element '
                                + 'allowed')
            element = element[0]
        return element

    def _parse_labels(self, labels, element):
        r"""
        This private method is used for converting \'labels\' to a proper
        format, including dealing with wildcards (\*).

        Parameters
        ----------
        labels : str or List[str]
            The label or list of labels to be parsed. Note that the \* can be
            used as a wildcard.

        Returns
        -------
        A list of label strings, with all wildcard matches included if
        applicable.

        """
        if labels is None:
            raise Exception('Labels cannot be None')
        if isinstance(labels, str):
            labels = [labels]
        # Parse the labels list
        parsed_labels = []
        for label in labels:
            # Remove element from label, if present
            if element in label:
                label = label.split('.', 1)[-1]
            # Deal with wildcards
            if '*' in label:
                Ls = [L.split('.', 1)[-1] for L in self.labels(element=element)]
                if label.startswith('*'):
                    temp = [L for L in Ls if L.endswith(label.strip('*'))]
                if label.endswith('*'):
                    temp = [L for L in Ls if L.startswith(label.strip('*'))]
                temp = [element+'.'+L for L in temp]
            elif element+'.'+label in self.keys():
                temp = [element+'.'+label]
            else:
                temp = [element+'.'+label]
            parsed_labels.extend(temp)
            # Remove duplicates if any
            _ = [parsed_labels.remove(L) for L in parsed_labels
                 if parsed_labels.count(L) > 1]
        return parsed_labels

    def _parse_mode(self, mode, allowed=None, single=False):
        r"""
        This private method is for checking the \'mode\' used in the calling
        method.

        Parameters
        ----------
        mode : str or List[str]
            The mode(s) to be parsed
        allowed : List[str]
            A list containing the allowed modes.  This list is defined by the
            calling method.  If any of the received modes are not in the
            allowed list an exception is raised.
        single : bool (default is False)
            Indicates if only a single mode is allowed.  If this argument is
            True than a string is returned rather than a list of strings, which
            makes it easier to work with in the caller method.

        Returns
        -------
        A list containing the received modes as strings, checked to ensure they
        are all within the allowed set (if provoided).  Also, if the ``single``
        argument was True, then a string is returned.

        """
        if isinstance(mode, str):
            mode = [mode]
        for item in mode:
            if (allowed is not None) and (item not in allowed):
                raise Exception('\'mode\' must be one of the following: '
                                + allowed.__str__())
        # Remove duplicates, if any
        _ = [mode.remove(L) for L in mode if mode.count(L) > 1]
        if single:
            if len(mode) > 1:
                raise Exception('Multiple modes received when only one mode '
                                + 'is allowed by this method')
            mode = mode[0]
        return mode

    def _parse_prop(self, propname, element):
        element = self._parse_element(element, single=True)
        if propname.split('.', 1)[0] in ['pore', 'throat']:
            propname = propname.split('.', 1)[-1]
        return element + '.' + propname


class LabelMixin:
    """r
    This mixin adds functionality to the Base2 class so that boolean arrays
    are treated as labels
    """

    def _get_labels(self, element, locations, mode):
        r"""
        This is the actual label getter method, but it should not be called
        directly.  Use ``labels`` instead.
        """
        # Parse inputs
        locations = self._parse_indices(locations)
        element = self._parse_element(element=element)
        # Collect list of all pore OR throat labels
        labels = [i for i in self.keys(mode='labels') if i.split('.', 1)[0] in element]
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

            ==============  ===================================================
            mode            meaning
            ==============  ===================================================
            'or'            Returns the labels that are assigned to *any* of
                            the given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the
                            given locations
            'xnor'          Labels that are present on *more than one* of the
                            given locations.
            ==============  ===================================================

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

        """
        # Short-circuit query when no pores or throats are given
        if (np.size(pores) == 0) and (np.size(throats) == 0):
            if element is None:
                element = ['pore', 'throat']
            if isinstance(element, str):
                element = [element]
            labels = PrintableList()
            for k, v in self.items():
                el, prop = k.split('.', 1)
                if (el in element) and (v.dtype == bool) and not prop.startswith('_'):
                    labels.append(k)
        elif (np.size(pores) > 0) and (np.size(throats) > 0):
            raise Exception('Cannot perform label query on pores and '
                            + 'throats simultaneously')
        elif np.size(pores) > 0:
            labels = self._get_labels(element='pore', locations=pores,
                                      mode=mode)
        elif np.size(throats) > 0:
            labels = self._get_labels(element='throat', locations=throats,
                                      mode=mode)
        return sorted(labels)

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

        if label.split('.', 1)[0] in ['pore', 'throat']:
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

    def _get_indices(self, element, labels, mode='or'):
        r"""
        This is the actual method for getting indices, but should not be called
        directly.  Use ``pores`` or ``throats`` instead.
        """
        # Parse and validate all input values.
        element = self._parse_element(element, single=True)
        labels = self._parse_labels(labels=labels, element=element)

        # Begin computing label array
        if mode in ['or', 'any', 'union']:
            union = np.zeros([self._count(element), ], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                union = union + self[element+'.'+item.split('.', 1)[-1]]
            ind = union
        elif mode in ['and', 'all', 'intersection']:
            intersect = np.ones([self._count(element), ], dtype=bool)
            for item in labels:  # Iterate over labels and collect all indices
                intersect = intersect*self[element+'.'+item.split('.', 1)[-1]]
            ind = intersect
        elif mode in ['xor', 'exclusive_or']:
            xor = np.zeros([self._count(element), ], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.', 1)[-1]]
                xor = xor + np.int8(info)
            ind = (xor == 1)
        elif mode in ['nor', 'not', 'none']:
            nor = np.zeros([self._count(element), ], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.', 1)[-1]]
                nor = nor + np.int8(info)
            ind = (nor == 0)
        elif mode in ['nand']:
            nand = np.zeros([self._count(element), ], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.', 1)[-1]]
                nand = nand + np.int8(info)
            ind = (nand < len(labels)) * (nand > 0)
        elif mode in ['xnor', 'nxor']:
            xnor = np.zeros([self._count(element), ], dtype=int)
            for item in labels:  # Iterate over labels and collect all indices
                info = self[element+'.'+item.split('.', 1)[-1]]
                xnor = xnor + np.int8(info)
            ind = (xnor > 1)
        else:
            raise Exception('Unsupported mode: '+mode)
        # Extract indices from boolean mask
        ind = np.where(ind)[0]
        ind = ind.astype(dtype=int)
        return ind

    def pores(self, labels=None, mode='or', asmask=False):
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

            ==============  ===================================================
            mode            meaning
            ==============  ===================================================
            'or'            Returns the labels that are assigned to *any* of
                            the given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the
                            given locations
            'xnor'          Labels that are present on *more than one* of the
                            given locations.
            ==============  ===================================================

        asmask : bool
            If ``True`` then a boolean array of length Np is returned with
            ``True`` values indicating the pores that satisfy the query.

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

        """
        if labels is None:
            labels = self.name
        ind = self._get_indices(element='pore', labels=labels, mode=mode)
        if asmask:
            ind = self.to_mask(pores=ind)
        return ind

    def throats(self, labels=None, mode='or', asmask=False):
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

            ==============  ===================================================
            mode            meaning
            ==============  ===================================================
            'or'            Returns the labels that are assigned to *any* of
                            the given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the
                            given locations
            'xnor'          Labels that are present on *more than one* of the
                            given locations.
            ==============  ===================================================

        asmask : bool
            If ``True`` then a boolean array of length Nt is returned with
            ``True`` values indicating the throats that satisfy the query.

        Returns
        -------
        A Numpy array containing throat indices filtered by the logic specified
        in ``mode``.

        See Also
        --------
        pores

        """
        if labels is None:
            labels = self.name
        ind = self._get_indices(element='throat', labels=labels, mode=mode)
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

            ==============  ===================================================
            mode            meaning
            ==============  ===================================================
            'or'            Returns the labels that are assigned to *any* of
                            the given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the
                            given locations
            'xnor'          Labels that are present on *more than one* of the
                            given locations.
            ==============  ===================================================

        Returns
        -------
        A list of pores (or throats) that have been filtered according the
        given criteria.  The returned list is a subset of the received list of
        pores (or throats).

        See Also
        --------
        pores
        throats

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
            return np.array([], dtype=int)
        labels = self._parse_labels(labels=labels, element=element)
        labels = [element+'.'+item.split('.', 1)[-1] for item in labels]
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

            ==============  ===================================================
            mode            meaning
            ==============  ===================================================
            'or'            Returns the labels that are assigned to *any* of
                            the given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the
                            given locations
            'xnor'          Labels that are present on *more than one* of the
                            given locations.
            ==============  ===================================================

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

            ==============  ===================================================
            mode            meaning
            ==============  ===================================================
            'or'            Returns the labels that are assigned to *any* of
                            the given locations. Also accepts 'union' and 'any'
            'and'           Labels that are present on all the given locations.
                            also accepts 'intersection' and 'all'
            'xor'           Labels that are present on *only one*
                            of the given locations.Also accepts 'exclusive_or'
            'nor'           Labels that are *not* present on any of
                            the given locations. Also accepts 'not' and 'none'
            'nand'          Labels that are present on *all but one* of the
                            given locations
            'xnor'          Labels that are present on *more than one* of the
                            given locations.
            ==============  ===================================================

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
