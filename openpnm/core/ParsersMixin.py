import scipy as sp
from openpnm.core import logging
logger = logging.getLogger()


class ParsersMixin():

    def _parse_indices(self, indices):
        r"""
        This private method accepts a list of pores or throats and returns a
        properly structured Numpy array of indices.

        Parameters
        ----------
        indices : multiple options
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
            indices = sp.array([], ndmin=1, dtype=int)
        locs = sp.array(indices, ndmin=1)
        if sp.issubdtype(locs.dtype, sp.number) or \
           sp.issubdtype(locs.dtype, sp.bool_):
            pass
        else:
            raise Exception('Invalid data type received as indices, ' +
                            ' must be boolean or numeric')
        if locs.dtype == bool:
            if sp.size(locs) == self.Np:
                locs = self.Ps[locs]
            elif sp.size(locs) == self.Nt:
                locs = self.Ts[locs]
            else:
                raise Exception('Boolean list of locations must be either ' +
                                'Np nor Nt long')
        locs = locs.astype(dtype=int)
        return locs

    def _parse_element(self, element, single=False):
        r"""
        This private method is used to parse the keyword \'element\' in many
        of the above methods.

        Parameters
        ----------
        element : string or list of strings
            The element argument to check.  If is None is recieved, then a list
            containing both \'pore\' and \'throat\' is returned.

        single : boolean (default is False)
            When set to True only a single element is allowed and it will also
            return a string containing the element.

        Returns
        -------
        When ``single`` is False (default) a list contain the element(s) is
        returned.  When ``single`` is True a bare string containing the element
        is returned.
        """
        if element is None:
            element = ['pore', 'throat']
        # Convert element to a list for subsequent processing
        if type(element) is str:
            element = [element]
        # Convert 'pore.prop' and 'throat.prop' into just 'pore' and 'throat'
        element = [item.split('.')[0] for item in element]
        # Make sure all are lowercase
        element = [item.lower() for item in element]
        # Deal with an plurals
        element = [item.rsplit('s', maxsplit=1)[0] for item in element]
        for item in element:
            if item not in ['pore', 'throat']:
                raise Exception('Invalid element received: '+item)
        # Remove duplicates if any
        [element.remove(L) for L in element if element.count(L) > 1]
        if single:
            if len(element) > 1:
                raise Exception('Both elements recieved when single element ' +
                                'allowed')
            else:
                element = element[0]
        return element

    def _parse_labels(self, labels, element):
        r"""
        This private method is used for converting \'labels\' to a proper
        format, including dealing with wildcards (\*).

        Parameters
        ----------
        labels : string or list of strings
            The label or list of labels to be parsed. Note that the \* can be
            used as a wildcard.

        Returns
        -------
        A list of label strings, with all wildcard matches included if
        applicable.
        """
        if labels is None:
            raise Exception('Labels cannot be None')
        if type(labels) is str:
            labels = [labels]
        # Parse the labels list
        parsed_labels = []
        for label in labels:
            # Remove element from label, if present
            if element in label:
                label = label.split('.')[-1]
            # Deal with wildcards
            if '*' in label:
                Ls = [L.split('.')[-1] for L in self.labels(element=element)]
                if label.startswith('*'):
                    temp = [L for L in Ls if L.endswith(label.strip('*'))]
                if label.endswith('*'):
                    temp = [L for L in Ls if L.startswith(label.strip('*'))]
                temp = [element+'.'+L for L in temp]
            elif element+'.'+label in self.keys():
                temp = [element+'.'+label]
            else:
                # TODO: The following Error should/could be raised but it
                # breaks the net-geom and phase-phys look-up logic
                # raise KeyError('\''+element+'.'+label+'\''+' not found')
                logger.warning('\''+element+'.'+label+'\''+' not found')
                temp = [element+'.'+label]
            parsed_labels.extend(temp)
            # Remove duplicates if any
            [parsed_labels.remove(L) for L in parsed_labels
             if parsed_labels.count(L) > 1]
        return parsed_labels

    def _parse_mode(self, mode, allowed=None, single=False):
        r"""
        This private method is for checking the \'mode\' used in the calling
        method.

        Parameters
        ----------
        mode : string or list of strings
            The mode(s) to be parsed

        allowed : list of strings
            A list containing the allowed modes.  This list is defined by the
            calling method.  If any of the received modes are not in the
            allowed list an exception is raised.

        single : boolean (default is False)
            Indicates if only a single mode is allowed.  If this argument is
            True than a string is returned rather than a list of strings, which
            makes it easier to work with in the caller method.

        Returns
        -------
        A list containing the received modes as strings, checked to ensure they
        are all within the allowed set (if provoided).  Also, if the ``single``
        argument was True, then a string is returned.
        """
        if type(mode) is str:
            mode = [mode]
        for item in mode:
            if (allowed is not None) and (item not in allowed):
                raise Exception('\'mode\' must be one of the following: ' +
                                allowed.__str__())
        # Remove duplicates, if any
        [mode.remove(L) for L in mode if mode.count(L) > 1]
        if single:
            if len(mode) > 1:
                raise Exception('Multiple modes received when only one mode ' +
                                'allowed')
            else:
                mode = mode[0]
        return mode

    def _parse_prop(self, propname, element):
        r"""
        """
        return element + '.' + propname.split('.')[-1]