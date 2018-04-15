from openpnm.algorithms import GenericAlgorithm
from openpnm.core import logging
logger = logging.getLogger(__name__)


class GenericPercolation(GenericAlgorithm):

    def set_inlets(self, pores=[], overwrite=False):
        r"""
        Set the locations from which the invader enters the network

        Parameters
        ----------
        pores : array_like
            Locations that are initially filled with invader, from which
            clusters grow and invade into the network

        overwrite : boolean
            If ``True`` then all existing inlet locations will be removed and
            then the supplied locations will be added.  If ``False`` (default),
            then supplied locations are added to any already existing inlet
            locations.


        """
        Ps = self._parse_indices(pores)
        if sum(self['pore.outlets'][Ps]) > 0:
            raise Exception('Some inlets are already defined as outlets')
        if overwrite:
            self['pore.inlets'] = False
        self['pore.inlets'][Ps] = True

    def set_outlets(self, pores=[], overwrite=False):
        r"""
        Set the locations through which defender exits the network.
        This is only necessary if 'trapping' was set to True when ``setup``
        was called.

        Parameters
        ----------
        pores : array_like
            Locations where the defender can exit the network.  Any defender
            that does not have access to these sites will be trapped.

        overwrite : boolean
            If ``True`` then all existing outlet locations will be removed and
            then the supplied locations will be added.  If ``False`` (default),
            then supplied locations are added to any already existing outlet
            locations.

        """
        if self.settings['trapping'] is False:
            logger.warning('Setting outlets is meaningless unless trapping ' +
                           'was set to True during setup')
        Ps = self._parse_indices(pores)
        if sum(self['pore.inlets'][Ps]) > 0:
            raise Exception('Some outlets are already defined as inlets')
        if overwrite:
            self['pore.outlets'] = False
        self['pore.outlets'][Ps] = True

    def set_residual(self, pores=[], throats=[], overwrite=False):
        r"""
        Specify locations of any residual invader.  These locations are set
        to invaded at the start of the simulation.

        Parameters
        ----------
        pores : array_like
            The pores locations that are to be filled with invader at the
            beginning of the simulation.

        throats : array_like
            The throat locations that are to be filled with invader at the
            beginning of the simulation.

        overwrite : boolean
            If ``True`` then all existing inlet locations will be removed and
            then the supplied locations will be added.  If ``False``, then
            supplied locations are added to any already existing locations.

        """
        Ps = self._parse_indices(pores)
        if overwrite:
            self['pore.residual'] = False
        self['pore.residual'][Ps] = True
        Ts = self._parse_indices(throats)
        if overwrite:
            self['throat.residual'] = False
        self['throat.residual'][Ts] = True
