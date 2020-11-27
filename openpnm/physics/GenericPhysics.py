import numpy as np
from openpnm.core import Subdomain, ModelsMixin
from openpnm.utils import Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericPhysics(Subdomain, ModelsMixin):
    r"""
    This generic class is meant as a starter for custom Physics objects

    It produces a blank object with no pore-scale models attached.  Users can
    add models from the ``models`` module (or create their own).

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this Physics should be attached

    phase : OpenPNM Phase object
        The Phase object to which this Physics applies

    geometry : OpenPNM Geometry object
        The Geometry object that defines the pores/throats where this Physics
        should be applied.

    name : str, optional
        A unique string name to identify the Physics object, typically same as
        instance name but can be anything.  If left blank, and name will be
        generated that includes the class name and an integer index.

    """

    def __init__(self, project=None, network=None, phase=None,
                 geometry=None, settings={}, **kwargs):

        # Define some default settings
        self.settings.update({'prefix': 'phys'})
        # Overwrite with user supplied settings, if any
        self.settings.update(settings)
        if phase is None:
            self.settings['freeze_models'] = True

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project

        super().__init__(project=project, **kwargs)

        network = self.project.network
        if network:
            if phase is None:
                logger.warning('No Phase provided, ' + self.name
                               + ' will not be associated with a phase')
            else:
                self.set_phase(phase=phase)
            if geometry is None:
                logger.warning('No Geometry provided, ' + self.name
                               + ' will not be associated with any locations')
            else:
                if (phase is None):
                    logger.warning('Cannot associate with a geometry unless '
                                   + 'a phase is also given')
                else:
                    self.set_geometry(geometry=geometry)

    def set_phase(self, phase=None, mode='swap'):
        r"""
        Sets the association between this physics and a phase.

        Parameters
        ----------
        phase : OpenPNM Phase object
            If mode is 'add' or 'swap', this must be specified so that
            associations can be recorded in the phase dictionary.  If the
            mode is 'drop', this is not needed since the existing association
            can be used to find it.
        mode : str
            Options are:

            'swap' - Associations will be made with the new phase, and
            the pore and throat locations from the current phase will be
            transferred to the new one.

            'drop' - Associations with the existing phase will be removed.

            'add' - If the physics does not presently have an associated
            phase, this will create associations, but no pore or throat
            locations will assigned.  This must be done using the
            ``set_geometry`` method.

        Notes
        -----
        In all cases the property data will be deleted since it will not
        be relevant to the new phase, so the ``regenerate_models`` method
        must be run.

        """
        if mode in ['add', 'swap']:
            if phase not in self.project:
                raise Exception(self.name + ' not in same project as given phase')
            try:
                old_phase = self.project.find_phase(self)
                phase['pore.'+self.name] = old_phase['pore.'+self.name]
                phase['throat.'+self.name] = old_phase['throat.'+self.name]
                old_phase.pop('pore.'+self.name, None)
                old_phase.pop('throat.'+self.name, None)
                self.clear()
            except Exception as e:
                logger.debug(e)
                phase['pore.'+self.name] = False
                phase['throat.'+self.name] = False
        elif mode in ['remove', 'drop']:
            self.update({'pore.all': np.array([], dtype=bool)})
            self.update({'throat.all': np.array([], dtype=bool)})
            phase = self.project.find_phase(self)
            phase.pop('pore.'+self.name, None)
            phase.pop('throat.'+self.name, None)
            self.clear()
        else:
            raise Exception('mode ' + mode + ' not understood')

    def set_geometry(self, geometry=None, mode='add'):
        r"""
        Sets the association between this physics and a geometry (i.e. a
        set of pores and throats that define a subdomain)

        Parameters
        ----------
        geometry : OpenPNM Geometry object
            The geometry defining the pores and throats to which this physics
            should be attached
        mode : str
            Options are:

            'swap' - Associations will be made with the new geometry, and
            the pore and throat locations from the current geometry will be
            transferred to the new one.

            'drop' - Associations with the current geometry will be removed.

            'add' - If the physics does not presently have an associated
            geometry, this will create associations.

        """
        if mode in ['add', 'swap']:
            if geometry not in self.project:
                raise Exception(self.name + ' not in same project as given geometry')
            try:
                old_geometry = self.project.find_geometry(self)
                Ps = self.network.pores(old_geometry.name)
                Ts = self.network.throats(old_geometry.name)
                self._set_locations(element='pore', indices=Ps, mode='drop')
                self._set_locations(element='throat', indices=Ts, mode='drop')
            except Exception as e:
                logger.debug(e)
            Ps = self.network.pores(geometry.name)
            Ts = self.network.throats(geometry.name)
            self._set_locations(element='pore', indices=Ps, mode='add')
            self._set_locations(element='throat', indices=Ts, mode='add')
        if mode in ['remove', 'drop']:
            phase = self.project.find_phase(self)
            phase['pore.'+self.name] = False
            phase['throat.'+self.name] = False
            self.update({'pore.all': np.array([], dtype=bool)})
            self.update({'throat.all': np.array([], dtype=bool)})
            self.clear()
