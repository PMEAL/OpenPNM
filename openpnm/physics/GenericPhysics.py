import numpy as np
from openpnm.core import Subdomain, ModelsMixin, ParamMixin
from openpnm.utils import Workspace, logging
logger = logging.getLogger(__name__)
ws = Workspace()


class GenericPhysics(ParamMixin, Subdomain, ModelsMixin):
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

    def __init__(self, phase=None, geometry=None, pores=None, throats=None,
                 settings={}, **kwargs):
        self.settings.update({'prefix': 'phys'})  # Define some default settings
        self.settings.update(settings)  # Overwrite with user supplied settings
        super().__init__(**kwargs)

        network = self.project.network
        if network:
            if phase is not None:
                self.set_phase(phase=phase, mode='add')
            if geometry is None:
                if (pores is None) and (throats is None):
                    logger.warning('No Geometry provided, ' + self.name
                                   + ' will not be associated with any locations')
                else:
                    self.set_locations(pores=pores, throats=throats, mode='add')
            else:
                if phase is None:
                    logger.warning('Cannot associate with a geometry unless '
                                   + 'a phase is also given')
                else:
                    self.set_geometry(geometry=geometry)

    def _set_phase(self, phase):
        if phase is None:
            self._del_phase()
            return
        try:
            self.set_phase(phase=phase, mode='move')
        except Exception:
            self.set_phase(phase=phase, mode='add')

    def _get_phase(self):
        return self.project.find_phase(self)

    def _del_phase(self):
        self.set_phase(phase=self.phase, mode='drop')

    phase = property(fget=_get_phase, fset=_set_phase, fdel=_del_phase)

    def set_phase(self, phase=None, mode='add'):
        r"""
        Sets the association between this physics and a phase.

        Parameters
        ----------
        phase : OpenPNM Phase object
            If mode is 'add' or 'move', this must be specified so that
            associations can be recorded in the phase dictionary.  If the
            mode is 'drop', this is not needed since the existing association
            can be used to find it.
        mode : str
            Options are:

            * 'add' (default)
                If the physics does not presently have an associated phase,
                this will create associations, but no pore or throat
                locations will assigned.  This must be done using the
                ``set_geometry`` method.
            * 'drop'
                Associations with the existing phase will be removed.
            * 'move'
                Associations will be made with the new phase, and the pore
                and throat locations from the current phase will be
                transferred to the new one.

        Notes
        -----
        In all cases the property data will be deleted since it will not
        be relevant to the new phase, so the ``regenerate_models`` method
        must be run.

        """
        self._parse_mode(mode=mode, allowed=['add', 'drop', 'move'])
        if (phase is not None) and (phase not in self.project):
            raise Exception(self.name + ' not in same project as given phase')

        if mode in ['move']:
            old_phase = self.project.find_phase(self)
            # Transfer locs from old to new phase
            phase['pore.'+self.name] = old_phase.pop('pore.'+self.name, False)
            phase['throat.'+self.name] = old_phase.pop('throat.'+self.name, False)
            self.clear()
        elif mode in ['add']:
            temp = None
            try:
                temp = self.phase
            except Exception:
                phase['pore.'+self.name] = False
                phase['throat.'+self.name] = False
                return
            if temp is not None:
                raise Exception(f"{self.name} is already associated with " +
                f"{self.phase.name}. Use mode 'move' instead.")
        elif mode in ['drop']:
            self.update({'pore.all': np.array([], dtype=bool)})
            self.update({'throat.all': np.array([], dtype=bool)})
            self.phase.pop('pore.' + self.name, None)
            self.phase.pop('throat.' + self.name, None)
            self.clear()

    def _set_geo(self, geo):
        if geo is None:
            self._del_geo()
            return
        self.set_geometry(geo, mode='add')

    def _get_geo(self):
        return self.project.find_geometry(physics=self)

    def _del_geo(self):
        self.set_geometry(mode='drop')

    geometry = property(fset=_set_geo, fget=_get_geo, fdel=_del_geo)

    def set_geometry(self, geometry=None, mode='add'):
        r"""
        Sets the association between this physics and a geometry

        This association is done by setting the pores and throats that define
        the Subdomain to match.

        Parameters
        ----------
        geometry : OpenPNM Geometry object
            The geometry defining the pores and throats to which this physics
            should be attached
        mode : str
            Controls how the assignment is done. Options are:

            * 'add' (default)
                If the physics does not presently have an associated
                geometry, this will create associations, otherwise an
                Exception is raised
            * 'drop'
                Associations with the current geometry will be removed
            * 'move'
                Associations will be made with the provided geometry, and
                the pore and throat locations from the current geometry will
                be transferred to the new one.

        See Also
        --------
        set_locations

        """
        self._parse_mode(mode=mode, allowed=['add', 'move', 'drop'])
        try:
            phase = self.project.find_phase(self)
        except Exception:
            raise Exception("Cannot set locations for a physics object" +
                            " that is not associated with a phase")

        if (geometry is not None) and (geometry not in self.project):
            raise Exception(self.name + ' not in same project as given geometry')

        if mode in ['add']:
            Ps = self.network.pores(geometry.name)
            Ts = self.network.throats(geometry.name)
            self.set_locations(pores=Ps, throats=Ts, mode='add')
            phase.set_label(label=self.name, pores=Ps, throats=Ts, mode='add')
        elif mode in ['drop']:
            phase.set_label(label=self.name, mode='clear')
            self.update({'pore.all': np.array([], dtype=bool)})
            self.update({'throat.all': np.array([], dtype=bool)})
            self.clear()
        elif mode in ['move']:
            phys = self.project.find_physics(phase=phase, geometry=geometry)
            phys.set_geometry(mode='drop')
            self.set_geometry(mode='drop')
            self.set_geometry(geometry=geometry, mode='add')
