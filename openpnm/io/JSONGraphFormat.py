from openpnm.core import logging
from openpnm.io import GenericIO
logger = logging.getLogger(__name__)


class JSONGraphFormat(GenericIO):
    r"""
    Class for reading and writing OpenPNM networks to JSON Graph Format (JGF).

    Notes
    -----
    The JGF standard must contain data formatted according to http://jsongraphformat.info and
    enforced by JSON schema validation.
    """

    @classmethod
    def save(cls, network, phases=[], filename=''):
        r"""
        Write the wetwork to disk as a JGF file.

        Parameters
        ----------
        network : OpenPNM Network Object

        filename : string
            Desired file name, defaults to network name if not given

        phases : list of phase objects ([])
            Phases that have properties we want to write to file

        """

    @classmethod
    def load(cls, filename, project=None):
        r"""
        Loads the JGF file onto the given project.

        Parameters
        ----------
        filename : string (optional)
            The name of the file containing the data to import.  The formatting
            of this file is outlined below.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project object is supplied then one will be created and
            returned.

        Returns
        -------
        If no project object is supplied then one will be created and returned.

        """
        return project
