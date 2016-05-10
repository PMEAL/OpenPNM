"""
###############################################################################
Controller:  Overall controller class
###############################################################################
"""
from OpenPNM.Base import Workspace
        for item in self._comments.values():
            if 'Using OpenPNM' in item:
                version = item.lstrip('Using OpenPNM ')
                if version < OpenPNM.__version__:
                    logger.warning('File was created with an earlier version ' +
                                   'OpenPNM: \n' +
                                   '--> File saved with version: ' +
                                   str(version) +
                                   '\n' +
                                   '--> Current version: ' +
                                   str(OpenPNM.__version__))
            if logger.level > 20:
                logger.error('To view comments set the loglevel to 20 or less')
                logger.info(key + ': ' + self._comments[key])


class Controller(Workspace):
    pass
