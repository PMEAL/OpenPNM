"""
###############################################################################
Controller:  A dummy class support backwards compatibility
###############################################################################
"""
import OpenPNM
from OpenPNM.Base import logging
logger = logging.getLogger()


class Controller(OpenPNM.Base.Workspace):
    # This is a dummy class to provide the functionality of the Workspace
    # class # under the name Controller to provide backwards compatibility.
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
