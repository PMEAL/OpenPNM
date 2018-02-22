from openpnm.core import logging
from openpnm.algorithms import TransientTransport
logger = logging.getLogger(__name__)


class CrankNicholson(TransientTransport):
    r"""
    """

    def build_B(self):
        r"""
        """
        print('CrankNicholsonMixin: build_B')

    def update_A(self):
        r"""
        """
        print('CrankNicholsonMixin: update_A')
        self.build_B()
