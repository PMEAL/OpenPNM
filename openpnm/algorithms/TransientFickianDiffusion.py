from openpnm.algorithms import TransientReactiveTransport, FickianDiffusion
from openpnm.utils import logging, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


class TransientFickianDiffusion(TransientReactiveTransport, FickianDiffusion):
    r"""
    A subclass of GenericTransport to simulate diffusion.

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            't_initial': None,
                                            't_final': None,
                                            't_step': None,
                                            't_output': None,
                                            't_tolerance': None,
                                            't_scheme': ''},
                           'set_IC':       {'values': None},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    @docstr.dedent
    def setup(self, **kwargs):
        r"""
        Other Parameters
        ----------------
        %(TransientReactiveTransport.setup.other_parameters)s

        """
        super().update(**kwargs)
