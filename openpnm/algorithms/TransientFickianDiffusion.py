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
            self.settings['phase'] = phase.name

    @docstr.dedent
    def setup(self, **kwargs):
        r"""
        Parameters
        ----------
        %(ReactiveTransport.setup.parameters)s

        ----

        The following settings are used to control the details of the transient
        simulation:

        %(TransientReactiveTransport.setup.other_parameters)s

        ----

        The following settings are used by the source term iterations:

        %(ReactiveTransport.setup.other_parameters)s

        %(ReactiveTransport.setup.notes)s

        ----

        The following settings are used to control the behavior of the solver:

        %(GenericTransport.setup.other_parameters)s
        """
