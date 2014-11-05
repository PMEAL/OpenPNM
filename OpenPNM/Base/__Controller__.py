'''
###############################################################################
Controller:  Overall simulation controller class
###############################################################################
'''
import pprint as pp

class Controller(dict):
    r"""

    """

    def __init__(self):
        pass

    def __str__(self):
        for item in self.keys():
            print(item)
        return ''

    def add(self,obj):
        self[obj.name] = obj

    def find(self,name):
        pass

    def save(self,filename,filetype='pnm'):
        pass

    def load(self,filename):
        pass

    @classmethod
    def _add_logger(cls,**kwargs):
        import logging as _logging
        # set up logging to file - see previous section for more details
        _logging.basicConfig(level=_logging.ERROR,
                             format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                             datefmt='%m-%d %H:%M',
                             )

        if 'loggername' in kwargs.keys():
            cls._logger = _logging.getLogger(kwargs['loggername'])
        else:
            cls._logger = _logging.getLogger(cls.__class__.__name__)
        if 'loglevel' in kwargs.keys():
            cls._loglevel = kwargs['loglevel']


if __name__ == '__main__':
    sim = Controller()












