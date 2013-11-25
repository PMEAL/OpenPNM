# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 22:09:10 2013

@author: Jeff
"""
import logging
try: print('Logger already instantiated, named: ', logger.name)
except:
    # create logger
    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(name)s: %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)

class foo:
    def __init__(self):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.test()
    def test(self):
        self._logger.debug('test4_debug')
        self._logger.info('test4_info')
        self._logger.warning('test4_warning')
        self._logger.error('test4_error')       
        self._logger.critical('test4_critical')
        
class spam:
    def __init__(self):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.test()
    def test(self):
        self._logger.debug('test5_debug')
        self._logger.info('test5_info')
        self._logger.warning('test5_warning')
        self._logger.error('test5_error')
        self._logger.critical('test5_critical')        
if __name__ =="__main__":
    y = foo()
    x = spam()