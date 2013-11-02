# -*- coding: utf-8 -*-

import os

r"""
This script is to be run before any commit to the development branch. 

All testing scripts should be housed at the base of the relevant .py files
"""

#move to the OpenPNM directory
os.chdir('OpenPNM')
#go into each folder present
OpenPNMfiles = os.listdir('.')
for itemname in OpenPNMfiles:
    if '.' not in itemname:
    
        os.chdir(itemname)
        curdir = os.listdir('.')
        for filename in curdir:
            if filename[-3:]=='.py':
                os.system('python '+filename)
        os.chdir('..')
    