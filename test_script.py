# -*- coding: utf-8 -*-

import os
import sys

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
    