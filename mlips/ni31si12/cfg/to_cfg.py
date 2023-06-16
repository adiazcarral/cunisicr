#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 15:41:31 2022

@author: angel
"""

from mlip import *

filename = 'POSCAR'

convert_to_db(filename, format='vasp', db_path="cfgs.db")
db_to_cfg('ni31si12.cfg', db_path='cfgs.db', yaml_path='system.yaml')
