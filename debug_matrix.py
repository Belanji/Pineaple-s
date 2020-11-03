# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np

mat=np.loadtxt("inspect",skiprows=2)
l1max=5;

Npc=mat[0:l1max+3,0:l1max+3]
Nmc=mat[l1max+3:2*l1max+6,l1max+3:2*l1max+6]