#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 11:17:37 2022

@author: irinastefanescu
"""


def initialize():
    
    global no_modules
    global no_modulesM
    no_modules = 1   # no of EndCap sectors
    no_modulesM = 10  # no of Mantle modules


def match2geantvars(args):
    """

    """
    dy1 = args[0]
    dx1 = args[1]
    dx2 = args[2]
    dy2 = args[3]
    dx3 = args[4]
    dx4 = args[5]
    k1 = (dx1 - dx2) / (dy1 + dy1)
    k2 = (dx3 - dx4) / (dy2 + dy2)
    k = (k1 + k2) / 2
    m1 = (dx1 + dx2) / 2
    m2 = (dx3 + dx4) / 2
    dx1 = m1 + k*dy1
    dx2 = m1 - k*dy1
    dx3 = m2 + k*dy2
    dx4 = m2 - k*dy2
    return dy1, dx1, dx2, dy2, dx3, dx4
