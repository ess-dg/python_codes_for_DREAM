#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 11:17:37 2022

@author: irinastefanescu
"""
no_modulesM =  7 # 14  # no of Mantle modules
# Settings for endcap detectors
index_rot = 17  # index to define first module of endcap: index_rot x 12 degrees
no_modules = 13 # 23   # no of EndCap sectors
bwd_keep_high = 13  # no of last module of endcap bwd to keep
bwd_keep_low = 3  # no of first module of endcap bwd to keep
fwd_keep_high = 10  # no of last module of endcap fwd to keep
fwd_keep_low = 6  # no of first module of endcap bwd to keep

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
