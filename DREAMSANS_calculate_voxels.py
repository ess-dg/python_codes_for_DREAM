#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sunday, January 28, 2024

@author: celinedurniak
"""

import numpy as np

np.set_printoptions(precision=4)

"""
*******The SANS for DREAM ******

 Code adapted from implementation for High-Resolution detector

data is in mm

"""

"""**************User parameters *****************"""

nSANS = 8             # no of segments per module
no_modulesSANS = 9    # no of modules in the frame
no_sectorsSANS = 4    # no of sectors

"""***********************************************"""

hSANS = 170       # segment height front, sample side
sensSANS = 313    # segment depth
wSANS = 20        # segment width front, sample side

margin = 10      # empty space inside the segment front-back and left-right
marginb = 1.5    # empty space inside the segment top-bottom
shieldz = 25     # length of the shielding block at the back of the segment
mod_dist = 2     # distance between neighboring segments in order to avoid overlap

AlthickH = 0.3      # thickness Al cathode material and housing
BthickSANS = 0.0013   # thickness Boron coating

# start position for placing the modules in the frame, integer number
# to multiply with 12*deg
init_angleSANS = 0     # init angle in deg

tilt_theta = -10     # tilt_angle in deg
tilt_thetaSANS = -1    # tilt_angle in deg
tilt_phiSANS = -1      # inclination angle module in deg

n_wiresSANS = 16     # no of wires
n_stripsSANS = 32    # no of strips

dist_detSANS = 2500   # distance from the sample in mm

# start calculations 

# Al housing (box)

AlySANS = (hSANS + 1) / 2   # 1/2 height Al housing at entrance window SANS
AlzSANS = (sensSANS + margin) / 2    # depth housing SANS
AlxSANS = wSANS / 2   # 1/2 width segment front in SANS

"""
 the next lines are for calculating the dimensions of
 the inner volume that will be filled with gas voxels
 the principle is similar to the Russian Matryoschka
 dolls (nested dolls)
 From the Al box representing the segment
 housing subtract the box
 representing the Boron coating on the side walls
 and the remaining volume is the gas volume,
 therefore the "G" letter used to name
 the variables relevant to the gas voxels
"""

# Boron box

BySANS = AlySANS - AlthickH
BzSANS = AlzSANS - AlthickH
BxSANS = AlxSANS - AlthickH

# Gas volume available to the gas voxels

GySANS = BySANS
GzSANS = BzSANS - BthickSANS
GxSANS = BxSANS - BthickSANS

CathSubstrYSANS = GySANS - marginb / 2
CathSubstrZSANS = BzSANS - margin / 2

CathConvYSANS = CathSubstrYSANS
CathConvZSANS = CathSubstrZSANS
CathConvXSANS = BthickSANS / 2
CathSubstrX1 = AlthickH / 2

# calculate the dimensions of the gas voxels

xxSANS = GxSANS - (2*CathConvXSANS + AlthickH)

izzSANS = sensSANS / n_stripsSANS  # strip pitch SANS in mm, all equal

GLzSANS = izzSANS / 2
GLySANS = hSANS / 2 / n_wiresSANS
GLxSANS = xxSANS / 2

# calculate the centers of the voxels 
shp = (n_wiresSANS, n_stripsSANS)

voxelXXhr = np.zeros(shp)
voxelXXchr = np.zeros(shp)
voxelYYhr = np.zeros(shp)
voxelZZhr = np.zeros(shp)

for strip in range(n_stripsSANS):
    for wire in range(n_wiresSANS):
        # fill the voxels from the top to bottom and from front to back
        voxelYYhr[wire, strip] = (n_wiresSANS - 2 * wire - 1) * GLySANS
        voxelZZhr[wire, strip] = (n_stripsSANS//2 - strip - 1) * izzSANS + margin/2
        voxelXXhr[wire, strip] = AlthickH
        voxelXXchr[wire, strip] = GLxSANS + 2 * CathConvXSANS + CathSubstrX1

# calculate the segment positions in the detector frame

segY = np.zeros(nSANS)

segZ = - dist_detSANS - AlzSANS
segX = 0

for seg_no in range(nSANS):
    segY[seg_no] = 80 + (2 * AlxSANS + 0.4) * (seg_no + 1)
    
# calculate the module positions in the detector frame

pos_mod = hSANS + mod_dist  # allow for 2 mm space between the modules
# Edit Feb 2024: to match the new design from CDT.
# The translation along the x-axis for Sector 1 has to be reversed

m_hrx = np.array([0, 0, 0, pos_mod, pos_mod, pos_mod, 2*pos_mod, 2*pos_mod, 3*pos_mod])

m_hry = np.array([0, pos_mod, 2*pos_mod, 0, pos_mod, 2*pos_mod, 0, pos_mod, 0])

m_hrz = 0

shr_angle = np.array([0, -90, 180, 90])

""" calculate the lookup table  """

shp = (no_sectorsSANS, 100 * no_modulesSANS + nSANS, n_wiresSANS, n_stripsSANS)

# voxel positions after placing the segment in the frame
sxhr_z = np.zeros(shp)
syhr_z = np.zeros(shp)
szhr_z = np.zeros(shp)

sxhr_zx = np.zeros(shp)
syhr_zx = np.zeros(shp)
szhr_zx = np.zeros(shp)

sxhr_zxy = np.zeros(shp)
syhr_zxy = np.zeros(shp)
szhr_zxy = np.zeros(shp)

# voxel positions after placing the module in the frame
mhrx = np.zeros(shp)
mhry = np.zeros(shp)
mhrz = np.zeros(shp)

secx = np.zeros(shp)
secy = np.zeros(shp)
secz = np.zeros(shp)

with open('SANS_temp.txt', "a") as f:

    # voxels in the 'top counter'
    for sec in range(no_sectorsSANS):

        secZ_s = np.sin(np.deg2rad(shr_angle[sec]))
        secZ_c = np.cos(np.deg2rad(shr_angle[sec]))

        for md in range(no_modulesSANS):

            for segment in range(nSANS):

                angZ = 90
                angX = tilt_thetaSANS
                angY = tilt_phiSANS

                segZ_s = np.sin(np.deg2rad(angZ))
                segZ_c = np.cos(np.deg2rad(angZ))
                segY_s = np.sin(np.deg2rad(angY))
                segY_c = np.cos(np.deg2rad(angY))
                segX_s = np.sin(np.deg2rad(angX))
                segX_c = np.cos(np.deg2rad(angX))

                md_segt_id = 100 * md + segment

                for strip in range(n_stripsSANS):
                    for wire in range(n_wiresSANS):
                        # rotation of each segment of the module by angZ
                        # around the Z-axis
                        sxhr_z[sec, md_segt_id, wire, strip] = -voxelYYhr[wire, strip] * segZ_s + \
                                                               (voxelXXhr[wire, strip] +
                                                                voxelXXchr[wire, strip]) * segZ_c

                        syhr_z[sec, md_segt_id, wire, strip] = \
                            voxelYYhr[wire, strip] * segZ_c + (voxelXXhr[wire, strip] +
                                                               voxelXXchr[wire, strip]) * segZ_s

                        szhr_z[sec, md_segt_id, wire, strip] = voxelZZhr[wire, strip]

                        # rotation of the segment around the X axis
                        sxhr_zx[sec, md_segt_id, wire, strip] = \
                            sxhr_z[sec, md_segt_id, wire, strip] * segX_c + \
                            szhr_z[sec, md_segt_id, wire, strip] * segX_s

                        syhr_zx[sec, md_segt_id, wire, strip] = syhr_z[sec, md_segt_id, wire, strip]

                        szhr_zx[sec, md_segt_id, wire, strip] = \
                            -sxhr_z[sec, md_segt_id, wire, strip] * segX_s + \
                            szhr_z[sec, md_segt_id, wire, strip] * segX_c

                        # rotation of the segment around the Y axis + translation
                        sxhr_zxy[sec, md_segt_id, wire, strip] = \
                            sxhr_zx[sec, md_segt_id, wire, strip] + segX

                        syhr_zxy[sec, md_segt_id, wire, strip] = \
                            syhr_zx[sec, md_segt_id, wire, strip] * segY_c - \
                            szhr_zx[sec, md_segt_id, wire, strip] * segY_s + segY[segment]

                        szhr_zxy[sec, md_segt_id, wire, strip] = \
                            syhr_zx[sec, md_segt_id, wire, strip] * segY_s + \
                            szhr_zx[sec, md_segt_id, wire, strip] * segY_c + segZ

                        # translation module
                        mhrx[sec, md_segt_id, wire, strip] = \
                            sxhr_zxy[sec, md_segt_id, wire, strip] + m_hrx[md]

                        mhry[sec, md_segt_id, wire, strip] = \
                            syhr_zxy[sec, md_segt_id, wire, strip] + m_hry[md]

                        mhrz[sec, md_segt_id, wire, strip] = \
                            szhr_zxy[sec, md_segt_id, wire, strip] + m_hrz

                        # rotation of the module around the Z-axis
                        secx[sec, md_segt_id, wire, strip] = \
                            -mhry[sec, md_segt_id, wire, strip] * secZ_s + \
                            mhrx[sec, md_segt_id, wire, strip] * secZ_c

                        secy[sec, md_segt_id, wire, strip] = \
                            mhry[sec, md_segt_id, wire, strip] * secZ_c + \
                            mhrx[sec, md_segt_id, wire, strip] * secZ_s

                        secz[sec, md_segt_id, wire, strip] = \
                            mhrz[sec, md_segt_id, wire, strip]

                        # Legend:
                        # 9 = 'SANS'
                        # sector no, module no, segment no, wire no, strip no, counter no
                        temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                            9, sec + 1, md + 1, segment + 1, wire + 1, strip + 1, 1
                        )

                        # Legend: x,y,z voxel centers
                        # the positions are first calculated for the High Resolution
                        # A rotation of 180 deg around Y axis is applied below
                        temp1 = '%.2f\t%.2f\t%.2f' % (
                            -secx[sec, md_segt_id, wire, strip],
                            secy[sec, md_segt_id, wire, strip],
                            -secz[sec, md_segt_id, wire, strip]
                        )

                        # Legend:
                        # voxel dimensions to be used to generate Nexus
                        temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                            0,
                            2 * GLxSANS,
                            2 * GLxSANS,
                            2 * GLySANS,
                            2 * GLySANS,
                            2 * GLzSANS)

                        temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                            angZ, -angX, -angY, shr_angle[sec], 0, 0, 0, 0, 0
                        )
                        stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                        f.writelines(stringa)

    # voxels in the 'bottom counter'
    for sec in range(no_sectorsSANS):

        secZ_s = np.sin(np.deg2rad(shr_angle[sec]))
        secZ_c = np.cos(np.deg2rad(shr_angle[sec]))

        for md in range(no_modulesSANS):

            for segment in range(nSANS):

                angZ = 90
                angX = tilt_thetaSANS
                angY = tilt_phiSANS

                segZ_s = np.sin(np.deg2rad(angZ))
                segZ_c = np.cos(np.deg2rad(angZ))
                segY_s = np.sin(np.deg2rad(angY))
                segY_c = np.cos(np.deg2rad(angY))
                segX_s = np.sin(np.deg2rad(angX))
                segX_c = np.cos(np.deg2rad(angX))

                md_segt_id = 100 * md + segment

                for strip in range(n_stripsSANS):
                    for wire in range(n_wiresSANS):
                        # rotation of each segment of the module by angZ
                        # around the Z-axis
                        sxhr_z[sec, md_segt_id, wire, strip] = \
                            -voxelYYhr[wire, strip] * segZ_s - \
                            (voxelXXhr[wire, strip] +
                             voxelXXchr[wire, strip]) * segZ_c

                        syhr_z[sec, md_segt_id, wire, strip] = \
                            voxelYYhr[wire, strip] * segZ_c - \
                            (voxelXXhr[wire, strip] +
                             voxelXXchr[wire, strip]) * segZ_s

                        szhr_z[sec, md_segt_id, wire, strip] = voxelZZhr[wire, strip]

                        # rotation of the segment around the X axis
                        sxhr_zx[sec, md_segt_id, wire, strip] = \
                            sxhr_z[sec, md_segt_id, wire, strip] * segX_c + \
                            szhr_z[sec, md_segt_id, wire, strip] * segX_s

                        syhr_zx[sec, md_segt_id, wire, strip] = \
                            syhr_z[sec, md_segt_id, wire, strip]

                        szhr_zx[sec, md_segt_id, wire, strip] = \
                            -sxhr_z[sec, md_segt_id, wire, strip] * segX_s + \
                            szhr_z[sec, md_segt_id, wire, strip] * segX_c

                        # rotation of the segment around the Y axis + translation
                        sxhr_zxy[sec, md_segt_id, wire, strip] = \
                            sxhr_zx[sec, md_segt_id, wire, strip] + segX

                        syhr_zxy[sec, md_segt_id, wire, strip] = \
                            syhr_zx[sec, md_segt_id, wire, strip] * segY_c - \
                            szhr_zx[sec, md_segt_id, wire, strip] * segY_s + segY[segment]

                        szhr_zxy[sec, md_segt_id, wire, strip] = \
                            syhr_zx[sec, md_segt_id, wire, strip] * segY_s + \
                            szhr_zx[sec, md_segt_id, wire, strip] * segY_c + segZ

                        # translation module
                        mhrx[sec, md_segt_id, wire, strip] = \
                            sxhr_zxy[sec, md_segt_id, wire, strip] + m_hrx[md]

                        mhry[sec, md_segt_id, wire, strip] = \
                            syhr_zxy[sec, md_segt_id, wire, strip] + m_hry[md]

                        mhrz[sec, md_segt_id, wire, strip] = \
                            szhr_zxy[sec, md_segt_id, wire, strip] + m_hrz

                        # rotation of the module around the Z-axis
                        secx[sec, md_segt_id, wire, strip] = \
                            -mhry[sec, md_segt_id, wire, strip] * secZ_s + \
                            mhrx[sec, md_segt_id, wire, strip] * secZ_c

                        secy[sec, md_segt_id, wire, strip] = \
                            mhry[sec, md_segt_id, wire, strip] * secZ_c + \
                            mhrx[sec, md_segt_id, wire, strip] * secZ_s

                        secz[sec, md_segt_id, wire, strip] = mhrz[sec, md_segt_id, wire, strip]

                        # Legend:
                        # 9 = 'SANS'
                        # sector no, module no, segment no, wire no, strip no, counter no
                        temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                            9, sec + 1, md + 1, segment + 1,
                            wire + 1, strip + 1, 2
                        )

                        # Legend: x,y,z voxel centers
                        # the positions are first calculated for the High Resolution
                        # A rotation of 180 deg around Y axis is applied below
                        temp1 = '%.2f\t%.2f\t%.2f' % (
                            -secx[sec, md_segt_id, wire, strip],
                            secy[sec, md_segt_id, wire, strip],
                            -secz[sec, md_segt_id, wire, strip]
                        )

                        # Legend:
                        # voxel dimensions to be used to generate Nexus
                        temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                            0,
                            2 * GLxSANS,
                            2 * GLxSANS,
                            2 * GLySANS,
                            2 * GLySANS,
                            2 * GLzSANS)

                        temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                            angZ, -angX, -angY, shr_angle[sec], 0, 0, 0, 0, 0
                        )

                        stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                        f.writelines(stringa)
