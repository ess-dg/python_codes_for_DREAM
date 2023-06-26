#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday May 11 10:45:32 2022

@author: irinastefanescu
"""

import numpy as np

np.set_printoptions(precision=4)

"""
*******The HR for EndCap ******

 Uses the engineering specifications from the company CDT and
 generates the centers and dimensions of the voxels for the
 DREAM HR module. A Mantle detector segment consists
 of 2 wire grids mounted on each side of a common segmented
 cathode. The whole assembly is enclosed in an Al
 housing with trapezoidal shape. The formalism used here to
 calculate the detector voxels is similar to the implementation
 of the Mantle segment in GEANT4.

First version of this script written by Irina Stefanescu, ESS DG.
May 2022.

data is in mm

"""

"""**************User parameters *****************"""

nHR = 8             # no of segments per module
no_modulesHR = 9    # no of modules in the frame
no_sectorsHR = 4    # no of sectors 

"""***********************************************"""

hHR = 170       # segment height front, sample side
sensHR = 313    # segment depth
wHR = 20        # segment width front, sample side

margin = 10      # empty space inside the segment front-back and left-right
marginb = 1.5    # empty space inside the segment top-bottom
shieldz = 25     # length of the shielding block at the back of the segment
mod_dist = 2     # distance between neighboring segments in order to avoid overlap

AlthickH = 0.3      # thickness Al cathode material and housing
BthickHR = 0.0013   # thickness Boron coating

# start position for placing the modules in the frame, integer number
# to multiply with 12*deg
init_angleHR = 0     # init angle in deg

tilt_theta = -10     # tilt_angle in deg
tilt_thetaHR = -1    # tilt_angle in deg
tilt_phiHR = -1      # inclination angle module in deg

n_wiresHR = 16     # no of wires
n_stripsHR = 32    # no of strips

dist_detHR = 2500   # distance from the sample in mm

# start calculations 

# Al housing (box)

AlyHR = (hHR + 1) / 2   # 1/2 height Al housing at entrance window HR
AlzHR = (sensHR + margin) / 2    # depth housing HR
AlxHR = wHR / 2   # 1/2 width segment front in HR

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

ByHR = AlyHR - AlthickH
BzHR = AlzHR - AlthickH
BxHR = AlxHR - AlthickH

# Gas volume available to the gas voxels

GyHR = ByHR
GzHR = BzHR - BthickHR
GxHR = BxHR - BthickHR

CathSubstrYHR = GyHR - marginb / 2
CathSubstrZHR = BzHR - margin / 2

CathConvYHR = CathSubstrYHR
CathConvZHR = CathSubstrZHR
CathConvXHR = BthickHR / 2
CathSubstrX1 = AlthickH / 2

# calculate the dimensions of the gas voxels

xxHR = GxHR - (2*CathConvXHR + AlthickH)

izzHR = sensHR / n_stripsHR         # strip pitch HR in mm, all equal

GLzHR = izzHR / 2
GLyHR = hHR / 2 / n_wiresHR
GLxHR = xxHR / 2

# calculate the centers of the voxels 
shp = (n_wiresHR, n_stripsHR)

voxelXXhr = np.zeros(shp)
voxelXXchr = np.zeros(shp)
voxelYYhr = np.zeros(shp)
voxelZZhr = np.zeros(shp)


for strip in range(n_stripsHR):
    for wire in range(n_wiresHR):
        # fill the voxels from the top to bottom and from front to back
        # voxelYYhr[wire, strip] = (n_wiresHR//2 - wire - 0.5) * 2 * GLyHR
        voxelYYhr[wire, strip] = (n_wiresHR - 2 * wire - 1) * GLyHR
        voxelZZhr[wire, strip] = (n_stripsHR//2 - strip - 1) * izzHR + margin/2
        voxelXXhr[wire, strip] = AlthickH
        voxelXXchr[wire, strip] = GLxHR + 2 * CathConvXHR + CathSubstrX1

# calculate the segment positions in the detector frame

segY = np.zeros(nHR)

segZ = - dist_detHR - AlzHR
segX = 0

for seg_no in range(nHR):
    segY[seg_no] = 80 + (2 * AlxHR + 0.4) * (seg_no + 1)
    
# calculate the module positions in the detector frame

pos_mod = hHR + mod_dist  # allow for 2 mm space between the modules

m_hrx = np.array([0, 0, 0, pos_mod, pos_mod, pos_mod, 2*pos_mod, 2*pos_mod, 3*pos_mod])

m_hry = np.array([0, pos_mod, 2*pos_mod, 0, pos_mod, 2*pos_mod, 0, pos_mod, 0])

m_hrz = 0

shr_angle = np.array([0, -90, 180, 90])

""" calculate the lookup table  """

shp = (no_sectorsHR, 100 * no_modulesHR + nHR, n_wiresHR, n_stripsHR)

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

with open('HR_temp.txt', "a") as f:

    # voxels in the 'top counter'
    for sec in range(no_sectorsHR):

        secZ_s = np.sin(np.deg2rad(shr_angle[sec]))
        secZ_c = np.cos(np.deg2rad(shr_angle[sec]))

        for md in range(no_modulesHR):

            for segment in range(nHR):

                angZ = 90
                angX = tilt_thetaHR
                angY = tilt_phiHR

                segZ_s = np.sin(np.deg2rad(angZ))
                segZ_c = np.cos(np.deg2rad(angZ))
                segY_s = np.sin(np.deg2rad(angY))
                segY_c = np.cos(np.deg2rad(angY))
                segX_s = np.sin(np.deg2rad(angX))
                segX_c = np.cos(np.deg2rad(angX))

                md_segt_id = 100 * md + segment

                for strip in range(n_stripsHR):
                    for wire in range(n_wiresHR):
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
                        # 8 = 'HR'
                        # sector no, module no, segment no, wire no, strip no, counter no
                        temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                            8, sec + 1, md + 1, segment + 1, wire + 1, strip + 1, 1
                        )

                        # Legend: x,y,z voxel centers
                        temp1 = '%.2f\t%.2f\t%.2f' % (
                            secx[sec, md_segt_id, wire, strip],
                            secy[sec, md_segt_id, wire, strip],
                            secz[sec, md_segt_id, wire, strip]
                        )

                        # Legend:
                        # voxel dimensions to be used to generate Nexus
                        temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                            0,
                            2 * GLxHR,
                            2 * GLxHR,
                            2 * GLyHR,
                            2 * GLyHR,
                            2 * GLzHR)

                        temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                            angZ, -angX, -angY, shr_angle[sec], 0, 0, 0, 0, 0
                        )
                        stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                        f.writelines(stringa)

    # voxels in the 'bottom counter'
    for sec in range(no_sectorsHR):

        secZ_s = np.sin(np.deg2rad(shr_angle[sec]))
        secZ_c = np.cos(np.deg2rad(shr_angle[sec]))

        for md in range(no_modulesHR):

            for segment in range(nHR):

                angZ = 90
                angX = tilt_thetaHR
                angY = tilt_phiHR

                segZ_s = np.sin(np.deg2rad(angZ))
                segZ_c = np.cos(np.deg2rad(angZ))
                segY_s = np.sin(np.deg2rad(angY))
                segY_c = np.cos(np.deg2rad(angY))
                segX_s = np.sin(np.deg2rad(angX))
                segX_c = np.cos(np.deg2rad(angX))

                md_segt_id = 100 * md + segment

                for strip in range(n_stripsHR):
                    for wire in range(n_wiresHR):
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
                        # 8 = 'HR'
                        # sector no, module no, segment no, wire no, strip no, counter no
                        temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                            8, sec + 1, md + 1, segment + 1,
                            wire + 1, strip + 1, 2
                        )

                        # Legend: x,y,z voxel centers
                        temp1 = '%.2f\t%.2f\t%.2f' % (
                            secx[sec, md_segt_id, wire, strip],
                            secy[sec, md_segt_id, wire, strip],
                            secz[sec, md_segt_id, wire, strip]
                        )

                        # Legend:
                        # voxel dimensions to be used to generate Nexus
                        temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                            0,
                            2 * GLxHR,
                            2 * GLxHR,
                            2 * GLyHR,
                            2 * GLyHR,
                            2 * GLzHR)

                        temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                            angZ, -angX, -angY, shr_angle[sec], 0, 0, 0, 0, 0
                        )

                        stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                        f.writelines(stringa)
