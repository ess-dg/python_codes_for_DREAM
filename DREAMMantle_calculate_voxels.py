#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Mar 14 12:19:10 2022

@author: irinastefanescu
"""

import numpy as np
import sys
import globals

globals.initialize()

np.set_printoptions(precision=4)

"""
 Uses the engineering specifications from the company CDT and
 generates the centers and dimensions of the voxels for the
 DREAM Mantle segment. A Mantle detector segment consists
 of 2 wire grids mounted on each side of a common segmented
 cathode. The whole assembly is enclosed in an Al
 housing with trapezoidal shape. The formalism used here to
 calculate the detector voxels is similar to the implementation
 of the Mantle segment in GEANT4.

First version of this script written by Irina Stefanescu, ESS DG.
March 2022.

Segment engineering specifications received from CDT

data is in mm

"""

"""*******User parameters ************************"""

# no_modulesM = 10  # no of modules in the frame

"""**************************"""

Aly1 = 2335 / 2      # segment length front, sample side
Aly2 = 3400 / 2      # segment length back, readout side
Alz = 368 / 2        # segment depth
Alx1 = 22.4 / 2      # segment width front, sample side
Alx2 = 30.2 / 2      # segment width back, readout side

Althick = 0.3       # thickness Al cathode material and housing
Bthick = 0.0011     # thickness Boron coating

tilt_theta = -10  # tilt_angle=0 when generating the IDF file!

margin = 10
shieldz = 25  # length of the shielding block at the back of the segment
dist_det = 1091  # distance detector front from the sample

dthetaM = np.deg2rad(0.35)  # cathode strip resolution

n_wiresM = 32   # no of wires
n_stripsM = 256  # no of strips
nM = 6  # no of segments per module

dist_between_segs = 2      # distance between neighboring segments in mm

# 0.5*thickness cathode substrate made of Al
CathSubstrX1 = Althick / 2
# X1 & X2 refer to the 2 sides of the cathode. Both sides are coated with Boron
CathSubstrX2 = Althick / 2
# 0.%*thickness Boron layer coated on the cathode
CathConvX1 = Bthick / 2
CathConvX2 = Bthick / 2

init_angleM = 30  # start angle in [deg] for placing the modules in the frame

"""
 the next lines are for calculating the dimensions of
 the inner volume that will be filled with gas voxels
 the principle is similar to the Russian Matryoschka
 dolls (nested dolls)
 From the Al trapezoid representing the segment
 housing subtract the trapezoid
 representing the Boron coating on the side walls
 and the remaining volume is the gas volume,
 therefore the "G" letter used to name
 the variables relevant to the gas voxels
"""

theta_y = np.arctan((Aly2 - Aly1) / 2 / Alz)  # the wings of the segment
theta_x = np.arctan((Alx2 - Alx1) / 2 / Alz)

# Boron trapezoid
By1 = Aly1 - Althick / np.cos(theta_y)
Bz = Alz - Althick
Bx1 = Alx1 - Althick / np.cos(theta_x)

By2 = Aly2 - Althick / np.cos(theta_y)
Bx2 = Alx2 - Althick / np.cos(theta_x)

# Gas volume available to the gas voxels

Gz = Bz - Bthick

Gy1 = By1 - Bthick / np.cos(theta_y)
Gx1 = Bx1 - Bthick / np.cos(theta_x)

Gy2 = By2 - Bthick / np.cos(theta_y)
Gx2 = Bx2 - Bthick / np.cos(theta_x)

CathSubstrY1 = By1 - margin
CathSubstrY2 = By2 - margin
CathSubstrZ = Bz - margin

CathConvY1 = CathSubstrY1
CathConvY2 = CathSubstrY2
CathConvZ = CathSubstrZ

# calculate the dimensions of the gas voxels

xx1 = Gx1 - 2 * CathConvX1 - 3 * CathSubstrX1
xx2 = Gx2 - 2 * CathConvX2 - 3 * CathSubstrX2

# wire pitch
ixx = (xx2 - xx1) / 2 / n_wiresM

# determined by wire pitch
# subtract a margin of 0.1 mm to avoid overlaps
# with housing/coating
izz = 2 * Gz / n_wiresM - 0.1

shp = (n_stripsM//2, n_wiresM)

GLz = np.zeros(shp)
GLy1 = np.zeros(shp)
GLy2 = np.zeros(shp)
GLx1 = np.zeros(shp)
GLx2 = np.zeros(shp)


for strip in range(n_stripsM//2):  # loop over 1/2*strips

    for wire in range(n_wiresM):  # loop over all 32 wires
        # voxel depth
        GLz[strip, wire] = izz / 2

        # voxel height front
        GLy1[strip, wire] = xx1 + wire * ixx

        # voxel height back
        GLy2[strip, wire] = xx1 + (wire + 1) * ixx

        # voxel width front
        GLx1[strip, wire] = (dist_det + wire * izz) * \
                            (np.tan((strip + 1) * dthetaM)
                             - np.tan(strip * dthetaM))

        # voxel width back
        GLx2[strip, wire] = (dist_det + (wire + 1) * izz) * \
                            (np.tan((strip + 1) * dthetaM)
                             - np.tan(strip * dthetaM))

# push up or down the gas voxels to avoid overlap with the cathode components
voxelYd = 4 * CathConvX1 + 2 * CathSubstrX1
voxelYu = 2 * CathConvX1 + 2 * CathSubstrX1

# calculate the centers of the voxels in initial position
voxelZZM = np.zeros(shp)
voxelXXM = np.zeros(n_wiresM)

voxelZZM[0, :] = GLx2[0, :] * 0.5

# not used, but keep it for historical reasons
voxelYYM = GLy1[0, :] + GLy2[0, :]

for j in range(n_wiresM):
    # segment depth
    voxelXXM[j] = Althick + Bthick + j * izz + izz/2

for wire in range(n_wiresM):
    for strip in range(n_stripsM//2-1):
        voxelZZM[strip+1, wire] = voxelZZM[strip, wire] + GLx2[strip+1, wire]

# calculate the solid angle coverage by the segment

distY_seg = 2 * GLy1[0][n_wiresM-1] + 4 * (Althick + Bthick) + voxelYd

seg_middle = dist_det + Alz + shieldz

# Alz and shieldz are already 'halved' when defined;
# added 2 mm between segments in order to avoid the
# volume overlap claimed by GEANT4

st_angle = np.arcsin((distY_seg + dist_between_segs) / seg_middle)    # solid angle

# calculate the segment positions in the detector frame
segX = np.zeros(nM)
segY = np.zeros(nM)

for seg_no in range(nM):
    segX[seg_no] = dist_det * np.cos(seg_no * st_angle)
    segY[seg_no] = dist_det * np.sin(seg_no * st_angle)


no_modulesM = globals.no_modulesM

# calculate the lookup table
shp = (100 * no_modulesM + nM, n_wiresM, n_stripsM//2)

# voxel positions after placing the segment in the frame
sxM = np.zeros(shp)
syM = np.zeros(shp)
szM = np.zeros(shp)
# voxel positions after placing the module in the frame
mxM = np.zeros(shp)
myM = np.zeros(shp)
mzM = np.zeros(shp)

file_path = sys.argv[1]

f = open(file_path, "a")

"""
 The following part of the code is a repetition of a block consisting
  of 4 'for' loops.
 The reason for this is that the segment contains 2 identical wire
 counters on each side of the segmented cathode.
 I call these counters 'top' and 'bottom'. Furthermore, in each
 counter there are 256 strips, which are mirrored  with respect to
 the line corresponding to the center of the segment.
 This means that 128 strips on the left side of the segment center
 are identical in size with the 128 strips on the right and differ
 by the sign of the twist angle dthetaM.
 
 data is saved in the file starting with the voxels from the 
 center of the segment, i.e., strip#128

"""

# strips top right

for md in range(no_modulesM):

    angM = -md * nM * np.rad2deg(st_angle) + init_angleM
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nM):

        ang = segment*np.rad2deg(st_angle) - tilt_theta
        segZ_s = np.sin(np.deg2rad(ang))
        segZ_c = np.cos(np.deg2rad(ang))

        md_segt_id = 100 * md + segment

        for wire in range(n_wiresM):
            for strip in range(n_stripsM//2):

                # rotation of each segment of the module by angM
                # around the Z-axis followed
                # by a translation by sin X and Y
                sxM[md_segt_id, wire, strip] = \
                    segX[segment] + voxelXXM[wire] * segZ_c \
                    - (-voxelYYM[wire] * 0.25 + voxelYu) * segZ_s

                syM[md_segt_id, wire, strip] = \
                    segY[segment] + voxelXXM[wire] * segZ_s \
                    + (-voxelYYM[wire] * 0.25 + voxelYu) * segZ_c

                szM[md_segt_id, wire, strip] = - voxelZZM[strip, wire]

                # rotation of the whole module by angM around the Z-axis
                mxM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_c \
                    - syM[md_segt_id, wire, strip] * mZ_s

                myM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_s \
                    + syM[md_segt_id, wire, strip] * mZ_c

                mzM[md_segt_id, wire, strip] = szM[md_segt_id, wire, strip]

                # Legend:
                # 7 = 'Mantle',
                # 1 = sectors number, always 1 for Mantle
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    7, 1, md + 1, segment + 1, wire + 1, n_stripsM//2 - strip, 1
                )

                # Legend: x,y,z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    mxM[md_segt_id, wire, strip],
                    myM[md_segt_id, wire, strip],
                    mzM[md_segt_id, wire, strip]
                )

                # Legend:
                # voxel dimensions to be used to generate Nexus and
                # IDF files for Mantle detector
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    strip * dthetaM,
                    GLx1[strip, wire],
                    GLx2[strip, wire],
                    GLy1[strip, wire],
                    GLy2[strip, wire],
                    2*GLz[strip, wire])
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    90, angM, ang, 0, 0, 0, 0, 0, 0
                    )
                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                f.writelines(stringa)

# strips top left

for md in range(no_modulesM):

    angM = -md * nM * np.rad2deg(st_angle) + init_angleM
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nM):

        ang = segment * np.rad2deg(st_angle) - tilt_theta
        segZ_s = np.sin(np.deg2rad(ang))
        segZ_c = np.cos(np.deg2rad(ang))

        md_segt_id = 100 * md + segment

        for wire in range(n_wiresM):
            for strip in range(n_stripsM//2):

                # rotation of each segment of the module by angM
                # around the Z-axis followed
                # by a translation by sin X and Y
                sxM[md_segt_id, wire, strip] = \
                    segX[segment] + voxelXXM[wire] * segZ_c \
                    - (-voxelYYM[wire] * 0.25 + voxelYu) * segZ_s

                syM[md_segt_id, wire, strip] = \
                    segY[segment] + voxelXXM[wire] * segZ_s \
                    + (-voxelYYM[wire] * 0.25 + voxelYu) * segZ_c

                szM[md_segt_id, wire, strip] = \
                    voxelZZM[strip, wire]

                # rotation of the whole module by angM around the Z-axis
                mxM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_c \
                    - syM[md_segt_id, wire, strip] * mZ_s

                myM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_s \
                    + syM[md_segt_id, wire, strip] * mZ_c

                mzM[md_segt_id, wire, strip] = \
                    szM[md_segt_id, wire, strip]

                # Legend:
                # 7 = 'Mantle',
                # 1 = sectors number, always 1 for Mantle
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    7, 1, md + 1, segment + 1, 
                    wire + 1, n_stripsM//2 + strip + 1, 1
                )

                # Legend: x,y,z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    mxM[md_segt_id, wire, strip],
                    myM[md_segt_id, wire, strip],
                    mzM[md_segt_id, wire, strip]
                    )

                # Legend:
                # voxel dimensions to be used to generate Nexus
                # and IDF files for Mantle detector
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    - strip*dthetaM,
                    GLx1[strip, wire],
                    GLx2[strip, wire],
                    GLy1[strip, wire],
                    GLy2[strip, wire],
                    2 * GLz[strip, wire]
                    )
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    90, angM, ang, 0, 0, 0, 0, 0, 0
                    )
                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                f.writelines(stringa)


# strips bottom right

for md in range(no_modulesM):

    angM = -md * nM * np.rad2deg(st_angle) + init_angleM
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))
 
    for segment in range(nM):

        ang = segment * np.rad2deg(st_angle) - tilt_theta
        segZ_s = np.sin(np.deg2rad(ang))
        segZ_c = np.cos(np.deg2rad(ang))

        md_segt_id = 100 * md + segment
 
        for wire in range(n_wiresM):
            for strip in range(n_stripsM//2):

                # rotation of each segment of the module by angM
                # around the Z-axis followed
                # by a translation by sin X and Y
                sxM[md_segt_id, wire, strip] = \
                    segX[segment] + voxelXXM[wire] * segZ_c \
                    - (voxelYYM[wire] * 0.25 - voxelYd) * segZ_s
                syM[md_segt_id, wire, strip] = \
                    segY[segment] + voxelXXM[wire] * segZ_s \
                    + (voxelYYM[wire] * 0.25 - voxelYd) * segZ_c
                szM[md_segt_id, wire, strip] = - voxelZZM[strip, wire]

                # rotation of the whole module by angM around the Z-axis
                mxM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_c \
                    - syM[md_segt_id, wire, strip] * mZ_s
                myM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_s \
                    + syM[md_segt_id, wire, strip] * mZ_c
                mzM[md_segt_id, wire, strip] = szM[md_segt_id, wire, strip]

                # Legend:
                # 7 = 'Mantle',
                # 1 = sectors number, always 1 for Mantle
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    7, 1, md + 1, segment + 1, wire + 1, n_stripsM//2 - strip, 2
                )

                # Legend: x, y, z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    mxM[md_segt_id, wire, strip],
                    myM[md_segt_id, wire, strip],
                    mzM[md_segt_id, wire, strip])

                # Legend:
                # voxel dimensions to be used to generate Nexus
                # and IDF files for Mantle detector
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    strip * dthetaM,
                    GLx1[strip, wire],
                    GLx2[strip, wire],
                    GLy1[strip, wire],
                    GLy2[strip, wire],
                    2 * GLz[strip, wire]
                    )

                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    90, angM, ang, 0, 0, 0, 0, 0, 0
                    )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                f.writelines(stringa)

# strips bottom left

for md in range(no_modulesM):

    angM = -md * nM * np.rad2deg(st_angle) + init_angleM
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nM):

        ang = segment * np.rad2deg(st_angle) - tilt_theta
        segZ_s = np.sin(np.deg2rad(ang))
        segZ_c = np.cos(np.deg2rad(ang))

        md_segt_id = 100 * md + segment

        for wire in range(n_wiresM):
            for strip in range(n_stripsM//2):
                # rotation of each segment of the module by angM
                # around the Z-axis followed by a
                # translation by sin X and Y
                sxM[md_segt_id, wire, strip] = \
                    segX[segment] + voxelXXM[wire] * segZ_c \
                    - (voxelYYM[wire] * 0.25 - voxelYd) * segZ_s

                syM[md_segt_id, wire, strip] = \
                    segY[segment] + voxelXXM[wire] * segZ_s \
                    + (voxelYYM[wire] * 0.25 - voxelYd) * segZ_c

                szM[md_segt_id, wire, strip] = \
                    voxelZZM[strip, wire]

                # rotation of the whole module by angM around the Z-axis
                mxM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_c \
                    - syM[md_segt_id, wire, strip] * mZ_s

                myM[md_segt_id, wire, strip] = \
                    sxM[md_segt_id, wire, strip] * mZ_s \
                    + syM[md_segt_id, wire, strip] * mZ_c

                mzM[md_segt_id, wire, strip] = \
                    szM[md_segt_id, wire, strip]

                # Legend:
                # 7 = 'Mantle',
                # 1 = sectors number, always 1 for Mantle
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    7, 1, md + 1, segment + 1, wire + 1, n_stripsM//2 + strip + 1, 2
                )

                # Legend: x, y, z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    mxM[md_segt_id, wire, strip],
                    myM[md_segt_id, wire, strip],
                    mzM[md_segt_id, wire, strip]
                    )

                # Legend:
                # voxel dimensions to be used to generate Nexus
                # and IDF files for Mantle detector
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    - strip*dthetaM,
                    GLx1[strip, wire],
                    GLx2[strip, wire],
                    GLy1[strip, wire],
                    GLy2[strip, wire],
                    2 * GLz[strip, wire]
                    )

                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    90, angM, ang, 0, 0, 0, 0, 0, 0
                    )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                f.writelines(stringa)

f.close()
