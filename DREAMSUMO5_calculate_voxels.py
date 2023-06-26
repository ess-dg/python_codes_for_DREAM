#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 10:13:30 2022

@author: irinastefanescu
"""

import numpy as np
import globals

globals.initialize()

np.set_printoptions(precision=4)

"""
*******The SUMO4 for EndCap ******

 Uses the engineering specifications from the company CDT and
 generates the centers and dimensions of the voxels for the
 DREAM SUMO5 module. A Mantle detector segment consists
 of 2 wire grids mounted on each side of a common segmented
 cathode. The whole assembly is enclosed in an Al
 housing with trapezoidal shape. The formalism used here to
 calculate the detector voxels is similar to the implementation
 of the Mantle segment in GEANT4.

First version of this script written by Irina Stefanescu, ESS DG.
March 2022.

Segment engineering specifications 
from Table 2, page 5, document "Description of work DREAM Milestone 2.1"

data is in mm

"""

hfS5 = 160        # segment height front, sample side
hbS5 = 201       # segment height back, readout side
sensS5 = 330      # segment depth
wbfS5 = 130       # segment width bottom front, sample side
wtfS5 = 161        # segment width top front, sample side
wbbS5 = 145      # segment width bottom back, readout side
wtbS5 = 195     # segment width top back, readout side

margin = 10      # empty space inside the segment front-back and left-right
marginb = 1.5    # empty space inside the segment top-bottom
shieldz = 25     # length of the shielding block at the back of the segment

Althick = 0.3       # thickness Al cathode material and housing
Bthick = 0.0011     # thickness Boron coating

# start position for placing the modules in the frame, integer number
# to multiply with 12*deg
index_rot = 10

tilt_theta = -10     # tilt_angle in deg
tilt_phiS5 = 23      # inclination angle module in deg
dphi = np.deg2rad(12)  # angular coverage in phi for the whole module
s_angle_S5 = 0.39    # rotation segment around the Y axis
s_offset_S5 = 0.18   # additional angle for the Z rotation of the segment
alpha_physS5 = 7.21  # I don't remember what this is :-(

n_wires = 16     # no of wires
n_strips = 16    # no of strips
nS5 = 8           # no of segments per module

offsetX_S5 = 40.       # x-distance SUMO5 from the sample
offsetY_S5 = 790.      # y-distance SUMO5 from the sample
offsetZ_S5 = -1320.    # z-distance SUMO5 from the sample

# start calculations 

# this function recalculates the shape parameters in order
# to avoid the crash of the G4Trap class
# I must have it here too, for debugging purposes

Aly1S5 = hfS5 / 2  # 1/2 height Al housing at entrance window SUMO5
Aly2S5 = hbS5 / 2  # 1/2 height Al housing at the back SUMO5
AlzS5 = (sensS5 + margin) / 2  # depth housing SUMO5
  
# subtract (no_of_segments - 1)*mm from the widths to account for the space between the segments  

# 1/2 width Al housing at the bottom front in SUMO5
Alx1bS5 = (wbfS5 - 15) / 2 / nS5  
# 1/2 width segment at the top front in SUMO5
Alx1tS5 = (wtfS5 - 15) / 2 / nS5 
#  1/2 width segment at the bottom back SUMO5
Alx2bS5 = (wbbS5 - 25) / 2 / nS5 
#  1/2 width segment at the top back SUMO5
Alx2tS5 = (wtbS5 - 12) / 2 / nS5 

# recalculated shape parameters, same names
Aly1S5, Alx1bS5, Alx1tS5, Aly2S5, Alx2bS5, Alx2tS5 = \
    globals.match2geantvars([Aly1S5, Alx1bS5, Alx1tS5, Aly2S5, Alx2bS5, Alx2tS5])

# 0.5*thickness cathode substrate made of Al
CathSubstrX1 = Althick / 2
# X1 & X2 refer to the 2 sides of the cathode. Both sides are coated with Boron
CathSubstrX2 = Althick / 2
# 0.5*thickness Boron layer coated on the cathode
CathConvX1 = Bthick / 2
CathConvX2 = Bthick / 2

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

# Boron trapezoid

By1S5 = Aly1S5 - Althick / np.cos(dphi/2)
By2S5 = Aly2S5 - Althick / np.cos(dphi/2)
BzS5 = AlzS5 - Althick
Bx1tS5 = Alx1tS5 - Althick / np.cos(dphi/2)
Bx2tS5 = Alx2tS5 - Althick / np.cos(dphi/2)
Bx1bS5 = Alx1bS5 - Althick / np.cos(dphi/2)
Bx2bS5 = Alx2bS5 - Althick / np.cos(dphi/2)

# recalculated shape parameters, same names

By1S5, Bx1bS5, Bx1tS5, By2S5, Bx2bS5, Bx2tS5 = \
    globals.match2geantvars([By1S5, Bx1bS5, Bx1tS5, By2S5, Bx2bS5, Bx2tS5])

# Gas volume available to the gas voxels

Gy1S5 = By1S5
Gy2S5 = By2S5
GzS5 = BzS5 - Bthick
Gx1tS5 = Bx1tS5 - Bthick / np.cos(dphi/2)
Gx2tS5 = Bx2tS5 - Bthick / np.cos(dphi/2)
Gx1bS5 = Bx1bS5 - Bthick / np.cos(dphi/2)
Gx2bS5 = Bx2bS5 - Bthick / np.cos(dphi/2)

Gy1S5, Gx1bS5, Gx1tS5, Gy2S5, Gx2bS5, Gx2tS5 = \
    globals.match2geantvars([Gy1S5, Gx1bS5, Gx1tS5, Gy2S5, Gx2bS5, Gx2tS5])

CathSubstrY1S5 = Gy1S5 - marginb / 2
CathSubstrY2S5 = Gy2S5 - marginb / 2
CathSubstrZS5 = BzS5 - margin / 2

CathConvY1S5 = CathSubstrY1S5
CathConvY2S5 = CathSubstrY2S5
CathConvZS5 = CathSubstrZS5

# calculate the dimensions of the gas voxels

xx1tS5 = Gx1tS5 - 2
xx2tS5 = Gx2tS5 - 2
xx1bS5 = Gx1bS5 - 2
xx2bS5 = Gx2bS5 - 2

eta_b = (2 * (xx1tS5 - xx1bS5) / Gy1S5)
eta_t = ((xx2bS5 - xx1bS5) / sensS5)
eta_w = ((Gy2S5 - Gy1S5) / sensS5)  # Gy1 and Gy2 already halved

izzS5 = sensS5 / n_strips  # strip pitch SUMO5 in mm, all equal

dthetaS5 = alpha_physS5 / n_wires  # wire pitch SUMO5 in deg, all equal

shp = (n_strips, n_wires)

GLzS5 = np.zeros(shp)
GLy1S5 = np.zeros(shp)
GLy2S5 = np.zeros(shp)
GLx1bS5 = np.zeros(shp)
GLx2bS5 = np.zeros(shp)
GLx1bbS5 = np.zeros(shp)
GLx2bbS5 = np.zeros(shp)

# loop goes from 0 to max-1!!
# wires run parallel to the beam axis and in fan out geometry
# wire pitch smaller at the segment front and larger at the back
# also, wire pitch is symmetric with respect to the segment center
# the calculations starts from the center of the segment, 
# so wire=1 is close to the center, wire = 8 is at the bottom (or top)

for strip in range(n_strips):  # loop over strips

    for wire in range(n_wires//2):  # loop over all 1/2*wires    
        # voxel depth
        GLzS5[wire, strip] = izzS5 / 2
        GLzS5[n_wires//2 + wire, strip] = izzS5 / 2
        # wire pitch front of the voxel
        # GLy1S5[wire, strip] = (strip * 2 * izzS5 * eta_w + 2 * (Gy1S5 - 4)) / 2 / n_wires
        GLy1S5[wire, strip] = (strip * izzS5 * eta_w + Gy1S5 - 4) / n_wires
        # wire pitch, back of the voxel
        # GLy2S5[wire, strip] = ((strip + 1) * 2 * izzS5 * eta_w +
        #                        2 * (Gy1S5 - 4)) / 2 / n_wires
        GLy2S5[wire, strip] = ((strip + 1) * izzS5 * eta_w + Gy1S5 - 4) / n_wires
        GLy1S5[n_wires//2 + wire, strip] = GLy1S5[wire, strip]
        GLy2S5[n_wires//2 + wire, strip] = GLy2S5[wire, strip] 

        # for the bottom 8 wires
        GLx1bbS5[wire, strip] = 0.5 * (2 * xx1bS5 +
                                       wire * GLy1S5[wire, strip] * eta_b)
        GLx2bbS5[wire, strip] = GLx1bbS5[wire, strip] + 0.5 * strip * izzS5 * eta_t

        # for the upper 8 wires

        GLx1bS5[wire, strip] = 0.5 * (2 * xx1bS5 +
                                      (n_wires//2 + wire) *
                                      GLy1S5[wire, strip] * eta_b)
        GLx2bS5[wire, strip] = GLx1bS5[wire, strip] + 0.5 * strip * izzS5 * eta_t
        GLx1bbS5[n_wires//2 + wire, strip] = GLx1bS5[wire, strip]
        GLx2bbS5[n_wires//2 + wire, strip] = GLx2bS5[wire, strip]

# calculate the centers of the voxels 
shp = (n_wires, n_strips)

voxelXX = np.zeros(shp)
voxelXXc = np.zeros(shp)
voxelYY = np.zeros(shp)
voxelZZ = np.zeros(shp)

for strip in range(n_strips):
    for wire in range(n_wires//2):
        # voxels created by the lowest 8 wires
        # fill the voxels from the bottom to the segment center
        # voxelYY[wire, strip] = -(n_wires//2 - wire - 0.5) * 2 * GLy2S5[n_wires//2 - wire, strip]
        voxelYY[wire, strip] = -(n_wires - 2 * wire - 1) * GLy2S5[n_wires // 2 - wire, strip]
        voxelZZ[wire, strip] = \
            (n_strips - strip - 1) * 2 * GLzS5[n_wires//2 - wire, 0] + GLzS5[n_wires//2 - wire, 0] + margin/2
        voxelXX[wire, strip] = Althick
        voxelXXc[wire, strip] = 0.5 * GLx1bbS5[wire, strip]
        # voxels created by the lowest 8 wires
        # fill the voxels from the bottom to the segment center
        voxelYY[n_wires//2 + wire, strip] = (wire + 0.5) * 2 * GLy2S5[wire, strip]
        voxelZZ[n_wires//2 + wire, strip] = \
            (n_strips - strip - 1) * 2 * GLzS5[wire, 0] + GLzS5[wire, 0] + margin/2
        voxelXX[n_wires//2 + wire, strip] = Althick
        voxelXXc[n_wires//2 + wire, strip] = 0.5 * GLx1bS5[wire, strip]

# calculate the segment positions in the detector frame
segX = np.zeros(nS5)
segZ = np.zeros(nS5)
segY = np.zeros(nS5)

for seg_no in range(nS5):
    segX[seg_no] = Alx2tS5 * (nS5 + 1 - 2 * (seg_no + 1))
    segZ[seg_no] = -sensS5 / 2
    
radS5 = np.sqrt(np.power(offsetX_S5, 2) + np.power(offsetY_S5, 2))

no_modules = globals.no_modules

modX = np.zeros(no_modules)
modY = np.zeros(no_modules)
modZ = np.zeros(no_modules)
modZ[:] = offsetZ_S5

for mod in range(no_modules):
    # modX[mod] = radS5 * np.sin(np.deg2rad(-(index_rot + mod) * np.rad2deg(dphi)))
    # modY[mod] = radS5 * np.cos(np.deg2rad(-(index_rot + mod) * np.rad2deg(dphi)))
    modX[mod] = radS5 * np.sin(-(index_rot + mod) * dphi)
    modY[mod] = radS5 * np.cos(-(index_rot + mod) * dphi)

""" calculate the lookup table  """

shp = (100 * no_modules + nS5, n_wires, n_strips)

# voxel positions after placing the segment in the frame
sx_z = np.zeros(shp)
sy_z = np.zeros(shp)
sz_z = np.zeros(shp)

sx_y = np.zeros(shp)
sy_y = np.zeros(shp)
sz_y = np.zeros(shp)

# voxel positions after placing the module in the frame
mxxx = np.zeros(shp)
myyy = np.zeros(shp)
mzzz = np.zeros(shp)

mxx = np.zeros(shp)
myy = np.zeros(shp)
mzz = np.zeros(shp)

mx = np.zeros(shp)
my = np.zeros(shp)
mz = np.zeros(shp)

VX = np.zeros(shp)
VY = np.zeros(shp)
VZ = np.zeros(shp)

XF = np.zeros(shp)
YF = np.zeros(shp)
ZF = np.zeros(shp)

VXF = np.zeros(shp)
VYF = np.zeros(shp)
VZF = np.zeros(shp)

fF = open('Forward_temp.txt', "a")
fB = open('Backward_temp.txt', "a")

mY_s = np.sin(np.deg2rad(-tilt_theta))
mY_c = np.cos(np.deg2rad(-tilt_theta))

mX_s = np.sin(np.deg2rad(tilt_phiS5))
mX_c = np.cos(np.deg2rad(tilt_phiS5))

fY_s = np.sin(np.deg2rad(180))
fY_c = np.cos(np.deg2rad(180))

fZ_s = np.sin(np.deg2rad(90))
fZ_c = np.cos(np.deg2rad(90))

# voxels in the left counter

for md in range(no_modules):
    
    angM = (index_rot + md) * np.rad2deg(dphi)
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nS5):

        angZ = (nS5 + 1 - 2 * (segment + 1)) * (s_angle_S5 + s_offset_S5)
        angY = (nS5 + 1 - 2 * (segment + 1)) * s_angle_S5
        segZ_s = np.sin(np.deg2rad(angZ))
        segZ_c = np.cos(np.deg2rad(angZ))
        segY_s = np.sin(np.deg2rad(angY))
        segY_c = np.cos(np.deg2rad(angY))

        # identification using module and segment numbers
        md_segt_id = 100 * md + segment

        for strip in range(n_strips):
            for wire in range(n_wires):

                # rotation of each segment of the module by angZ
                # around the Z-axis followed
                # by a rotation of each segment by angY around Y-axis
                sx_z[md_segt_id, wire, strip] = \
                      (segX[segment] - voxelXX[wire, strip] -
                       voxelXXc[wire, strip]) + voxelYY[wire, strip] * segZ_s

                sy_z[md_segt_id, wire, strip] = \
                    (segY[segment] + voxelYY[wire, strip]) * segZ_c

                sz_z[md_segt_id, wire, strip] = \
                    segZ[segment] + voxelZZ[wire, strip]

                sx_y[md_segt_id, wire, strip] = \
                    voxelYY[wire, strip] * segZ_s * segY_c - \
                    voxelZZ[wire, strip] * segY_s + segX[segment] - \
                    voxelXX[wire, strip] - voxelXXc[wire, strip]

                sy_y[md_segt_id, wire, strip] = sy_z[md_segt_id, wire, strip]

                sz_y[md_segt_id, wire, strip] = \
                    voxelYY[wire, strip] * segZ_s * segY_s + \
                    voxelZZ[wire, strip] * segY_c + segZ[segment]
                    
                # rotation of the module around the Y-axis
                mxxx[md_segt_id, wire, strip] = \
                    sx_y[md_segt_id, wire, strip] * mY_c + \
                    sz_y[md_segt_id, wire, strip] * mY_s

                myyy[md_segt_id, wire, strip] = sy_y[md_segt_id, wire, strip]

                mzzz[md_segt_id, wire, strip] = \
                    -sx_y[md_segt_id, wire, strip] * mY_s + \
                    sz_y[md_segt_id, wire, strip] * mY_c

                # rotation of the module around the X-axis
                mxx[md_segt_id, wire, strip] = mxxx[md_segt_id, wire, strip]

                myy[md_segt_id, wire, strip] = \
                    myyy[md_segt_id, wire, strip] * mX_c - \
                    mzzz[md_segt_id, wire, strip] * mX_s

                mzz[md_segt_id, wire, strip] = \
                    myyy[md_segt_id, wire, strip] * mX_s + \
                    mzzz[md_segt_id, wire, strip] * mX_c

                # rotation of the module around the Z-axis
                mx[md_segt_id, wire, strip] = \
                    mxx[md_segt_id, wire, strip] * mZ_c - \
                    myy[md_segt_id, wire, strip] * mZ_s

                my[md_segt_id, wire, strip] = \
                    mxx[md_segt_id, wire, strip] * mZ_s + \
                    myy[md_segt_id, wire, strip] * mZ_c

                mz[md_segt_id, wire, strip] = mzz[md_segt_id, wire, strip]
                    
                # translation of the module
                VX[md_segt_id, wire, strip] = \
                    modX[md] + mx[md_segt_id, wire, strip]

                VY[md_segt_id, wire, strip] = \
                    modY[md] + my[md_segt_id, wire, strip]

                VZ[md_segt_id, wire, strip] = \
                    modZ[md] + mz[md_segt_id, wire, strip]
                
                # Forward detector, rotation around Y-axis by 180 deg
                XF[md_segt_id, wire, strip] = \
                    VX[md_segt_id, wire, strip] * fY_c + \
                    VZ[md_segt_id, wire, strip] * fY_s

                YF[md_segt_id, wire, strip] = VY[md_segt_id, wire, strip]

                ZF[md_segt_id, wire, strip] = \
                    VZ[md_segt_id, wire, strip] * fY_c - \
                    VX[md_segt_id, wire, strip] * fY_s

                # Forward detector, rotation around Z-axis by 180 deg
                VXF[md_segt_id, wire, strip] = \
                    XF[md_segt_id, wire, strip] * fZ_c - \
                    YF[md_segt_id, wire, strip] * fZ_s

                VYF[md_segt_id, wire, strip] = \
                    XF[md_segt_id, wire, strip] * fZ_s + \
                    YF[md_segt_id, wire, strip] * fZ_c

                VZF[md_segt_id, wire, strip] = ZF[md_segt_id, wire, strip]
                    
                # Legend:
                # 5 = 'SUMO5 Backward', 15 = 'SUMO5 Forward'
                # 1 = sectors number, always 1 for SUMO5 & Forward
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    5, 1, md + 1, segment + 1, wire + 1, strip + 1, 1
                )

                tempF = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    15, 1, md + 1, segment + 1, wire + 1, strip + 1, 1
                )

                # Legend: x,y,z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    VX[md_segt_id, wire, strip],
                    VY[md_segt_id, wire, strip],
                    VZ[md_segt_id, wire, strip]
                )

                tempF1 = '%.2f\t%.2f\t%.2f' % (
                    VXF[md_segt_id, wire, strip],
                    VYF[md_segt_id, wire, strip],
                    VZF[md_segt_id, wire, strip]
                )

                # Legend:
                # voxel dimensions to be used to generate Nexus 
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    np.deg2rad(wire*dthetaS5)*0,
                    GLx1bbS5[wire, strip],
                    GLx2bbS5[wire, strip],
                    2*GLy1S5[wire, strip],
                    2*GLy2S5[wire, strip],
                    2*GLzS5[wire, strip])

                # Legend:
                # rotation angles to put the voxels into the right positions, Backward EndCap
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS5, angM, 0, 0, 0, 0
                )

                # rotation angles to put the voxels into the right positions, Forward EndCap
                tempF3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS5, angM, 180, 90, 0, 0
                )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                stringaf = tempF + '\t' + tempF1 + '\t' + temp2 + '\t' + tempF3

                fB.writelines(stringa)
                fF.writelines(stringaf)

# voxels in the right counter

for md in range(no_modules):
    
    angM = (index_rot + md) * np.rad2deg(dphi)
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nS5):

        md_segt_id = 100 * md + segment

        angZ = (nS5 + 1 - 2 * (segment + 1)) * (s_angle_S5 + s_offset_S5)
        angY = (nS5 + 1 - 2 * (segment + 1)) * s_angle_S5
        segZ_s = np.sin(np.deg2rad(angZ))
        segZ_c = np.cos(np.deg2rad(angZ))
        segY_s = np.sin(np.deg2rad(angY))
        segY_c = np.cos(np.deg2rad(angY))

        for strip in range(n_strips):
            for wire in range(n_wires):

                # rotation of each segment of the module by angZ
                # around the Z-axis followed
                # by a rotation of each segment by angY around Y-axis
                sx_z[md_segt_id, wire, strip] = \
                    (segX[segment] + voxelXX[wire, strip] +
                     voxelXXc[wire, strip]) + voxelYY[wire, strip] * segZ_s

                sy_z[md_segt_id, wire, strip] = \
                    (segY[segment] + voxelYY[wire, strip]) * segZ_c

                sz_z[md_segt_id, wire, strip] = \
                    segZ[segment] + voxelZZ[wire, strip]

                sx_y[md_segt_id, wire, strip] = \
                    voxelYY[wire, strip] * segZ_s * segY_c - \
                    voxelZZ[wire, strip] * segY_s + segX[segment] + \
                    voxelXX[wire, strip] + voxelXXc[wire, strip]

                sy_y[md_segt_id, wire, strip] = sy_z[md_segt_id, wire, strip]

                sz_y[md_segt_id, wire, strip] = \
                    voxelYY[wire, strip] * segZ_s * segY_s + \
                    voxelZZ[wire, strip] * segY_c + segZ[segment]
                    
                # rotation of the module around the Y-axis
                mxxx[md_segt_id, wire, strip] = \
                    sx_y[md_segt_id, wire, strip] * mY_c + \
                    sz_y[md_segt_id, wire, strip] * mY_s

                myyy[md_segt_id, wire, strip] = sy_y[md_segt_id, wire, strip]

                mzzz[md_segt_id, wire, strip] = \
                    -sx_y[md_segt_id, wire, strip] * mY_s + \
                    sz_y[md_segt_id, wire, strip] * mY_c

                # rotation of the module around the X-axis
                mxx[md_segt_id, wire, strip] = mxxx[md_segt_id, wire, strip]

                myy[md_segt_id, wire, strip] = \
                    myyy[md_segt_id, wire, strip] * mX_c - \
                    mzzz[md_segt_id, wire, strip] * mX_s

                mzz[md_segt_id, wire, strip] = \
                    myyy[md_segt_id, wire, strip] * mX_s + \
                    mzzz[md_segt_id, wire, strip] * mX_c

                # rotation of the module around the Z-axis
                mx[md_segt_id, wire, strip] = \
                    mxx[md_segt_id, wire, strip] * mZ_c - \
                    myy[md_segt_id, wire, strip] * mZ_s

                my[md_segt_id, wire, strip] = \
                    mxx[md_segt_id, wire, strip] * mZ_s + \
                    myy[md_segt_id, wire, strip] * mZ_c

                mz[md_segt_id, wire, strip] = mzz[md_segt_id, wire, strip]
                    
                # translation of the module
                VX[md_segt_id, wire, strip] = \
                    modX[md] + mx[md_segt_id, wire, strip]

                VY[md_segt_id, wire, strip] = \
                    modY[md] + my[md_segt_id, wire, strip]

                VZ[md_segt_id, wire, strip] = \
                    modZ[md] + mz[md_segt_id, wire, strip]
                
                # Forward detector, rotation around Y-axis by 90 deg
                XF[md_segt_id, wire, strip] = \
                    VX[md_segt_id, wire, strip] * fY_c + \
                    VZ[md_segt_id, wire, strip] * fY_s

                YF[md_segt_id, wire, strip] = VY[md_segt_id, wire, strip]

                ZF[md_segt_id, wire, strip] = \
                    VZ[md_segt_id, wire, strip] * fY_c - \
                    VX[md_segt_id, wire, strip] * fY_s

                # Forward detector, rotation around Z-axis by 90 deg
                VXF[md_segt_id, wire, strip] = \
                    XF[md_segt_id, wire, strip] * fZ_c - \
                    YF[md_segt_id, wire, strip] * fZ_s

                VYF[md_segt_id, wire, strip] = \
                    XF[md_segt_id, wire, strip] * fZ_s + \
                    YF[md_segt_id, wire, strip] * fZ_c

                VZF[md_segt_id, wire, strip] = ZF[md_segt_id, wire, strip]
                    
                # Legend:
                # 5 = 'SUMO5 Backward', 15 = 'SUMO5 Forward'
                # 1 = sectors number, always 1 for SUMO5
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    5, 1, md + 1, segment + 1, wire + 1, strip + 1, 2
                )

                tempF = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    15, 1, md + 1, segment + 1, wire + 1, strip + 1, 2
                )

                # Legend: x,y,z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    VX[md_segt_id, wire, strip],
                    VY[md_segt_id, wire, strip],
                    VZ[md_segt_id, wire, strip]
                )

                tempF1 = '%.2f\t%.2f\t%.2f' % (
                    VXF[md_segt_id, wire, strip],
                    VYF[md_segt_id, wire, strip],
                    VZF[md_segt_id, wire, strip]
                )

                # Legend:
                # voxel dimensions to be used to generate Nexus 
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    -np.deg2rad(wire * dthetaS5)*0,
                    GLx1bbS5[wire, strip],
                    GLx2bbS5[wire, strip],
                    2 * GLy1S5[wire, strip],
                    2 * GLy2S5[wire, strip],
                    2 * GLzS5[wire, strip]
                )

                # Legend:
                # rotation angles to put the voxels into the right positions, Backward EndCap
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS5, angM, 0, 0, 0, 0
                )

                # rotation angles to put the voxels into the right positions, Forward EndCap
                tempF3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS5, angM, 180, 90, 0, 0
                )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                stringaf = tempF + '\t' + tempF1 + '\t' + temp2 + '\t' + tempF3

                fB.writelines(stringa)
                fF.writelines(stringaf)

fB.close()
fF.close()
