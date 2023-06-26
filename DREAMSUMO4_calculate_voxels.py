#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Mar 12 10:13:30 2022

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
 DREAM SUMO4 module. A Mantle detector segment consists
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

hfS4 = 157       # segment height front, sample side
hbS4 = 199       # segment height back, readout side
sensS4 = 315     # segment depth
wbfS4 = 93       # segment width bottom front, sample side
wtfS4 = 124      # segment width top front, sample side
wbbS4 = 107      # segment width bottom back, readout side
wtbS4 = 147      # segment width top back, readout side

margin = 10      # empty space inside the segment front-back and left-right
marginb = 1.5    # empty space inside the segment top-bottom
shieldz = 25     # length of the shielding block at the back of the segment

Althick = 0.3     # thickness Al cathode material and housing
Bthick = 0.0011   # thickness Boron coating

# start position for placing the modules in the frame, integer number
# to multiply with 12*deg
index_rot = 10

tilt_theta = -10       # tilt_angle in deg
tilt_phiS4 = 16        # inclination angle module in deg
dphi = np.deg2rad(12)  # angular coverage in phi for the whole module
s_angle_S4 = 0.35      # rotation segment around the Y axis
s_offset_S4 = 0.17     # additional angle for the Z rotation of the segment
alpha_physS4 = 7.64    # I don't remember what this is :-(

n_wires = 16    # no of wires
n_strips = 16   # no of strips
nS4 = 6         # no of segments per module

offsetX_S4 = 30.       # x-distance SUMO4 from the sample
offsetY_S4 = 592.      # y-distance SUMO4 from the sample
offsetZ_S4 = -1310.    # z-distance SUMO4 from the sample

# start calculations 

# this function recalculates the shape parameters in order
# to avoid the crash of the G4Trap class
# I must have it here too, for debugging purposes

Aly1S4 = hfS4 / 2   # 1/2 height Al housing at entrance window SUMO4
Aly2S4 = hbS4 / 2   # 1/2 height Al housing at the back SUMO4
AlzS4 = (sensS4 + margin) / 2    # depth housing SUMO4
  
# subtract (no_of_segments - 1)*mm from the widths to account for the space between the segments  

# 1/2 width Al housing at the bottom front in SUMO4
Alx1bS4 = (wbfS4 - 12) / 2 / nS4  
# 1/2 width segment at the top front in SUMO4
Alx1tS4 = (wtfS4 - 10) / 2 / nS4 
# 1/2 width segment at the bottom back SUMO4
Alx2bS4 = (wbbS4 - 12) / 2 / nS4 
# 1/2 width segment at the top back SUMO4
Alx2tS4 = (wtbS4 - 23) / 2 / nS4 

# recalculated shape parameters, same names
Aly1S4, Alx1bS4, Alx1tS4, Aly2S4, Alx2bS4, Alx2tS4 = \
    globals.match2geantvars([Aly1S4, Alx1bS4, Alx1tS4, Aly2S4, Alx2bS4, Alx2tS4])

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

By1S4 = Aly1S4 - Althick / np.cos(dphi/2)
By2S4 = Aly2S4 - Althick / np.cos(dphi/2)
BzS4 = AlzS4 - Althick
Bx1tS4 = Alx1tS4 - Althick / np.cos(dphi/2)
Bx2tS4 = Alx2tS4 - Althick / np.cos(dphi/2)
Bx1bS4 = Alx1bS4 - Althick / np.cos(dphi/2)
Bx2bS4 = Alx2bS4 - Althick / np.cos(dphi/2)

# recalculated shape parameters, same names

By1S4, Bx1bS4, Bx1tS4, By2S4, Bx2bS4, Bx2tS4 = \
    globals.match2geantvars([By1S4, Bx1bS4, Bx1tS4, By2S4, Bx2bS4, Bx2tS4])

# Gas volume available to the gas voxels

Gy1S4 = By1S4
Gy2S4 = By2S4
GzS4 = BzS4 - Bthick
Gx1tS4 = Bx1tS4 - Bthick / np.cos(dphi/2)
Gx2tS4 = Bx2tS4 - Bthick / np.cos(dphi/2)
Gx1bS4 = Bx1bS4 - Bthick / np.cos(dphi/2)
Gx2bS4 = Bx2bS4 - Bthick / np.cos(dphi/2)

Gy1S4, Gx1bS4, Gx1tS4, Gy2S4, Gx2bS4, Gx2tS4 = \
    globals.match2geantvars([Gy1S4, Gx1bS4, Gx1tS4, Gy2S4, Gx2bS4, Gx2tS4])

CathSubstrY1S4 = Gy1S4 - marginb / 2
CathSubstrY2S4 = Gy2S4 - marginb / 2
CathSubstrZS4 = BzS4 - margin / 2

CathConvY1S4 = CathSubstrY1S4
CathConvY2S4 = CathSubstrY2S4
CathConvZS4 = CathSubstrZS4

# calculate the dimensions of the gas voxels

xx1tS4 = Gx1tS4 - 2
xx2tS4 = Gx2tS4 - 2
xx1bS4 = Gx1bS4 - 2
xx2bS4 = Gx2bS4 - 2

eta_b = 2 * (xx1tS4 - xx1bS4) / Gy1S4
eta_t = (xx2bS4 - xx1bS4) / sensS4
eta_w = (Gy2S4 - Gy1S4) / sensS4    # Gy1 and Gy2 already halved

izzS4 = sensS4 / n_strips         # strip pitch SUMO4 in mm, all equal

dthetaS4 = alpha_physS4 / n_wires  # wire pitch SUMO4 in deg, all equal

shp = (n_strips, n_wires)

GLzS4 = np.zeros(shp)
GLy1S4 = np.zeros(shp)
GLy2S4 = np.zeros(shp)
GLx1bS4 = np.zeros(shp)
GLx2bS4 = np.zeros(shp)
GLx1bbS4 = np.zeros(shp)
GLx2bbS4 = np.zeros(shp)

# loop goes from 0 to max-1!!
# wires run parallel to the beam axis and in fan out geometry
# wire pitch smaller at the segment front and larger at the back
# also, wire pitch is symmetric with respect to the segment center
# the calculations starts from the center of the segment, 
# so wire=1 is close to the center, wire = 8 is at the bottom (or top)

for strip in range(n_strips):  # loop over strips

    for wire in range(n_wires//2):  # loop over all 1/2*wires    
        # voxel depth
        GLzS4[wire, strip] = izzS4 / 2
        GLzS4[n_wires//2 + wire, strip] = izzS4 / 2
        # wire pitch front of the voxel
        # GLy1S4[wire, strip] = (strip * 2 * izzS4 * eta_w + 2 * (Gy1S4 - 4)) / 2 / n_wires
        GLy1S4[wire, strip] = (strip * izzS4 * eta_w + Gy1S4 - 4) / n_wires
        # wire pitch, back of the voxel
        # GLy2S4[wire, strip] = ((strip + 1) * 2 * izzS4 * eta_w + 2 * (Gy1S4 - 4)) / 2 / n_wires
        GLy2S4[wire, strip] = ((strip + 1) * izzS4 * eta_w + Gy1S4 - 4) / n_wires
        GLy1S4[n_wires//2 + wire, strip] = GLy1S4[wire, strip]
        GLy2S4[n_wires//2 + wire, strip] = GLy2S4[wire, strip] 

        # for the bottom 8 wires

        GLx1bbS4[wire, strip] = 0.5 * (2 * xx1bS4 +
                                       wire * GLy1S4[wire, strip] * eta_b)
        GLx2bbS4[wire, strip] = GLx1bbS4[wire, strip] + 0.5 * strip * izzS4 * eta_t
        # for the upper 8 wires

        GLx1bS4[wire, strip] = 0.5 * (2 * xx1bS4 +
                                      (n_wires//2 + wire) *
                                      GLy1S4[wire, strip] * eta_b)
        GLx2bS4[wire, strip] = GLx1bS4[wire, strip] + 0.5 * strip * izzS4 * eta_t
        GLx1bbS4[n_wires//2 + wire, strip] = GLx1bS4[wire, strip]
        GLx2bbS4[n_wires//2 + wire, strip] = GLx2bS4[wire, strip]

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
        # voxelYY[wire, strip] = -(n_wires//2 - wire - 0.5) * 2 * GLy2S4[n_wires//2 - wire, strip]
        voxelYY[wire, strip] = -(n_wires - 2 * wire - 1.) * GLy2S4[n_wires // 2 - wire, strip]
        voxelZZ[wire, strip] = \
            (n_strips - strip - 1) * 2 * GLzS4[n_wires//2 - wire, 0] \
            + GLzS4[n_wires//2 - wire, 0] + margin/2
        voxelXX[wire, strip] = Althick
        voxelXXc[wire, strip] = 0.5 * GLx1bbS4[wire, strip]
        # voxels created by the lowest 8 wires
        # fill the voxels from the bottom to the segment center
        # voxelYY[n_wires//2 + wire, strip] = (wire + 0.5) * 2 * GLy2S4[wire, strip]
        voxelYY[n_wires // 2 + wire, strip] = (2 * wire + 1.) * GLy2S4[wire, strip]
        voxelZZ[n_wires//2 + wire, strip] = \
            (n_strips - strip - 1) * 2 * GLzS4[wire, 0] \
            + GLzS4[wire, 0] + margin/2
        voxelXX[n_wires//2 + wire, strip] = Althick
        voxelXXc[n_wires//2 + wire, strip] = 0.5 * GLx1bS4[wire, strip]


# calculate the segment positions in the detector frame
segX = np.zeros(nS4)
segZ = np.zeros(nS4)
segY = np.zeros(nS4)

for seg_no in range(nS4):
    segX[seg_no] = Alx2tS4 * (nS4 + 1 - 2 * (seg_no + 1))
    segZ[seg_no] = -sensS4 / 2
    
radS4 = np.sqrt(np.power(offsetX_S4, 2) + np.power(offsetY_S4, 2))

no_modules = globals.no_modules

modX = np.zeros(no_modules)
modY = np.zeros(no_modules)
modZ = np.zeros(no_modules)
modZ[:] = offsetZ_S4

for mod in range(no_modules):
    # modX[mod] = radS4 * np.sin(np.deg2rad(-(index_rot + mod) * np.rad2deg(dphi)))
    # modY[mod] = radS4 * np.cos(np.deg2rad(-(index_rot + mod) * np.rad2deg(dphi)))
    modX[mod] = radS4 * np.sin(-(index_rot + mod) * dphi)
    modY[mod] = radS4 * np.cos(-(index_rot + mod) * dphi)

""" calculate the lookup table  """

shp = (100 * no_modules + nS4, n_wires, n_strips)

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

mX_s = np.sin(np.deg2rad(tilt_phiS4))
mX_c = np.cos(np.deg2rad(tilt_phiS4))

fY_s = np.sin(np.deg2rad(180))
fY_c = np.cos(np.deg2rad(180))

fZ_s = np.sin(np.deg2rad(90))
fZ_c = np.cos(np.deg2rad(90))

# voxels in the left counter

for md in range(no_modules):
    
    angM = (index_rot + md) * np.rad2deg(dphi)
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nS4):

        angZ = (nS4 + 1 - 2 * (segment + 1)) * (s_angle_S4 + s_offset_S4)
        angY = (nS4 + 1 - 2 * (segment + 1)) * s_angle_S4
        segZ_s = np.sin(np.deg2rad(angZ))
        segZ_c = np.cos(np.deg2rad(angZ))
        segY_s = np.sin(np.deg2rad(angY))
        segY_c = np.cos(np.deg2rad(angY))

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

                sz_z[md_segt_id, wire, strip] = segZ[segment] + voxelZZ[wire, strip]

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
                VX[md_segt_id, wire, strip] = modX[md] + mx[md_segt_id, wire, strip]
                VY[md_segt_id, wire, strip] = modY[md] + my[md_segt_id, wire, strip]
                VZ[md_segt_id, wire, strip] = modZ[md] + mz[md_segt_id, wire, strip]
                
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
                # 4 = 'SUMO4 Backward', 14 = 'SUMO4 Forward'
                # 1 = sectors number, always 1 for SUMO4 & Forward
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    4, 1, md + 1, segment + 1, wire + 1, strip + 1, 1
                )

                tempF = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    14, 1, md + 1, segment + 1, wire + 1, strip + 1, 1
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
                    np.deg2rad(wire*dthetaS4)*0,
                    GLx1bbS4[wire, strip],
                    GLx2bbS4[wire, strip],
                    2*GLy1S4[wire, strip],
                    2*GLy2S4[wire, strip],
                    2*GLzS4[wire, strip])

                # Legend:
                # rotation angles to put the voxels into the right positions, Backward EndCap
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS4, angM, 0, 0, 0, 0
                    )

                # rotation angles to put the voxels into the right positions, Forward EndCap
                tempF3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS4, angM, 180, 90, 0, 0
                    )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3

                stringaf = tempF + '\t' + tempF1 + '\t' + \
                    temp2 + '\t' + tempF3

                fB.writelines(stringa)
                fF.writelines(stringaf)

# voxels in the right counter

for md in range(no_modules):
    
    angM = (index_rot + md) * np.rad2deg(dphi)
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nS4):

        angZ = (nS4 + 1 - 2 * (segment + 1)) * (s_angle_S4 + s_offset_S4)
        angY = (nS4 + 1 - 2 * (segment + 1)) * s_angle_S4
        segZ_s = np.sin(np.deg2rad(angZ))
        segZ_c = np.cos(np.deg2rad(angZ))
        segY_s = np.sin(np.deg2rad(angY))
        segY_c = np.cos(np.deg2rad(angY))

        md_segt_id = 100 * md + segment

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

                sz_z[md_segt_id, wire, strip] = segZ[segment] + voxelZZ[wire, strip]

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
                # 4 = 'SUMO4 Backward', 14 = 'SUMO4 Forward'
                # 1 = sectors number, always 1 for SUMO4
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    4, 1, md + 1, segment + 1, wire + 1, strip + 1, 2
                )
                tempF = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    14, 1, md + 1, segment + 1, wire + 1, strip + 1, 2
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
                    -np.deg2rad(wire * dthetaS4) * 0,
                    GLx1bbS4[wire, strip],
                    GLx2bbS4[wire, strip],
                    2*GLy1S4[wire, strip],
                    2*GLy2S4[wire, strip],
                    2*GLzS4[wire, strip])

                # Legend:
                # rotation angles to put the voxels into the right positions, Backward EndCap
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS4, angM, 0, 0, 0, 0
                    )

                # rotation angles to put the voxels into the right positions, Forward EndCap
                tempF3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS4, angM, 180, 90, 0, 0
                    )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                stringaf = tempF + '\t' + tempF1 + '\t' + temp2 + '\t' + tempF3

                fB.writelines(stringa)
                fF.writelines(stringaf)

fB.close()
fF.close()
