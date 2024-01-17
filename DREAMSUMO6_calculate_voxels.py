#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 10:13:30 2022

@author: irinastefanescu
"""

import numpy as np
import globals


np.set_printoptions(precision=4)

"""
*******The SUMO6 for EndCap ******

 Uses the engineering specifications from the company CDT and
 generates the centers and dimensions of the voxels for the
 DREAM SUMO6 module. A Mantle detector segment consists
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

hfS6 = 165      # segment height front, sample side
hbS6 = 214      # segment height back, readout side
sensS6 = 350    # segment depth
wbfS6 = 171     # segment width bottom front, sample side
wtfS6 = 201     # segment width top front, sample side
wbbS6 = 204     # segment width bottom back, readout side
wtbS6 = 243     # segment width top back, readout side

margin = 10      # empty space inside the segment front-back and left-right
marginb = 1.5    # empty space inside the segment top-bottom
shieldz = 25     # length of the shielding block at the back of the segment

Althick = 0.3       # thickness Al cathode material and housing
Bthick = 0.0011     # thickness Boron coating

# start position for placing the modules in the frame, integer number
# to multiply with 12*deg
index_rot = globals.index_rot

tilt_theta = -10       # tilt_angle in deg
tilt_phiS6 = 30        # inclination angle module in deg
dphi = np.deg2rad(12)  # angular coverage in phi for the whole module
s_angle_S6 = 0.39      # rotation segment around the Y axis
s_offset_S6 = 0.       # additional angle for the Z rotation of the segment
alpha_physS6 = 7.96    # I don't remember what this is :-(

n_wires = 16     # no of wires
n_strips = 16    # no of strips
nS6 = 10         # no of segments per module

offsetX_S6 = 50.      # x-distance SUMO6 from the sample
offsetY_S6 = 1011.    # y-distance SUMO6 from the sample
offsetZ_S6 = -1340.   # z-distance SUMO6 from the sample

# start calculations 

# this function recalculates the shape parameters in order
# to avoid the crash of the G4Trap class
# I must have it here too, for debugging purposes

Aly1S6 = hfS6 / 2   # 1/2 height Al housing at entrance window SUMO6
Aly2S6 = hbS6 / 2   # 1/2 height Al housing at the back SUMO6
AlzS6 = (sensS6 + margin) / 2    # depth housing SUMO6
  
# subtract (no_of_segments - 1)*mm from the widths to account for the space between the segments  

# 1/2 width Al housing at the bottom front in SUMO6
Alx1bS6 = (wbfS6 - 17) / 2 / nS6  
# 1/2 width segment at the top front in SUMO6
Alx1tS6 = (wtfS6 - 17) / 2 / nS6 
#  1/2 width segment at the bottom back SUMO6
Alx2bS6 = (wbbS6 - 28) / 2 / nS6 
#  1/2 width segment at the top back SUMO6
Alx2tS6 = (wtbS6 - 12) / 2 / nS6 

# recalculated shape parameters, same names
Aly1S6, Alx1bS6, Alx1tS6, Aly2S6, Alx2bS6, Alx2tS6 = \
    globals.match2geantvars([Aly1S6, Alx1bS6, Alx1tS6, Aly2S6, Alx2bS6, Alx2tS6])

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

By1S6 = Aly1S6 - Althick / np.cos(dphi/2)
By2S6 = Aly2S6 - Althick / np.cos(dphi/2)
BzS6 = AlzS6 - Althick
Bx1tS6 = Alx1tS6 - Althick / np.cos(dphi/2)
Bx2tS6 = Alx2tS6 - Althick / np.cos(dphi/2)
Bx1bS6 = Alx1bS6 - Althick / np.cos(dphi/2)
Bx2bS6 = Alx2bS6 - Althick / np.cos(dphi/2)

# recalculated shape parameters, same names

By1S6, Bx1bS6, Bx1tS6, By2S6, Bx2bS6, Bx2tS6 = \
    globals.match2geantvars([By1S6, Bx1bS6, Bx1tS6, By2S6, Bx2bS6, Bx2tS6])

# Gas volume available to the gas voxels

Gy1S6 = By1S6
Gy2S6 = By2S6
GzS6 = BzS6 - Bthick
Gx1tS6 = Bx1tS6 - Bthick / np.cos(dphi/2)
Gx2tS6 = Bx2tS6 - Bthick / np.cos(dphi/2)
Gx1bS6 = Bx1bS6 - Bthick / np.cos(dphi/2)
Gx2bS6 = Bx2bS6 - Bthick / np.cos(dphi/2)

Gy1S6, Gx1bS6, Gx1tS6, Gy2S6, Gx2bS6, Gx2tS6 = \
    globals.match2geantvars([Gy1S6, Gx1bS6, Gx1tS6, Gy2S6, Gx2bS6, Gx2tS6])

CathSubstrY1S6 = Gy1S6 - marginb / 2
CathSubstrY2S6 = Gy2S6 - marginb / 2
CathSubstrZS6 = BzS6 - margin / 2

CathConvY1S6 = CathSubstrY1S6
CathConvY2S6 = CathSubstrY2S6
CathConvZS6 = CathSubstrZS6

# calculate the dimensions of the gas voxels

xx1tS6 = Gx1tS6 - 2
xx2tS6 = Gx2tS6 - 2
xx1bS6 = Gx1bS6 - 2
xx2bS6 = Gx2bS6 - 2

eta_b = 2 * (xx1tS6 - xx1bS6) / Gy1S6
eta_t = (xx2bS6 - xx1bS6) / sensS6
eta_w = (Gy2S6 - Gy1S6) / sensS6   # Gy1 and Gy2 already halved

izzS6 = sensS6 / n_strips         # strip pitch SUMO6 in mm, all equal

dthetaS6 = alpha_physS6 / n_wires  # wire pitch SUMO6 in deg, all equal

shp = (n_strips, n_wires)

GLzS6 = np.zeros(shp)
GLy1S6 = np.zeros(shp)
GLy2S6 = np.zeros(shp)
GLx1bS6 = np.zeros(shp)
GLx2bS6 = np.zeros(shp)
GLx1bbS6 = np.zeros(shp)
GLx2bbS6 = np.zeros(shp)

# loop goes from 0 to max-1!!
# wires run parallel to the beam axis and in fan out geometry
# wire pitch smaller at the segment front and larger at the back
# also, wire pitch is symmetric with respect to the segment center
# the calculations starts from the center of the segment, 
# so wire=1 is close to the center, wire = 8 is at the bottom (or top)

for strip in range(n_strips):  # loop over strips

    for wire in range(n_wires//2):  # loop over all 1/2*wires    
        # voxel depth
        GLzS6[wire, strip] = izzS6 / 2
        GLzS6[n_wires//2 + wire, strip] = izzS6 / 2
        # wire pitch front of the voxel
        GLy1S6[wire, strip] = (strip * 2 * izzS6 * eta_w +
                               2 * (Gy1S6 - 4)) / 2 / n_wires
        # wire pitch, back of the voxel
        GLy2S6[wire, strip] = ((strip + 1) * 2 * izzS6 * eta_w +
                               2 * (Gy1S6 - 4)) / 2 / n_wires
        GLy1S6[n_wires//2 + wire, strip] = GLy1S6[wire, strip]
        GLy2S6[n_wires//2 + wire, strip] = GLy2S6[wire, strip] 

        # for the bottom 8 wires

        GLx1bbS6[wire, strip] = 0.5 * (2 * xx1bS6 +
                                       wire * GLy1S6[wire, strip] * eta_b)
        GLx2bbS6[wire, strip] = GLx1bbS6[wire, strip] + 0.5 * strip * izzS6 * eta_t

        # for the upper 8 wires

        GLx1bS6[wire, strip] = 0.5 * (2 * xx1bS6 +
                                      (n_wires//2 + wire) *
                                      GLy1S6[wire, strip] * eta_b)
        GLx2bS6[wire, strip] = GLx1bS6[wire, strip] + 0.5 * strip * izzS6 * eta_t
        GLx1bbS6[n_wires//2 + wire, strip] = GLx1bS6[wire, strip]
        GLx2bbS6[n_wires//2 + wire, strip] = GLx2bS6[wire, strip]

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
        # voxelYY[wire, strip] = -(n_wires//2 - wire - 0.5) * \
        #                         2 * GLy2S6[n_wires//2 - wire, strip]
        voxelYY[wire, strip] = -(n_wires - 2 * wire - 1) * GLy2S6[n_wires // 2 - wire, strip]
        voxelZZ[wire, strip] = \
            (n_strips - strip - 1) * 2 * GLzS6[n_wires//2 - wire, 0] \
            + GLzS6[n_wires//2 - wire, 0] + margin/2
        voxelXX[wire, strip] = Althick
        voxelXXc[wire, strip] = 0.5 * GLx1bbS6[wire, strip]
        # voxels created by the lowest 8 wires
        # fill the voxels from the bottom to the segment center
        voxelYY[n_wires//2 + wire, strip] = (wire + 0.5) * 2 * GLy2S6[wire, strip]
        voxelZZ[n_wires//2 + wire, strip] = \
            (n_strips - strip - 1) * 2 * GLzS6[wire, 0] + GLzS6[wire, 0] + margin/2
        voxelXX[n_wires//2 + wire, strip] = Althick
        voxelXXc[n_wires//2 + wire, strip] = 0.5 * GLx1bS6[wire, strip]

# calculate the segment positions in the detector frame
segX = np.zeros(nS6)
segZ = np.zeros(nS6)
segY = np.zeros(nS6)

for seg_no in range(nS6):
    segX[seg_no] = Alx2tS6 * (nS6 + 1 - 2 * (seg_no + 1))
    segZ[seg_no] = -sensS6 / 2
    
radS6 = np.sqrt(np.power(offsetX_S6, 2) + np.power(offsetY_S6, 2))

no_modules = globals.no_modules

modX = np.zeros(no_modules)
modY = np.zeros(no_modules)
modZ = np.zeros(no_modules)
modZ[:] = offsetZ_S6

for mod in range(no_modules):
    # modX[mod] = radS6 * np.sin(np.deg2rad(-(index_rot + mod) * np.rad2deg(dphi)))
    # modY[mod] = radS6 * np.cos(np.deg2rad(-(index_rot + mod) * np.rad2deg(dphi)))
    modX[mod] = radS6 * np.sin(-(index_rot + mod) * dphi)
    modY[mod] = radS6 * np.cos(-(index_rot + mod) * dphi)

""" calculate the lookup table  """

shp = (100 * no_modules + nS6, n_wires, n_strips)

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

mX_s = np.sin(np.deg2rad(tilt_phiS6))
mX_c = np.cos(np.deg2rad(tilt_phiS6))

fY_s = np.sin(np.deg2rad(180))
fY_c = np.cos(np.deg2rad(180))

fZ_s = np.sin(np.deg2rad(90))
fZ_c = np.cos(np.deg2rad(90))

# voxels in the left counter

for md in range(no_modules):
    
    angM = (index_rot + md) * np.rad2deg(dphi)
    mZ_s = np.sin(np.deg2rad(angM))
    mZ_c = np.cos(np.deg2rad(angM))

    for segment in range(nS6):

        angZ = ((nS6 + 1) - 2 * (segment + 1)) * (s_angle_S6 + s_offset_S6)
        angY = ((nS6 + 1) - 2 * (segment + 1)) * s_angle_S6
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
                    (segX[segment] - voxelXX[wire, strip] - voxelXXc[wire, strip]) + \
                    voxelYY[wire, strip] * segZ_s

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

                # Forward detector: mirror reflection of Backward % x,y plane
                VXF[md_segt_id, wire, strip] = VX[md_segt_id, wire, strip]
                VYF[md_segt_id, wire, strip] = VY[md_segt_id, wire, strip]
                VZF[md_segt_id, wire, strip] = -VZ[md_segt_id, wire, strip]

                # Legend:
                # 6 = 'SUMO6 Backward', 16 = 'SUMO6 Forward'
                # 1 = sectors number, always 1 for SUMO6 & Forward
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    6, 1, md + 1, segment + 1, wire + 1, strip + 1, 1
                )

                tempF = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    16, 1, md + 1, segment + 1, 
                    wire + 1, strip + 1, 1
                )

                # Legend: x,y,z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    VX[md_segt_id, wire, strip],
                    VY[md_segt_id, wire, strip],
                    VZ[md_segt_id, wire, strip]
                )

               # Legend: x,y,z voxel centers
                tempF1 = '%.2f\t%.2f\t%.2f' % (
                    VXF[md_segt_id, wire, strip],
                    VYF[md_segt_id, wire, strip],
                    VZF[md_segt_id, wire, strip]
                )

                # Legend:
                # voxel dimensions to be used to generate Nexus 
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    np.deg2rad(wire*dthetaS6)*0,
                    GLx1bbS6[wire, strip],
                    GLx2bbS6[wire, strip],
                    2 * GLy1S6[wire, strip],
                    2 * GLy2S6[wire, strip],
                    2 * GLzS6[wire, strip])

                # Legend:
                # rotation angles to put the voxels into the right positions, Backward EndCap
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS6, angM, 0, 0, 0, 0
                    )

                # rotation angles to put the voxels into the right positions, Forward EndCap
                tempF3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS6, angM, 180, 90, 0, 0
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

    for segment in range(nS6):

        md_segt_id = 100 * md + segment

        angZ = (nS6 + 1 - 2 * (segment + 1)) * (s_angle_S6 + s_offset_S6)
        angY = (nS6 + 1 - 2 * (segment + 1)) * s_angle_S6
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

                # Forward detector: mirror reflection of Backward % x,y plane
                VXF[md_segt_id, wire, strip] = VX[md_segt_id, wire, strip]
                VYF[md_segt_id, wire, strip] = VY[md_segt_id, wire, strip]
                VZF[md_segt_id, wire, strip] = -VZ[md_segt_id, wire, strip]

                # Legend:
                # 6 = 'SUMO6 Backward', 16 = 'SUMO6 Forward'
                # 1 = sectors number, always 1 for SUMO6
                #     only relevant for the HR detector
                # module no, segment no, wire no, strip no, counter no 
                temp = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    6, 1, md + 1, segment + 1, wire + 1, strip + 1, 2
                )

                tempF = '%d\t%d\t%d\t%d\t%d\t%d\t%d' % (
                    16, 1, md + 1, segment + 1, wire + 1, strip + 1, 2
                )

                # Legend: x,y,z voxel centers
                temp1 = '%.2f\t%.2f\t%.2f' % (
                    VX[md_segt_id, wire, strip],
                    VY[md_segt_id, wire, strip],
                    VZ[md_segt_id, wire, strip]
                )

               # Legend: x,y,z voxel centers
                tempF1 = '%.2f\t%.2f\t%.2f' % (
                    VXF[md_segt_id, wire, strip],
                    VYF[md_segt_id, wire, strip],
                    VZF[md_segt_id, wire, strip]
                )

                # Legend:
                # voxel dimensions to be used to generate Nexus 
                temp2 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (
                    -np.deg2rad(wire*dthetaS6)*0,
                    GLx1bbS6[wire, strip],
                    GLx2bbS6[wire, strip],
                    2 * GLy1S6[wire, strip],
                    2 * GLy2S6[wire, strip],
                    2 * GLzS6[wire, strip])

                # Legend:
                # rotation angles to put the voxels into the right positions, Backward EndCap
                temp3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS6, angM, 0, 0, 0, 0
                )

                # rotation angles to put the voxels into the right positions, Forward EndCap
                tempF3 = '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n' % (
                    -angY, angZ, -tilt_theta, tilt_phiS6, angM, 180, 90, 0, 0
                )

                stringa = temp + '\t' + temp1 + '\t' + temp2 + '\t' + temp3
                stringaf = tempF + '\t' + tempF1 + '\t' + temp2 + '\t' + tempF3

                fB.writelines(stringa)
                fF.writelines(stringaf)

fB.close()
fF.close()
