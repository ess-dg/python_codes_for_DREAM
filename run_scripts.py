#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wednesday Mar 16 16:49:10 2022

@author: irinastefanescu
"""

import os
import datetime
import globals

now = datetime.datetime.now()
print(' ')
print("Current date and time: ")
print(now.strftime("%Y-%m-%d %H:%M:%S \n"))

# print the number of EndCap sectors
print(' ')
print('Number of EndCap sectors considered in the calculation: ',
      globals.no_modules) 
# print the number of Mantle modules
print('Number of Mantle modules considered in the calculation: ',
      globals.no_modulesM)

out_file = 'DREAMAll_voxels.txt'
tempM_file = 'Mantle_temp.txt'
tempF_file = 'Forward_temp.txt'
tempB_file = 'Backward_temp.txt'
tempHR_file = 'HR_temp.txt'


try:
    os.remove(out_file)
    print(f"\nOld output file {out_file} removed.\n")
except FileNotFoundError:
    print(f"Old file {out_file} not found. Creating it during this run.\n")

try:
    os.remove(tempF_file)
    print(f"\nOld output file {tempF_file} removed.\n")
except FileNotFoundError:
    print(f"Old file {tempF_file} not found. Creating it during this run.\n")

try:
    os.remove(tempM_file)
    print(f"\nOld output file {tempM_file} removed.\n")
except FileNotFoundError:
    print(f"Old file {tempM_file} not found. Creating it during this run.\n")

try:
    os.remove(tempB_file)
    print(f"\nOld output file {tempB_file} removed.\n")
except FileNotFoundError:
    print(f"Old file {tempB_file} not found. Creating it during this run.\n")

try:
    os.remove(tempHR_file)
    print(f"\nOld output file {tempHR_file} removed.\n")
except FileNotFoundError:
    print(f"Old file {tempHR_file} not found. Creating it during this run.\n")

os.system(f'python DREAMMantle_calculate_voxels.py {tempM_file}')
print('Mantle done!\n')

os.system(f'python DREAMHR_calculate_voxels.py {tempHR_file}')
print('High-Resolution done!\n')

os.system('python DREAMSUMO3_calculate_voxels.py')
print('SUMO3 Backward & SUMO3 Forward done!\n')

os.system('python DREAMSUMO4_calculate_voxels.py')
print('SUMO4 Backward & SUMO4 Forward done!\n')

os.system('python DREAMSUMO5_calculate_voxels.py')
print('SUMO5 Backward & SUMO5 Forward done!\n')

os.system('python DREAMSUMO6_calculate_voxels.py')
print('SUMO6 Backward & SUMO6 Forward done!\n')

os.system(f'cat {tempM_file} {tempHR_file} {tempB_file} {tempF_file} > {out_file}')
#os.system(f'cat {tempB_file} {tempF_file} > {out_file}')
