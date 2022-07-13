## python_codes_for_DREAM
contains the python scripts to generate the voxel table for the DREAM detectors

Scripts written by Irina Stefanescu, ESS NSS Detector Group and Celine Durniak, ESS DMSC DRAM Group

irina.stefanescu@ess.eu

celine.durniak@ess.eu

## Description

**To be completed**

One Python script per detectors' part ,i.e., EndCap FORWARD (SUMO3, SUMO4, SUMO5, SUMO6), EndCap BACKWARD (SUMO3, SUMO4, SUMO5, SUMO6), Mantle, High-resolution.  

File lists:

- DREAMHR_calculate_voxels.py
- DREAMMantle_calculate_voxels.py
- DREAMSUMO3_calculate_voxels.py
- DREAMSUMO4_calculate_voxels.py
- DREAMSUMO5_calculate_voxels.py
- DREAMSUMO6_calculate_voxels.py
- dream.py
- globals.py 
- run_scripts.py

The result of the calculation with *run_scripts.py* is a number of txt files (one per detector (sub)system) containing the information on the location (x, y, z coordinates with respect to the sample position at (0, 0, 0)) and shape parameters (trapezoidal) of each individual detector voxel along with some hardware information (wire number, strip number, segment number, module number) that will make it possible to match the calculated detector voxels to the real ones (when available). The number of Mantle and EndCap detector modules included in the calculation can be controled from the *globals.py* file. The file *DREAMAll_voxels.txt* is obtained through the concatanation of the (sub)system files and it is used by the script *dream.py* to generate the off- and nxs-files of the DREAM detector. 

## Installation and usage

- Create and activate a virtual environment in the folder containing the scripts using Python 3 
  (tested with Python 3.9). In the following we assume that `python` refers to the version of 
  Python you want to use.

   ```
   python -m venv .venv
   source .venv/bin/activate
   ```
   Executing this last command should change the prompt of your terminal. It should now starts with `(.venv)`. 

- Upgrade pip (optional)

  ```
  python -m pip install --upgrade pip
  ```

- Install the required libraries

  ```
  python -m pip install -r requirements.txt
  ```

- Generate the tables of voxels for all detector sub-systems for DREAM.

  ```
  python run_scripts.py
  ```
  You can comment out the tables not needed for your studies before running the above command.
  
- To deactivate the virtual environment, simply type `deactivate` in the terminal. The prompt should change back to 
  its initial state. And if you do not need this environment, you can simply delete the `.venv` folder.0
