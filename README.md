# python_codes_for_DREAM
contains python scripts to generate the voxel table for the DREAM detectors

# Codes written by Irina Stefanescu, ESS NSS Detector Group and Celine Durniak, ESS DMSC DRAM Group
irina.stefanescu@ess.eu
celine.durniak@ess.eu

## Description

**To be completed**
One Python script per detectors' part ,i.e., SUMO3, 4, 5 ,6, Mantle, High-resolution 

File lists:

- DREAMHR_calculate_voxels.py
- DREAMMantle_calculate_voxels.py
- DREAMSUMO3_calculate_voxels.py
- DREAMSUMO4_calculate_voxels.py
- DREAMSUMO5_calculate_voxels.py
- DREAMSUMO6_calculate_voxels.py
- dream_v7.py
- globals.py 
- run_scripts.py


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
