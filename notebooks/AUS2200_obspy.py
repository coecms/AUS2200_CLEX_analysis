import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os

# ---- Input ---- #
pathData    = '/scratch/public/pjs548/aus2200_v2/obs_int/' # Path to interpolated data
dateStart   = '20220222T0600Z' # Start date in format yyyymmThhmmZ
dateEnd     = '' # End date in format yyyymmThhmmz (included)
projectName = 'aus2200' # Project name (first part of date filenames)
varModel = ''
varObs   = ''

# Find required files and concatenate variables
filenames = [f"{projectName}_{startDate}.nc"]
# for file in os.listdir(pathData):
#     if(file.startswith(projectName and file.endswith('.nc')):
nc = xr.open_dataset(filenames[0])
lat = nc['lat']
lon = nc['lon']