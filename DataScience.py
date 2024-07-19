# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:18:51 2024

@author: YAMEOGO Philippe
"""

import pandas as pd
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Load the NetCDF4 file
data = "data.nc4"
dat = xr.open_dataset(data)



help(xr.open_dataset)


# Print the dataset structure to check dimensions
print(dat)

# Access the precipitation variable
precipitation = dat['precipitation']
print(precipitation)