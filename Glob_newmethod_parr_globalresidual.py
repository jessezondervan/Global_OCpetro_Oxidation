#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from time import time
start = time()

##Load the required packages
#Visualization
import rasterio as rio                                                                                                                  
from rasterio.plot import show  
import rasterio.mask
from rasterio.rio import stack
from rasterio.rio.stack import stack
import matplotlib.pyplot as plt                                                                                                         

#Working with Shapefiles
import numpy as np
import fiona

#Projections and Transformations
from pyproj import Transformer

#Others
import gc
import pandas as pd
import os
import glob
import scipy.interpolate
import matplotlib as mpl

#Parallelisation
from concurrent.futures import ProcessPoolExecutor

#computation
import numba
from numba import jit
from numba import float32, types
from numba.experimental import jitclass

import csv

import os

quantile = float(os.getenv('SLURM_ARRAY_TASK_ID'))/100.

glob_denud_tif = rasterio.open('input_global/glob_den_lor3.tif')
profile = glob_denud_tif.profile
del glob_denud_tif

problist = np.array([0.05, 0.10, 0.25, 0.50, 0.75, 0.80, 0.90, 0.95], dtype = 'float32')

with open('input_global/all_vals_conc.csv') as f:
    all_vals_conc = np.array(list(csv.reader(f,quoting=csv.QUOTE_NONNUMERIC)), dtype = 'float32')

with open('input_global/all_vals_denud.csv') as f:
    all_vals_denud = np.array(list(csv.reader(f,quoting=csv.QUOTE_NONNUMERIC)), dtype = 'float32')

glob_conc = rasterio.open('input_global/glob_con_1d_2.tif').read(1)
#somehow the nodata value for this raster is given as -128.0 which is a problem since
#the values are integers... so I have to define the mask myself.
glob_conc = np.ma.masked_array(glob_conc, mask=(glob_conc == -128))
glob_denud =  rasterio.open('input_global/glob_den_1d_2.tif').read(1, masked = True)

# Ravel the denud array to a 1-D shape
shp = glob_denud.shape
glob_denud_1d = glob_denud.reshape((shp[0]*shp[1], 1))
del glob_denud
gc.collect()
#do the same with the conc raster
shp = glob_conc.shape
glob_conc_1d = glob_conc.reshape((shp[0]*shp[1], 1))
del glob_conc
gc.collect()

glob_conc_process = glob_conc_1d[glob_conc_1d.mask == False].data

indices = np.where(glob_conc_1d.mask == False)[0]

del glob_conc_1d
gc.collect()

glob_denud_process = glob_denud_1d[glob_denud_1d.mask == False].data

del glob_denud_1d
gc.collect()

shp = glob_denud_process.shape
glob_denud_process = glob_denud_process.reshape((shp[0], 1))

shp = glob_conc_process.shape
glob_conc_process = glob_conc_process.reshape((shp[0], 1))

glob_process = np.concatenate((np.array(glob_conc_process, dtype = 'uint16'), glob_denud_process), axis = 1)

del glob_conc_process
del glob_denud_process

@jit(nopython = True)
def OCpetrodenudation_jit(input):
        conc_vals = all_vals_conc[input[0]]
        denud_vals = all_vals_denud[input[1]]
        OCpetro_vals = conc_vals*denud_vals
        q_OCpetro_vals = np.quantile(OCpetro_vals, quantile)
        return np.float32(q_OCpetro_vals)


@jit(nopython = True)
def OCpetrodenudation_par_jit(array):
    return [*map(OCpetrodenudation_jit, array)]

N = glob_process.shape[0]
P = 49 # Number of breaks (number of partitions + 1)

# Break up the indices into (roughly) equal parts
partitions = list(zip(np.linspace(0, N, P, dtype=int)[:-1],
    np.linspace(0, N, P, dtype=int)[1:]))

# Final range of indices should end +1 past last index for completeness
work = partitions[:-1]
work.append((partitions[-1][0], partitions[-1][1] + 1))

# Split the master array, base_array, into subarrays defined by the
#   starting and ending, i and j, indices

#OCpetrodenudation
with ProcessPoolExecutor(max_workers = 48) as executor:
    result = executor.map(OCpetrodenudation_par_jit, [
        glob_process[i:j,...] for i, j in work
    ])

results = list(result)    
del result
result_glob_OCpetro = np.concatenate(results, axis=0)
del results

#Record the global flux
MC_sample = np.sum(result_glob_OCpetro)*1000*1000/(10**6)/1000 #in tC y-1

output_glob_OCpetro = np.full(565251360, 0.0, dtype = 'float32')
output_glob_OCpetro[indices] = result_glob_OCpetro

del result_glob_OCpetro

output_glob_OCpetro_2d = output_glob_OCpetro.reshape((15684, 36040))

del output_glob_OCpetro

#Now evaluate the output:

with rasterio.open('glob_petrodenud.tif', 'w', **profile) as dst:
    dst.write(output_glob_OCpetro_2d, 1)
del output_glob_OCpetro_2d
    
#Mask out alluvial domain
glob_petrodenud_tif = r'glob_petrodenud.tif'
glob_petrodenud_rast = rasterio.open(glob_petrodenud_tif)
with fiona.open(r'input_global/Alluv.geojson', 'r') as source:
    Alluv = [feature["geometry"] for feature in source]
glob_petrodenud_wo_alluv = rasterio.mask.mask(glob_petrodenud_rast, Alluv, invert = True, nodata = 0.0)
del glob_petrodenud_rast
del Alluv
MC_sample_wo_alluv = np.sum(glob_petrodenud_wo_alluv[0][0])*1000*1000/(10**6)/1000 #in tC y-1

with rasterio.open('glob_petrodenud_wo_alluv.tif', 'w', **profile) as dst:
    dst.write(glob_petrodenud_wo_alluv[0][0], 1)
    
del glob_petrodenud_wo_alluv

##### Calibration #####
#######################

from xrspatial import zonal_stats
#to install this I needed to revert to python 3.9 and install numpy 1.21.0
#This does not seem to be working well in the linux environment, and so I think
#I will need to rebuild a new environment there.
import xarray as xr
import rioxarray

#Define filenames, both are rasters prepared in ArcGIS. See MC_global notebook for more notes.
basins_tif = r'input_global/Re_catch20.tif'
glob_petrodenud_wo_alluv_tif = r'glob_petrodenud_wo_alluv.tif'

# read the rasters
basins = rioxarray.open_rasterio(basins_tif)
glob_petrodenud = rioxarray.open_rasterio(glob_petrodenud_tif)

#perform the zonal statistics
zstats_df = zonal_stats(basins.sel(band=1), glob_petrodenud.sel(band=1), zone_ids = [
    1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
       18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
       35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52,
       53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
       71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83
], stats_funcs=['sum'])
del glob_petrodenud

#Add sub-basin sums to the large basins sums
#Amazon
#79: add 78
zstats_df.loc[zstats_df['zone'] == 79, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 79])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 78])[0]
#80: add new 79
zstats_df.loc[zstats_df['zone'] == 80, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 80])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 79])[0]
#82: add new 80
zstats_df.loc[zstats_df['zone'] == 82, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 82])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 80])[0]
#81: add new 82
zstats_df.loc[zstats_df['zone'] == 81, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 81])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 82])[0]
#40: add 83 and new 81
zstats_df.loc[zstats_df['zone'] == 40, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 40])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 83])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 81])[0]

#Himalayas
#44: add 69, 66, 68, 67
zstats_df.loc[zstats_df['zone'] == 44, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 44])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 69])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 66])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 68])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 67])[0]
 
#Mackenzie
#31: add 74 and 75 ## no, don't add Arctic red (74) and adjust the Re value calc in Re-compile
zstats_df.loc[zstats_df['zone'] == 31, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 31])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 75])[0]

#Taiwan
#6: add 7
zstats_df.loc[zstats_df['zone'] == 6, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 6])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 7])[0]

#Yukon:
#64: add 65 (and redo Re-compile basin calc)
zstats_df.loc[zstats_df['zone'] == 64, 'sum'] = np.array(zstats_df['sum'][zstats_df['zone'] == 64])[0] + np.array(zstats_df['sum'][zstats_df['zone'] == 65])[0]
    
#Calculate the model oxidation rate
zstats_df['OC_Petro_ox'] = zstats_df['sum']*1000*1000/(10**6)/1000 #in tC y-1

#Read the csv file of basins and their Re-estimated oxidation rates
Re_ests = pd.read_csv('input_global/basins_re_est.csv')

#Merge the data together
zstats_df = pd.merge(zstats_df, Re_ests, left_on='zone', right_on='ID')

#subset the table to exclude rows with NA values
zstats_df_nona = zstats_df.loc[-np.isnan(zstats_df['OCpetro_ox_Re']),]

#turning this off for now
#calculate the residual in tC y-1 for each basin
#zstats_df_nona['residual'] = zstats_df_nona['OC_Petro_ox'] - zstats_df_nona['OCpetro_ox_Re']

#Calculate the global residual as a ratio between 0 and 1
Global_resid = (np.sum(zstats_df_nona['OC_Petro_ox'])-np.sum(zstats_df_nona['OCpetro_ox_Re']))/np.sum(zstats_df_nona['OCpetro_ox_Re'])

#Calculate the basin specific residual, also as a ratio
#zstats_df_nona['residual_perc'] = zstats_df_nona['residual'] / zstats_df_nona['OCpetro_ox_Re']

#Next we need to save the data table and the value of the global residual in a file (CSV?)
#We'll do the global residual only for the first run, so we can define a threshold

###Second calibration with no alluv.

glob_petrodenud_wo_alluv = rioxarray.open_rasterio(glob_petrodenud_wo_alluv_tif)

zstats_df_2 = zonal_stats(basins.sel(band=1), glob_petrodenud_wo_alluv.sel(band=1), zone_ids = [
    1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
       18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
       35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52,
       53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
       71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83
], stats_funcs=['sum'])
del glob_petrodenud_wo_alluv

#Add sub-basin sums to the large basins sums
#Amazon
#79: add 78
zstats_df_2.loc[zstats_df_2['zone'] == 79, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 79])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 78])[0]
#80: add new 79
zstats_df_2.loc[zstats_df_2['zone'] == 80, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 80])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 79])[0]
#82: add new 80
zstats_df_2.loc[zstats_df_2['zone'] == 82, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 82])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 80])[0]
#81: add new 82
zstats_df_2.loc[zstats_df_2['zone'] == 81, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 81])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 82])[0]
#40: add 83 and new 81
zstats_df_2.loc[zstats_df_2['zone'] == 40, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 40])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 83])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 81])[0]

#Himalayas
#44: add 69, 66, 68, 67
zstats_df_2.loc[zstats_df_2['zone'] == 44, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 44])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 69])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 66])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 68])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 67])[0]
 
#Mackenzie
#31: add 74 and 75 ## no, don't add Arctic red (74) and adjust the Re value calc in Re-compile
zstats_df_2.loc[zstats_df_2['zone'] == 31, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 31])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 75])[0]

#Taiwan
#6: add 7
zstats_df_2.loc[zstats_df_2['zone'] == 6, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 6])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 7])[0]

#Yukon:
#64: add 65 (and redo Re-compile basin calc)
zstats_df_2.loc[zstats_df_2['zone'] == 64, 'sum'] = np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 64])[0] + np.array(zstats_df_2['sum'][zstats_df_2['zone'] == 65])[0]
    
#Calculate the model oxidation rate
zstats_df_2['OC_Petro_ox'] = zstats_df_2['sum']*1000*1000/(10**6)/1000 #in tC y-1

#Merge the data together
zstats_df_2 = pd.merge(zstats_df_2, Re_ests, left_on='zone', right_on='ID')

#subset the table to exclude rows with NA values
zstats_df_2_nona = zstats_df_2.loc[-np.isnan(zstats_df_2['OCpetro_ox_Re']),]

#turning this off for now
#calculate the residual in tC y-1 for each basin
#zstats_df_2_nona['residual'] = zstats_df_2_nona['OC_Petro_ox'] - zstats_df_2_nona['OCpetro_ox_Re']

#replacing the calculation as one where we first sum the model and Re oxidation rates and then calc the residual.
#Calculate the global residual as a ratio between 0 and 1
Global_resid_2 = (np.sum(zstats_df_2_nona['OC_Petro_ox'])-np.sum(zstats_df_2_nona['OCpetro_ox_Re']))/np.sum(zstats_df_2_nona['OCpetro_ox_Re'])

print(MC_sample)
print(Global_resid)
print(zstats_df_nona)
print(MC_sample_wo_alluv)
print(Global_resid_2)
print(zstats_df_2_nona)

end = time()
runtime = float(end - start)/60.
print('Took %.3f minutes' % runtime)
