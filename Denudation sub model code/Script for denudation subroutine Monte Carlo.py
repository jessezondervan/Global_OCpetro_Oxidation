
#Import dependencies, install if not on machine yet
import pandas as pd
import os
import glob
import scipy.interpolate
import matplotlib as mpl
import csv

#Read in all the csv files with the quantile functions
path = os.getcwd()+'\quantiles_glob'
csv_files = glob.glob(os.path.join(path, "*.csv"))

dataframes_names = {}

for f in csv_files:
      
    # read the csv file
    df = pd.read_csv(f)
    dataframes_names[os.path.basename(f).replace('.csv','')] = df

#function to take a quantile function produced in the R script, to run
#the denudation subroutine Monte Carlo
def sample_denud(name):
    df = dataframes_names[name]
    taus = np.arange(0, 1.01, 0.01).tolist()
    quants = df['x'].tolist()
    values = []
    for i in range(0,10000):
        x = np.random.uniform(low = 0, high = 1)
        x_interp = scipy.interpolate.interp1d(taus, quants)
        values.append(x_interp(x).tolist())
    return values
 
#Array which will contain the Monte Carlo simulation results for 
#each unique value of denudation and lithology 
all_vals = []

#Looping for all unique combinations and running the sample_denud function
#Collecting the results in the array defined above
for key in dataframes_names:
    if key == 'Am_Such_rat':
        continue
    elif 'ecdf' in key:
        continue
    else:
        all_vals.append(sample_denud(key))

#Write the array of the denudation model Monte Carlo results to a csv, to be read
#in each successive OCpetro oxidation simulation on the computer cluster
with open('all_vals_denud.csv', 'w') as f:
    write = csv.writer(f)
    write.writerow(all_vals_denud)
    
#We also need a global raster where each cell value corresponds to the right
#index in the Monte Carlo simulation output file. This is important for 
#computational efficiency:

#first we open the same grids as we did in the R script,
#of 1 km resolution lithology and denudation:
glob_lith = rasterio.open('glob_lit_lor2.tif').read(1)
glob_denud =  rasterio.open('glob_den_lor3.tif').read(1)
profile = glob_denud.profile

# Ravel the denud array to a 1-D shape
shp = glob_denud.shape
arr_flat = glob_denud.reshape((shp[0]*shp[1], 1))
base_array_glob_denud = arr_flat
#do the same with the lith map and stack the arrays
shp = glob_lith.shape
arr_flat = glob_lith.reshape((shp[0]*shp[1], 1))
base_array_glob_denud = np.concatenate((base_array_glob_denud, arr_flat), axis = 1)

del arr_flat
del glob_denud
gc.collect()

#function to derive the index value for the subroutine Monte Carlo outputs 
def denud_array(input):
    den = input[0]
    lit = input[1]
    if den > 0. and lit > 0.:
        name = (str(int(den))+'_'+str(int(lit)))
        value = np.where(all_names_denud == name)[0][0]
        return value
    else:
        return -32768.0

#Map the above function over the spatial raster files
denud_array_1d = [*map(denud_array, base_array_glob_denud)]
output_denud_array_1d = np.array(denud_array_1d).reshape((15684, 36040))

#output the spatial raster file which will form an input for the 
#OCpetro oxidation model simulations.
with rasterio.open('glob_den_1d.tif', 'w', **profile) as dst:
    dst.write(output_denud_array_1d, 1)
    
#This output tif file was subsequently modified to an 16bit tif, with
#pixel type = unsigned short, using the Copy Raster tool in ArcGIS Pro, changing 
#the nodata value to 65535. 
#the subsequent output file is glob_den_1d_2.tif