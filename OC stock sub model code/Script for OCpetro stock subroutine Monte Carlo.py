#Import dependencies, install if not on machine yet
import gc
import pandas as pd
import os
import glob
import scipy.interpolate
import matplotlib as mpl
import csv

import rasterio as rio                                                                                                                  
from rasterio.plot import show  
import rasterio.mask
from rasterio.rio import stack
from rasterio.rio.stack import stack

#Read in all the csv files with the quantile functions
path = os.getcwd()+'\quantiles_glob'
csv_files = glob.glob(os.path.join(path, "*.csv"))

dataframes_names = {}

for f in csv_files:
      
    # read the csv file
    df = pd.read_csv(f)
    dataframes_names[os.path.basename(f).replace('.csv','')] = df

#load all the unique combinations of continent and lithology, as output by the R script
unique_conc_names = pd.read_csv(os.getcwd() + r'\all_names_conc.csv')
unique_conc_names = unique_conc_names['x'].tolist()

#Create lookup tables of the numbers in these rasters, corresponding to lithological units
#and continents. The rasters use integer numbers only as this will make the global OCpetro 
#oxidation simulations on the cluster a lot quicker!

#First we make the lookup table for lithological units made up of one lithology:
lith_list = [1, 2, 3, 6, 8, 9, 15]
cdf_liths = ['ecdf_ig', 'ecdf_m', 'ecdf_s', 'ecdf_sh', 'ecdf_ma', 'ecdf_c', 'ecdf_shb']
lith_table = pd.DataFrame(list(zip(lith_list,cdf_liths)), columns = ['unit_id', 'ecdf'])

#These lookup tables are for units that are made of several lithologies:
unit_list = [4, 5, 10, 11]
cdf_units = [['ecdf_s', 'ecdf_sh'], ['ecdf_s', 'ecdf_sh', 'ecdf_c'], ['ecdf_sh', 'ecdf_c'], ['ecdf_s', 'ecdf_c']]
ratio_list = [['ratio_sands_shales'],['ratio_sands_shales_carbs', 'ratio_shales_carbs_sands'],['ratio_shales_carbs'],
              ['ratio_sands_carbs']]
unit_table = pd.DataFrame(list(zip(unit_list,cdf_units, ratio_list)), columns = ['unit_id', 'ecdfs', 'ratios'])

#Lookup tables for continent names:
cont_id = pd.DataFrame(list(zip([1,2,3,4,5,6,7,8],['Oceania', 'South America','Antarctica', 'Australia',
                                                  'Asia', 'Africa', 'North America', 'Europe'])),
                       columns = ['cont_id', 'cont_name'])

#function to take a quantile function produced in the R script, to run
#the OCpetro stock subroutine Monte Carlo
def sample_conc(name):
    lit = int(name.split('_')[1])
    cont = int(name.split('_')[0])
    values = []
    if lit > 0 and cont > 0:
        if lit in [1, 2, 3, 6, 8, 9, 15]:
            df = dataframes_names[lith_table['ecdf'][lith_table['unit_id']==lit].tolist()[0]]
            for i in range(0,10000):
                x = np.random.uniform(low = 0, high = 1)
                x_interp = scipy.interpolate.interp1d(df['Fnx'], df['x'])
                values.append(round(x_interp(x).tolist(),2))
            return values
        
        elif lit in [4, 5, 10, 11]:
            cdfs = unit_table['ecdfs'][unit_table['unit_id']==lit].tolist()[0]
            ratio_names = unit_table['ratios'][unit_table['unit_id']==lit].tolist()[0]
            ratios = []
            continent = cont_id['cont_name'][cont_id['cont_id']==cont].tolist()[0]
            for each_ratio in ratio_names:
                    ratios.append(dataframes_names['Am_Such_rat'][each_ratio][dataframes_names['Am_Such_rat']['Continent']
                                                                      ==continent].tolist()[0])
            for i in range(0,10000):
                conc_list = []
                for each_cdf in cdfs:
                    df = dataframes_names[each_cdf]
                    x = np.random.uniform(low = 0, high = 1)
                    x_interp = scipy.interpolate.interp1d(df['Fnx'], df['x'])
                    conc_list.append(x_interp(x).tolist())
                    
                if lit in [4, 10, 11]:
                    values.append(round(ratios[0]*conc_list[0]+conc_list[1]*(1-ratios[0]),2))
                else:
                    values.append(round(ratios[0]*conc_list[0]+ratios[1]*conc_list[1]+conc_list[2]*(1-ratios[0]-ratios[1]),2))
            return values
    
        elif lit in [7, 12, 13, 14]:
            for i in range(0,10000):
                values.append(0.00)
            return values
    else:
        for i in range(0,10000):
            values.append(0.00)
        return values

#Array which will contain the Monte Carlo simulation results for 
#each unique value of denudation and lithology
all_vals_conc = []
all_names_conc = []

#Loop through all unique combinations of lithology and continent, running the
#above function
#Collecting the results in the arrays defined above
for combination in unique_conc_names:
    all_vals_conc.append(sample_conc(combination))
    all_names_conc.append(combination)

#Write the array of the subroutine Monte Carlo simulations to a csv, to be read
#in each successive OCpetro oxidation simulation on the computer cluster
with open('all_vals_conc.csv', 'w', newline='') as f:
    write = csv.writer(f)
    write.writerows(all_vals_conc)
    
#We also need a global raster where each cell value corresponds to the right
#index in the Monte Carlo simulation output file. This is important for 
#computational efficiency:

#first we open the same grids as we did in the R script,
#of 1 km resolution lithology and continent:
glob_lith = rasterio.open('glob_lit_lor2.tif').read(1)
glob_cont = rasterio.open('glob_con_lors.tif').read(1)
profile = glob_lith.profile

# Ravel the lith array to a 1-D shape
shp = glob_lith.shape
arr_flat = glob_lith.reshape((shp[0]*shp[1], 1))
base_array_glob_conc = arr_flat
#do the same with the cont map and stack the arrays
shp = glob_cont.shape
arr_flat = glob_cont.reshape((shp[0]*shp[1], 1))
base_array_glob_conc = np.concatenate((base_array_glob_conc, arr_flat), axis = 1)

del arr_flat
del glob_lith
del glob_cont
gc.collect()

#function to derive the index value for the subroutine Monte Carlo outputs 
def conc_array(input):
   lit = input[0]
   cont = input[1]
   if lit > 0 and cont > 0:
       name = (str(int(cont))+'_'+str(int(lit)))
       value = np.where(all_names_conc == name)[0][0]
       return value
   else:
       return -32768.0

#Map the above function over the spatial raster files    
conc_array_1d = [*map(conc_array, base_array_glob_conc)]
output_conc_array_1d = np.array(conc_array_1d).reshape((15684, 36040))

#output the spatial raster file which will form an input for the 
#OCpetro oxidation model simulations.
with rasterio.open('glob_con_1d.tif', 'w', **profile) as dst:
    dst.write(output_conc_array_1d, 1)
    
#This output tif file was subsequently modified to an 8bit tif, with
#pixel type = char, using the Copy Raster tool in ArcGIS Pro, changing 
#the nodata value to -128. 
#the subsequent output file is glob_con_1d_2.tif