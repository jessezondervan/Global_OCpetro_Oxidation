#The purpose of this script is to find all the unique combinations of lithology 
#and continent, so we can generate CSVs of quantiles against
#OCpetro stock values (quantile function). These CSVs can then be read into a 
#python script where called upon, and used to interpolate an OCpetro stock value 
#based on a randomly generated quantile.

setwd("C:/Users/zonde")

##Required Libraries (use install.packages to install these)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(GISTools)
library(earth)
library(reshape2)
library(ggplot2) 
library(gridExtra)
library(quantreg)

####### read in the USGS Rock data files. This is data from the National
#Geochemical Database: Rock. Accessible from https://mrdata.usgs.gov/ngdb/rock/ 
#[last accessed 2 Feb 2023]
setwd('C:/Users/zonde/OneDrive - Durham University/ROC-CO2-Jesse/USGS/tab_delimited')

Geospatial <-  read.csv("tblRockGeoData.csv", header=TRUE, na.string = "NULL")
MajorChem <-  read.csv("xtbMajorChem.csv", header=TRUE, na.string = "NULL")

####### merge the chemistry with the geospatial data
MajorChemGs <- merge(Geospatial, MajorChem, by = "LAB_ID")

#read in the lookup table. This table was generated as follows:
#USGS_specnames <- table(MajorChemGs_OM$SPEC_NAME) # this line outputs all the 
#specimen names as reported by workers and how often these occur
#write.csv(USGS_specnames, "USGS_specnames.csv", row.names=FALSE)
#I then opened this csv in Excel and line by line put the relevant lithological 
#class code in a new column, following classification in table S5.
class_codes <-  read.csv("USGS_specnames.csv", header=TRUE, na.string = "NULL")

#append a new column in the USGS dataframe with the class code for each sample
MajorChemGs$class_code <- lapply(MajorChemGs$SPEC_NAME, function(x) 
  class_codes$class.code[match(x, class_codes$SPEC_NAME)])
MajorChemGs$class_code <- unlist(MajorChemGs$class_code)

#THe organic C is reported in weight percent, which we convert to units of
#kg/m3 (OCpetro stock) using density. This is required as OCpetro denudation
#rate can then be calculated by multiplying the OCpetro stock with denudation
#(in m/yr).
#First append a new column in the USGS dataframe with the representative 
#density for each sample, see Table S5.
units <- c('s', 'c', 'sh', 'shb', 'e', 'ig', 'm', 'ma')
density <- c(2200, 2600, 2350, 2350, 2200, 2850, 2750, 2800)
density_lkp <- data.frame(units, density)

MajorChemGs$density <- lapply(MajorChemGs$class_code, function(x) 
  density_lkp$density[match(x, density_lkp$units)])
MajorChemGs$density <- unlist(MajorChemGs$density)

MajorChemGs$OM_conc <- MajorChemGs$C_ORG_DFF/100*MajorChemGs$density

setwd('C:/Users/zonde/OneDrive - Nexus365/ROC-CO2-Jesse')

#OC stock quantiles:

for(i in c('s', 'c', 'sh', 'shb', 'ig', 'm', 'ma'))
{
  ecdf_x = ecdf((MajorChemGs$OM_conc[which(MajorChemGs$class_code == i & !is.na(MajorChemGs$OM_conc))]))
  knots_x = knots(ecdf_x)
  quants_x = ecdf_x(knots(ecdf_x))
  cdf_x = do.call(rbind, Map(data.frame, x=knots_x, Fnx=quants_x))
  cdf_x= rbind(data.frame(x = min(knots_x),  Fnx = 0), cdf_x)
  name = paste0("ecdf_", i)
  write.csv(cdf_x, paste0('~/OneDrive - University College London/ROC-CO2-Jesse/Monte Carlo/quantiles_glob/', name, '.csv'))
}

#As there are lithological classes with mixes of the above lithologies (see Table S5),
#we now also need a table which lists the ratios of these mixes, which are
#calibrated by continent using the dataset by Amiotte Suchet et al. (2013):
#https://doi.org/10.1029/2002GB001891

#let's read in the coverage percentage of lithology types from a simplified version
#of table 3 in Amiotte Suchet et al. (2013):
Amiotte_Suchet <- read.csv('C:/Users/zonde/OneDrive - University College London/ROC-CO2-Jesse/Amiotte-Suchet simple.csv', header = TRUE)
colnames(Amiotte_Suchet)

#calculate the ratio of sandstones to carbonates in each continent
Amiotte_Suchet$ratio_sands_carbs <- Amiotte_Suchet[,2]/(Amiotte_Suchet[,2]+Amiotte_Suchet[,4])
#and make sure that, since in the GLiM classification none of the lithologies dominates
#(see Hartmann & Moosdorf, 2012), we will have the minority lithology at at least 25%:
Amiotte_Suchet$ratio_sands_carbs <-pmax(Amiotte_Suchet$ratio_sands_carbs, 0.25)
Amiotte_Suchet$ratio_sands_carbs <-pmin(Amiotte_Suchet$ratio_sands_carbs, 0.75)

Amiotte_Suchet$ratio_shales_carbs <- Amiotte_Suchet[,3]/(Amiotte_Suchet[,3]+Amiotte_Suchet[,4])
Amiotte_Suchet$ratio_shales_carbs <-pmax(Amiotte_Suchet$ratio_shales_carbs, 0.25)
Amiotte_Suchet$ratio_shales_carbs <-pmin(Amiotte_Suchet$ratio_shales_carbs, 0.75)

Amiotte_Suchet$ratio_sands_shales <- Amiotte_Suchet[,2]/(Amiotte_Suchet[,2]+Amiotte_Suchet[,3])
Amiotte_Suchet$ratio_sands_shales <-pmax(Amiotte_Suchet$ratio_sands_shales, 0.25)
Amiotte_Suchet$ratio_sands_shales <-pmin(Amiotte_Suchet$ratio_sands_shales, 0.75)

Amiotte_Suchet$ratio_sands_shales_carbs <- Amiotte_Suchet[,2]/(Amiotte_Suchet[,2]+Amiotte_Suchet[,3]+Amiotte_Suchet[,4])
Amiotte_Suchet$ratio_shales_carbs_sands <- Amiotte_Suchet[,3]/(Amiotte_Suchet[,2]+Amiotte_Suchet[,3]+Amiotte_Suchet[,4])
Amiotte_Suchet$ratio_carbs_sands_shales <- Amiotte_Suchet[,4]/(Amiotte_Suchet[,2]+Amiotte_Suchet[,3]+Amiotte_Suchet[,4])


Am_Such_conc = Amiotte_Suchet[,-c(2:7)]
write.csv(Am_Such_conc, '~/OneDrive - University College London/ROC-CO2-Jesse/Monte Carlo/quantiles_glob/Am_Such_rat.csv')
#I added Oceania in by hand in Excel, as a copy of the Australia row.

#Now we have csv files containing a quantile function for each combination of lithology
#and continent, and a file with ratios of relevant lithologies in some mapped units.

#These files will be used for a subroutine Monte Carlo simulation in a python script

#To run the Monte Carlo subroutine we also need a list with all the unique combinations
#of continent and lithology.

#first load the rasters: 
#These are 1-km grids with values of lithology and continent
#Lith is a rasterized version of the GLiM dataset, with numbers 
#corresponding to the lithology classes that we defined (See Table S5)
#1 equals to igneous rocks, 2 equals crystalline metamorphic rocks, and all other numbers
#equal to sedimentary rock types.
#Denud was generated by using the median regression line through basins of each rock type
#(as see below in lines 99-105), the function of which was applied to the 
#Geomorpho90m slope dataset (Amatulli et al. 2020) using the ArcGIS Pro raster calculator.
#Geomorpho90m dataset: https://doi.org/10.1038/s41597-020-0479-6
#Make sure both raster are 1 km resolution, Mollweide WGS84 (equal area) projection,
#with the lithology and continent grids snapped to the denudation grid (see subroutine
#for denudation rates).
continent = raster("~/OneDrive - University College London/ROC-CO2-Jesse/Monte Carlo/glob_con_lors.tif")
lith = raster("~/OneDrive - University College London/ROC-CO2-Jesse/Monte Carlo/glob_lit_lor2.tif")

#Create a raster stack
Conc = stack(continent, lith)

#find the unique combination of values in the raster stacks
#We start with the erosion rates
unique_Conc <- unique(Conc)

Getname <- function(cont,litho){
  name <- paste0(as.integer(cont),'_', as.integer(litho))
}

all_names_con <- mapply(FUN = Getname, unique_Conc[complete.cases(unique_Conc),][,1], unique_Conc[complete.cases(unique_Conc),][,2])
write.csv(all_names_con,'~/OneDrive - University College London/ROC-CO2-Jesse/Monte Carlo/all_names_conc.csv')

