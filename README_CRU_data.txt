###############################################################################
#
# This readme explains how to use the CRU_data_extraction_v2.R function.
#
# As with any function for obtaining data make sure the data makes sense prior
# to using it in case anything went wrong extracting the data.
#
###############################################################################

The CRU data extraction function requires several inputs. These include the
location (latitude and longitude), a starting year 1901 or greater an ending 
year no later than 2016 (unless the data files are updated) a vector of 
site id's if not in the location file and a name for output files.

Location file:
The easiest way to use the function will be to have a file with the location
and site data. This file should be a csv with a minimum of latitude, longitude
and SiteNum fields. Latitude can be either latitude or lat, and longitude can
be longitude, long, or lon. The site number must be SiteNum or the siteid
variable to the function must be provided a vector of site numbers equal
in length to the number of locations provided. If SiteNum is included in the 
locations file then the siteid variable does not need to be specified or 
can be set to the default value of 1.

Other function inputs:
crufile - the CRU data file to extract data from
ystart - the starting year
yend - the ending year (must be greater than ystart)]
varout - a string containing the variable name being extracted
siteid - a vector of site ids or set to the value 1 if the site ids are in
	the locations file (see above).
infile - the file with the location data
outfile - the name of the file for the output.

EXAMPLES:
An example of how to use this function to extract data for all 10 CRU
variables can be seen by opening GetData_v2.R. The script loads the function
finds the cru data files has a vector of variables and automatically extracts
data for each of the variables for the locations in MortalitySites.csv.

The function used for extraction is in the file "CRU_data_extraction_v2.R" if
there are unanticipated errors with the function or modification is needed. If
the function is modified it should be saved with a new version number PRIOR TO
MAKING CHANGES.

CRU variable info:
# CRU climate variable notes:
#	cld - cloud cover (percentage %)
#	dtr - diurnal temperature range (degrees Celsius)
#	frs - frost day frequency (days)
#	pet - potential evapotranspiration (millimeters per day)
#	pre - precipitation (millimeters per month).
#	rhm - relative humidity (percentage %)
#	ssh - sunshine duration (hours)
#	tmp - daily mean temperature (degrees Celsius)
#	tmn - monthly average daily minimum temperature (degrees Celsius)
#	tmx - monthly average daily maximum temperature (degrees Celsius)
#	vap - vapor pressure (hectopascals hPa)
#	wet - wed day frequency (days)
#	wnd - wind speed (metres per second m/s)
#
#	source: https://crudata.uea.ac.uk/cru/data/hrg/#info

# Mike email about files:

Each file contains one CRU variable. 
I organized the data in a structure where the sites are stacked. 
Each site has a unique site number taken from “MortalitySites.csv,” which I think you sent me.  
I thought this would work well if you need an index in BUGS for each site.

# UPDATED email from Nov 28, 2017:
I updated the files I created last week.  
I went back and did additional checks and found an error in the code.  
This has been corrected and I did the checks I didn’t get to last week.  
I spot checked output against a netcdf viewer to ensure its values at locations matched my output.  
I also plotted some of the tmp values for each month at some of the sites to make sure they made sense. 
I checked a few sites from the northern and a few from the southern hemisphere 
(there is a map with points for sites in the folder I sent last week “location_points.png”).
 
Along with checking the output I fully commented the file with the extraction function. 
I also generalized the function so it should be easy to use to extract additional data in the future.
 
NEW FILES:
wet_2017_11_28.csv
vap_2017_11_28.csv
tmx_2017_11_28.csv
tmp_2017_11_28.csv
tmn_2017_11_28.csv
pre_2017_11_28.csv
pet_2017_11_28.csv
frs_2017_11_28.csv
dtr_2017_11_28.csv
cld_2017_11_28.csv
 
CRU DATA:
I added some new data files.  There are now data for 10 variables. 
All of the below descriptions are also in the R script for the extraction function.
 
# CRU climate variable notes:
#    cld - cloud cover (percentage %)
#    dtr - diurnal temperature range (degrees Celsius)
#    frs - frost day frequency (days)
#    pet - potential evapotranspiration (millimeters per day)
#    pre - precipitation (millimeters per month).
#    rhm - relative humidity (percentage %)
#    ssh - sunshine duration (hours)
#    tmp - daily mean temperature (degrees Celsius)
#    tmn - monthly average daily minimum temperature (degrees Celsius)
#    tmx - monthly average daily maximum temperature (degrees Celsius)
#    vap - vapor pressure (hectopascals hPa)
#    wet - wed day frequency (days)
#    wnd - wind speed (metres per second m/s)