"""  runForecastsAA.py
	Script to run the various forecast scripts for generating SIPN South forecasts
	Run with e.g. python runForecastsAA.py 
	Author:
		Alek Petty

	Update history:
		01/30/2022: Version 1
"""


#


import subprocess
import shlex
import sys

sys.path.append('./forecasts/')
sys.path.append('./gridding/')
sys.path.append('./plotting/')
import calcForecastYears
import grid_iceconcA
import grid_iceconcAA
import plotForecasts


# STEP 1: MAKE SURE WE HAVE UPDATED ICE EXTENT DATA FROM THE NSIDC.
# - Manually replace what's in the Data/SeaIceIndex/ folder with latest INDEX data.
# - See here: https://masie_web.apps.nsidc.org/pub//DATASETS/NOAA/G02135/south/monthly/data/

# STEP 2: Get updated NRT SIC DATA AND GRID
# - see automated functions below

# STEP 3: run forecasts of upcoming December, January, February and March
# - see automated functions below. 

# STEP 4: convert the monthly forecasts to daily forecasts (requested by SIPN South)
# - see the intMonthstoDayAA.py script and m,anually input values. 


pmonth=3 #9= September
month=11 # here month starts from 1
monthStr='%02d' %(month)
year=2023
yearStr=str(year)
sensor='f18'
hemisphere='s'


if (hemisphere=='s'):
	hemStr='S'
	poleStr='AA'
else:
	hemStr='N'
	poleStr='A'

region=0
iceType='area'

# ADD A CHECK TO SEE IF DATA EXISTS FIRST!
#subprocess.check_call(['./getData/get_nrt_cdr_monthly.sh', monthStr, yearStr, sensor, hemisphere])

#grid_iceconc_cdr.main(year, month-1, poleStr=poleStr)

#if (hemisphere=='s'):
#	grid_iceconcAA.main(year, month-1)
#else:
#	grid_iceconcA.main(year, month-1)

# Run forecasts
calcForecastYears.main(year, month, pmonth, hemStr=hemStr, iceType=iceType, region=region)

# Plot forecast
#plotForecasts.main(year, month, pmonth, iceType=iceType, hemStr=hemStr, minval=9, maxval=15, region=region)




