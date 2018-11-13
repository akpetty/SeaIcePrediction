"""  runForecasts.py
	Script to run the various forecast scripts
	Run with e.g. python runForecasts.py 
	Author:
		Alek Petty

	Update history:
		04/20/2018: Version 1
"""

import subprocess
import shlex

pmonth=9 #9= September
month=7 # here month starts from 1
monthStr='%02d' %(month)
year=2018
yearStr=str(year)
sensor='f18'
hemisphere='s'
region=0

# ADD A CHECK TO SEE IF DATA EXISTS FIRST!

#subprocess.check_call(['./getData/wgetDataNRT.sh', monthStr, yearStr, sensor, hemisphere])


# ADD A CHECK TO SEE IF DATA EXISTS FIRST!


import sys
sys.path.append('./gridding/')
import grid_iceconcAA

#grid_iceconcAA.main(year, month-1)

# Run forecasts
sys.path.append('./forecasts/')
import calcForecastYears

iceType='area'
calcForecastYears.main(year, month, pmonth, hemStr='N', iceType=iceType, region=region)

# Plot forecast
sys.path.append('./plotting/')
import plotForecasts

plotForecasts.main(year, month, pmonth, iceType=iceType, hemStr='N', minval=1, maxval=7, region=region)




