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
import sys

sys.path.append('./forecasts/')
sys.path.append('./gridding/')
sys.path.append('./plotting/')
import calcForecastYears
import grid_iceconcA
import grid_iceconcAA
import plotForecasts

pmonth=12 #9= September
month=11 # here month starts from 1
monthStr='%02d' %(month)
year=2022
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
subprocess.check_call(['./getData/wgetDataNRT.sh', monthStr, yearStr, sensor, hemisphere])

#grid_iceconc_cdr.main(year, month-1, poleStr=poleStr)

#if (hemisphere=='s'):
	grid_iceconcAA.main(year, month-1)
#else:
#	grid_iceconcA.main(year, month-1)

# Run forecasts
calcForecastYears.main(year, month, pmonth, hemStr=hemStr, iceType=iceType, region=region)

# Plot forecast
plotForecasts.main(year, month, pmonth, iceType=iceType, hemStr=hemStr, minval=9, maxval=15, region=region)




