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
month=5 # here month starts from 1
monthStr='%02d' %(month)
year=2018
yearStr=str(year)
sensor='f18'
hemisphere='n'



# ADD A CHECK TO SEE IF DATA EXISTS FIRST!

subprocess.check_call(['./getData/wgetDataNRT.sh', monthStr, yearStr, sensor, hemisphere])


# ADD A CHECK TO SEE IF DATA EXISTS FIRST!


import sys
sys.path.append('./gridding/')
import grid_iceconcA

grid_iceconcA.main(year, month-1)

# Run forecasts
sys.path.append('./forecasts/')
import calcForecastYears

calcForecastYears.main(year, month, pmonth)

# Plot forecast
sys.path.append('./plotting/')
import plotForecasts

plotForecasts.main(year, month, pmonth)




