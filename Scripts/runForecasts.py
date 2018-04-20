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

month=2 # here month starts from 1
year=2017
sensor=f17

# DO A CHECK TO SEE IF DATA EXISTS FIRST!


subprocess.call('./getData/wgetData.sh month year sensor %s %s %s' % (month, year, sensor) )

# DO A CHECK TO SEE IF DATA EXISTS FIRST!
import sys
sys.path.append('/gridding/')
import grid_iceconcA

grid_iceconcA.main(year, month-1)

# Run forecasts
import calcForecastYears

calcForecastYears.main(year, 6, 9))


# Plot forecast
import plotForecasts

plotForecasts.main(year, 6, 9))




