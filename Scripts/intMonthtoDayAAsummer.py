"""  runForecasts.py
	Script to run the various forecast scripts
	Run with e.g. python runForecasts.py 
	Author:
		Alek Petty

	Update history:
		04/20/2018: Version 1
"""


import numpy as np
from pylab import *
# Apr observed, May, June, July Aug Sep predicted

forecastDays=[-15, 15, 46, 76, 107, 138]
# 
forecasts=[4.16, 7.52, 10.43, 12.72, 14.16, 14.44]
days=np.arange(123)

# jan, feb march predicted
#forecastDays=[15, 45, 75]
#forecasts=[2.07]
#days=np.arange(90)


dailyForecasts= np.interp(days, forecastDays, forecasts)


#Quadratic fit
fit = np.polyfit(forecastDays, forecasts, 2)
fitPoly=np.poly1d(fit)
quadfitForecasts=fitPoly(days)

savetxt('./NASAGSFC_001_total-area2022-sh.txt', [quadfitForecasts], fmt='%.4f', delimiter=',')
plot(forecastDays, forecasts, marker='x')
plot(days, quadfitForecasts)
show()