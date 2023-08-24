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
# nov observed, dec, jan, feb march predicted
forecastDays=[-15, 15, 45, 75, 105]
# late 2019 forcastsforecasts=[11.11, 6.84, 3.33, 2.07, 2.73]
forecasts=[10.71, 6.39, 2.75, 1.67, 2.29]
days=np.arange(90)

# jan, feb march predicted
#forecastDays=[15, 45, 75]
#forecasts=[2.07]
#days=np.arange(90)


dailyForecasts= np.interp(days, forecastDays, forecasts)


#Quadratic fit
fit = np.polyfit(forecastDays, forecasts, 2)
fitPoly=np.poly1d(fit)
quadfitForecasts=fitPoly(days)

savetxt('./NASAGSFC_001_total-area2023.txt', [quadfitForecasts], fmt='%.4f', delimiter=',')
plot(forecastDays, forecasts, marker='x')
plot(days, quadfitForecasts)
show()