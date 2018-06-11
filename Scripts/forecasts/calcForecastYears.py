"""  plotForecasts.py
	Script to generate sea ice forecasts using a given variable (or multiple variables) for a given predictor and predictand month
	Run with e.g. python calcForecastYears.py 
	Author:
		Alek Petty

	Update history:
		04/20/2018: Version 1
"""

import matplotlib
matplotlib.use("AGG")

import sys
sys.path.append('../')
import forecast_funcs as ff
from pylab import *


def main(year, fmonth, pmonth, fvars=['conc'], iceType='extent', hemStr='N', siiVersion='v3.0', startYear=1980, weight=1, region=0, numYearsReq=5, plotForecast=1, plotSkill=1, outSkill=1, outLine=1, outWeights=1):
	""" 
	Main sea ice forecast script. Can be run here (and looped )

	Args:
		startYear=1980 # Start of forecast training
		endYear=2018 # End of forecast
		fmonth=11 #6=June, 9=Sep #  Forecast month
		pmonth=2 #9=SEP # Predicted month

		fvars: forecast variables (e.g. 'conc'). Need to be in brackets.
		weight: Spatially weighting the data (1=True)
		iceType: Ice type being forecast ('extent' or 'area')
		hemStr: Hemipshere (N or S)
		siiVersion: version of the NSIDC sea ice index (if used as the observed ice state)
		startYear: Start year of training data (default = 1980)
		region: Region being forecast (0 is pan Arctic/Antartctic)
			#0 implies pan-Arctic or Antarctic
			#2 Weddell Sea
			#3 Indian Ocean
			#4 Pacific Ocean
			#5 Ross Sea
			#6 Amundsen/BHausen Sea
			#A Alaskan

		numYearsReq: (defaul 5) Number of years required in a grid cell for it to count towards the training data
		plotSkill: =1 for plotting the skill
		outSkill: =1 far saving the skill values
		outLine: =1 for saving the forecast time series
		outWeights: =1 for saving the weightings

	Returns:
		A dumped Python array including the folowing variables:
			Observed ice extent
			Observed detrended ice extent from linear trend persistence (LTP)
			LTP extent
			Absolute forecast of sea ice extent (added LTP)
			Detrended (forecast) ice extent from LTP
			Forecast anomaly from observed
			Estimated forecast error (1 SD)

	"""
	
	rawDataPath = '/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Data/' 
	derivedDataPath = '/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/DataOutput/'
	

	if (hemStr=='S'):
		saveDataPath=derivedDataPath+'/Antarctic/'
		figPath='../../Figures/'+'/Antarctic/YearlyPredictions/'
	elif (hemStr=='N'):
		saveDataPath=derivedDataPath+'/Arctic/'
		figPath='../../Figures/'+'/Arctic/YearlyPredictions/'
	
	if (region=='A'):
		saveDataPath=derivedDataPath+'/Alaska/'
		figPath='../../Figures/'+'/Alaska/YearlyPredictions/'

	print ('Forecast year:', year)
	print ('Forecast data month:', fmonth, 'Predicted month:', pmonth)
	print ('Variables:', fvars)
	print ('Hemisphere:', hemStr)
	print ('Ice type predicted', iceType)
	print ('Weighted:', weight)

	varStrsOut=''.join(fvars)
	outStr='forecastDump'+iceType+varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(year)+'W'+str(weight)+'SII'+siiVersion

	if (year<2018):
		anomObsT=1
	else:
		anomObsT=0

	
	forecastVals=ff.CalcForecastMultiVar(rawDataPath, derivedDataPath, year, startYear, fvars, fmonth, pmonth=pmonth,
			region=region, anomObs=anomObsT, numYearsReq=numYearsReq, weight=weight, outWeights=outWeights, 
			icetype=iceType, hemStr=hemStr, siiVersion=siiVersion)
	
	print ('Observed ice state:', forecastVals[0])
	print ('Linear trend presistence:', forecastVals[2])
	print ('Forecast ice state:', forecastVals[3])
	print ('Detrended observed ice state :', forecastVals[1])
	print ('Detrended observed ice state :', forecastVals[4])
	print ('Forecast anomaly :', forecastVals[5])
	
	array(forecastVals).dump(saveDataPath+outStr)

	if (plotForecast==1):
		if (region==0):
			years, extent = ff.get_ice_extentN(rawDataPath, pmonth, startYear, year, 
					icetype=iceType, version=siiVersion, hemStr=hemStr)
			
			ff.plotForecastOneYear(figPath, years, extent, year, forecastVals, outStr, iceType, minval=15, maxval=20)


		elif (region=='A'):
			poleStr='A'

			extent=loadtxt(derivedDataPath+'/Extent/'+'ice_'+iceType+'_M'+str(pmonth)+'R'+str(region)+'_'+str(startYear)+'2017'+poleStr)
			extent=extent[0:year-startYear+1]
			years=np.arange(startYear, startYear+size(extent), 1)

			ff.plotForecastOneYear(figPath, years, extent, year, forecastVals, outStr, iceType, minval=0, maxval=2)

		

#-- run main program
if __name__ == '__main__':
	#main(2015, 6, 9)
	for y in range(1990, 2017+1, 1):
		main(y, 7, 9, hemStr='N', startYear=1979, region=0)
	#	for m in range(startMonth, endMonth+1):
	#		print (y, m)
	#		
				













