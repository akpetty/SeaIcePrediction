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
import os

def main(year, fmonth, pmonth, fvars=['conc'], iceType='extent', hemStr='N', siiVersion='v3.0', anomObsT=1, startYear=1979, weight=1, region=0, numYearsReq=5, plotForecast=1, plotSkill=1, outSkill=1, outLine=1, outWeights=1):
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
		anomObsT: flag for if we have data for the given forecast month (1 = True)
		startYear: Start year of training data (default = 1980)
		region: Region being forecast (0 is pan Arctic/Antartctic)
			#0 implies pan-Arctic or Antarctic
			#2 Weddell Sea
			#3 Indian Ocean
			#4 Pacific Ocean
			#5 Ross Sea
			#6 Amundsen/BHausen Sea
			#A Alaskan
			#S Siberian
			#At Atlantic
			#C Canadian

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
	repoPath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/'
	rawDataPath = repoPath+'/Data/' 
	derivedDataPath = repoPath+'/DataOutput/'
	

	if (hemStr=='S'):
		saveDataPath=derivedDataPath+'forecasts/Antarctic/'
		figPath=repoPath+'Figures/forecasts/Antarctic/YearlyPredictions/'
	elif (hemStr=='N'):
		saveDataPath=derivedDataPath+'forecasts/Arctic/'
		figPath=repoPath+'Figures/forecasts/Arctic/YearlyPredictions/'
	
	if (region=='A'):
		saveDataPath=derivedDataPath+'forecasts/Alaska/'
		figPath=repoPath+'/Figures/Alaska/YearlyPredictions/'

	if (region=='Alaska_v2'):
		saveDataPath=derivedDataPath+'forecasts/Alaska/'
		figPath=repoPath+'/Figures/Alaska/YearlyPredictions/'

	if (region=='Siberia_v2'):
		saveDataPath=derivedDataPath+'forecasts/Siberia/'
		figPath=repoPath+'/Figures/Siberia/YearlyPredictions/'
	
	if (region=='Atlantic_v2'):
		saveDataPath=derivedDataPath+'forecasts/Atlantic/'
		figPath=repoPath+'/Figures/Atlantic/YearlyPredictions/'

	if (region=='Canadian_v2'):
		saveDataPath=derivedDataPath+'forecasts/Canadian/'
		figPath=repoPath+'/Figures/Canadian/YearlyPredictions/'

	if (region=='Alaska_v3'):
		saveDataPath=derivedDataPath+'forecasts/Alaska/'
		figPath=repoPath+'/Figures/Alaska/YearlyPredictions/'

	if (region=='Siberia_v3'):
		saveDataPath=derivedDataPath+'forecasts/Siberia/'
		figPath=repoPath+'/Figures/Siberia/YearlyPredictions/'
	
	if (region=='Atlantic_v3'):
		saveDataPath=derivedDataPath+'forecasts/Atlantic/'
		figPath=repoPath+'/Figures/Atlantic/YearlyPredictions/'

	if (region=='Canadian_v3'):
		saveDataPath=derivedDataPath+'forecasts/Canadian/'
		figPath=repoPath+'/Figures/Canadian/YearlyPredictions/'

	print('save path:', saveDataPath)
	print('figure path:', figPath)

	if not os.path.exists(saveDataPath):
		os.makedirs(saveDataPath)

	if not os.path.exists(figPath):
		os.makedirs(figPath)


	print ('Forecast year:', year)
	print ('Forecast data month:', fmonth, 'Predicted month:', pmonth)
	print ('Variables:', fvars)
	print ('Hemisphere:', hemStr)
	print ('Ice type predicted', iceType)
	print ('Weighted:', weight)

	varStrsOut=''.join(fvars)
	outStr='forecastDump'+iceType+varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(year)+'W'+str(weight)+'SII'+siiVersion+''


	print('Do we have observational evidence of this year?', anomObsT)
	
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
			#print(years)
			#print(extent)
			
			ff.plotForecastOneYear(figPath, years, extent, year, forecastVals, outStr, iceType, minval=np.floor(np.percentile(extent, 1)), maxval=np.ceil(np.percentile(extent, 99)))


		elif isinstance(region, str):

			try:
				extent=loadtxt(derivedDataPath+'/Extent/'+'ice_'+iceType+'_M'+str(pmonth)+'_19792021'+region, delimiter=',')
			except:
				extent=loadtxt(derivedDataPath+'/Extent/'+'ice_'+iceType+'_M'+str(pmonth)+'_19792021'+region, delimiter=',') 
		
			extent=extent[0:year-startYear+1]
			years=np.arange(startYear, startYear+size(extent), 1)

			print('Plotting and saving to: ', figPath)
			ff.plotForecastOneYear(figPath, years, extent, year, forecastVals, outStr, iceType, minval=np.floor(np.percentile(extent, 1)), maxval=np.ceil(np.percentile(extent, 99)))

	return year, forecastVals[3], extent[-1]

#-- run main program
if __name__ == '__main__':
	#main(2015, 6, 9)
	years=[]
	forecasts=[]
	extents=[]

	region='Canadian_v3'
	if isinstance(region, str):
		region_str = region
	else:
		region_str=str(region)

	start_year=1990
	end_year=2021
	m=8
	p=9
	iceType='extent'

	for y in range(start_year, end_year+1, 1):
		year, forecast, extent = main(y, m, p, iceType=iceType, hemStr='N', startYear=1979, region=region)
	#	for m in range(startMonth, endMonth+1):
		years.append(year)
		forecasts.append(forecast)
		extents.append(extent)
	print(years)
	print(forecasts)
	print(extents)

	savetxt('../../DataOutput/forecasts/GSFC_Petty_time'+str(m)+str(p)+str(start_year)+'-'+str(end_year)+region_str+'.txt', np.asarray(years).reshape(1, np.asarray(years).shape[0]), fmt='%i', delimiter=',') 
	savetxt('../../DataOutput/forecasts/GSFC_Petty_init'+str(m)+str(p)+region_str+iceType+'.txt', np.asarray(forecasts).reshape(1, np.asarray(forecasts).shape[0]),  fmt='%2.3f', delimiter=',') 
	savetxt('../../DataOutput/forecasts/GSFC_Petty_sep_obs_'+iceType+str(m)+str(p)+region_str+'.txt', np.asarray(extents).reshape(1, np.asarray(extents).shape[0]),  fmt='%2.3f', delimiter=',') 














