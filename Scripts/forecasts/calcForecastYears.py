""" 
 get_multivar_forecast_months_exp.py
 Date: 02/03/2017
 Author: Alek Petty
 Description: Script to generate sea ice extent forecast using a given variable (or multiple variables) for a given predictor and predictand month

 Run with e.g. python -i get_multivar_forecast_months_exp.py 
 If you want to run with options, uncomment out the parser lines, then run as:
 python -i get_multivar_forecast_months_exp.py  --fmonth 6 --pmonth 9 --fvars conc --weight 1

"""
import matplotlib
matplotlib.use("AGG")

import sys
sys.path.append('../')
import forecast_funcs as ff
from pylab import *

def plotForecastOneYear(figPath, years, extent, year, forecastVars, outVarStr, iceType):
	"""Plot forecast data """
	rcParams['xtick.major.size'] = 2
	rcParams['ytick.major.size'] = 2
	rcParams['axes.linewidth'] = .5
	rcParams['lines.linewidth'] = .5
	rcParams['patch.linewidth'] = .5
	rcParams['axes.labelsize'] = 8
	rcParams['xtick.labelsize']=8
	rcParams['ytick.labelsize']=8
	rcParams['legend.fontsize']=8
	rcParams['font.size']=7
	rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
    

	fig = figure(figsize=(3.5,2.2))
	ax1=subplot(1, 1, 1)
	im1 = plot(years, extent, 'k')
	#im2 = plot(Years[start_year_pred-start_year:], lineT[start_year_pred-start_year:]+ExtentG, 'r')
	
	im3 = plot(years[-1], extent[-1], marker='o', markersize=2, color='k')
	im3 = plot(year, forecastVars[2], marker='x', markersize=2, color='k')
	im3 = plot(year, forecastVars[3], marker='o', markersize=2, color='r')
	#errorbar(YearsP, array(lineTP)+array(ExtentG) , yerr=prederror, color='r',fmt='',linestyle='',lw=0.4,capsize=0.5, zorder = 2)
	#if (np.isfinite(forecastVars[4])):
	
	ax1.errorbar(year, forecastVars[3] , yerr=forecastVars[6], color='r',fmt='',linestyle='',lw=0.6,capsize=0.5, zorder = 2)
	#ax1.errorbar(yearsP, extentPredAbs , yerr=[1.96*x for x in perr], color='r',fmt='',linestyle='',lw=0.3,capsize=0.5, zorder = 2)

	forecastStr='%.2f' %(forecastVars[3])
	observedStr='%.2f' %(extent[-1])

	ax1.annotate('Year: '+str(year)+'\nObserved: '+observedStr+r' M km$^2$',
 		xy=(0.7, 1.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')

	ax1.annotate('\nForecast: '+forecastStr+r' M km$^2$',
 		xy=(0.7, 0.9), xycoords='axes fraction', color='r', horizontalalignment='left', verticalalignment='top')

	ax1.annotate('June forecasts of September sea ice / @alekpetty / alekpetty.com', fontsize=5, 
 		xy=(0.02, 0.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

	ax1.set_ylabel(iceType+r' (Million km$^2$)')
	#ax1.set_xlabel('Years')
	ax1.set_xlim(1980, 2020)
	ax1.set_xticks(np.arange(1980, 2021, 10))
	ax1.set_xticks(np.arange(1980, 2021, 5), minor=True)
	#ax1.set_xticklabels([])
	ax1.set_ylim(3, 8)

	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)

	plt.tight_layout()
	#subplots_adjust(left=0.15, right=0.90, bottom=0.17, top=0.96, hspace=0)

	savefig(figPath+'/forecast'+outVarStr+'.png', dpi=300)
	close(fig)



def main(year, fmonth, pmonth, fvars=['conc'], iceType='extent', hemStr='N', siiVersion='v3.0', startYear=1980, weight=1, region=0, numYearsReq=5, plotForecast=1, plotSkill=1, outSkill=1, outLine=1, outWeights=1):
	""" 
	Main forecast script:

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
	extentObsDT, extentyr, extentForrDT, extentForrAbs, anom, prstd[0]
		A dumped Python array including the folowing variables:
			Observed ice extent
			Observed detrended ice extent from linear trend persistence (LTP)
			LTP extent
			Absolute forecast of sea ice extent (added LTP)
			Detrended (forecast) ice extent from LTP
			Forecast anomaly from observed
			Estimated forecast error (1 SD)

	"""
	
	rawDataPath = '../../Data/' 
	derivedDataPath = '../../DataOutput/'
	

	if (hemStr=='S'):
		saveDataPath=derivedDataPath+'/Antarctic/'
		figPath='../../Figures/'+'/Antarctic/'
	elif (hemStr=='N'):
		saveDataPath=derivedDataPath+'/Arctic/'
		figPath='../../Figures/'+'/Arctic/'
	elif (region=='A'):
		saveDataPath=derivedDataPath+'/Alaska/'
		figPath='../../Figures/'+'/Alaska/'

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
			region=region, anomObs=1, numYearsReq=numYearsReq, weight=weight, outWeights=outWeights, 
			icetype=iceType, hemStr=hemStr, siiVersion=siiVersion)
	
	print ('Observed ice state:', forecastVals[0])
	print ('Linear trend presistence:', forecastVals[2])
	print ('Forecast ice state:', forecastVals[3])
	print ('Detrended observed ice state :', forecastVals[1])
	print ('Detrended observed ice state :', forecastVals[4])
	print ('Forecast anomaly :', forecastVals[5])
	
	array(forecastVals).dump(saveDataPath+outStr)

	if (plotForecast==1):
		years, extent = ff.get_ice_extentN(rawDataPath, pmonth, startYear, year, 
				icetype=iceType, version=siiVersion, hemStr=hemStr)

		plotForecastOneYear(figPath, years, extent, year, forecastVals, outStr, iceType)


#-- run main program
if __name__ == '__main__':
	#main(2015, 6, 9)
	for y in range(1990, 2017+1, 1):
		main(y, 6, 9)
	#	for m in range(startMonth, endMonth+1):
	#		print (y, m)
	#		
				













