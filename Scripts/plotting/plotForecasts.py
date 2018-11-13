"""  plotForecasts.py
	Script to plot the various forecast scripts
	Run with e.g. python plotForecasts.py 
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


def main(endYear, fmonth, pmonth, fvars=['conc'], iceType='extent', hemStr='N', siiVersion='v3.0', startYear=1979, minval=3, maxval=8, region=0):
	
	monthStrs=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'June', 'July', 'August', 'September']
	textStr=monthStrs[fmonth-1]+' forecasts of '+monthStrs[pmonth-1]+' Arctic sea ice '+iceType

	repoPath='/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/'
	rawDataPath = repoPath+'/Data/' 
	derivedDataPath = repoPath+'/DataOutput/'
	if (hemStr=='S'):
		saveDataPath=derivedDataPath+'/Antarctic/'
		figPath=repoPath+'/Figures/Forecasts/'
	elif (region=='A'):
		saveDataPath=derivedDataPath+'/Alaska/'
		figPath=repoPath+'/Figures/Forecasts/'

	else:
		saveDataPath=derivedDataPath+'/Arctic/'
		figPath=repoPath+'/Figures/Forecasts/'
	
	
	startYearP=1990
	#endYear=2017
	weight=1

	varStrsOut=''.join(fvars)

	forecastVarsM=np.array([]).reshape(0,7)


	for year in range(startYearP, endYear+1):
		outStr='forecastDump'+iceType+varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(year)+'W'+str(weight)+'SII'+siiVersion

		forecastVars=load(saveDataPath+outStr)

		forecastVarsM=np.vstack([forecastVarsM, forecastVars])

	if (region==0):
		years, extent = ff.get_ice_extentN(rawDataPath, pmonth, startYear, year, 
			icetype=iceType, version=siiVersion, hemStr=hemStr)
	elif (region=='A'):
		poleStr='A'

		extent=loadtxt(derivedDataPath+'/Extent/'+'ice_'+iceType+'_M'+str(pmonth)+'R'+str(region)+'_'+str(startYear)+'2017'+poleStr)
		#extent=extent[0:year-startYear+1]
		years=np.arange(startYear, 2017+1, 1)



	yearsP=np.arange(startYearP, endYear+1)

	#extentyr, extentObsDT, extTrendP, extentForrAbs, extentForrDT, anom, prstd[0]
	print ('Prediction int:', forecastVars[-1])

	skill = '%.2f' %(1 - (ff.rms(forecastVarsM[0:-1, 5]))/(ff.rms(forecastVarsM[0:-1, 1])))

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


	fig = figure(figsize=(4.,2.2))
	ax1=subplot(1, 1, 1)
	im1 = plot(years, extent, 'k')
	#im2 = plot(Years[start_year_pred-start_year:], lineT[start_year_pred-start_year:]+ExtentG, 'r')

	#im3 = plot(years[-1], extent[-1], marker='o', markersize=2, color='k')
	im3 = plot(yearsP, forecastVarsM[:, 2], marker='x', markersize=1, alpha=0.5, color='k')
	im3 = plot(yearsP, forecastVarsM[:, 3], marker='o', markersize=1, alpha=0.5, color='r')
	#ax1.errorbar(yearsP, forecastVarsM[:, 3] , yerr=forecastVarsM[:, 6], color='r',fmt='',linestyle='',alpha=0.5, lw=0.6,capsize=0.5, zorder = 2)

	im3 = plot(yearsP[-1], forecastVarsM[-1, 2], marker='x', markersize=2, color='k')
	im3 = plot(yearsP[-1], forecastVarsM[-1, 3], marker='o', markersize=2, color='r')
	ax1.errorbar(yearsP[-1], forecastVarsM[-1, 3] , yerr=forecastVarsM[-1, 6], color='r',fmt='',linestyle='',lw=0.6,capsize=0.5, zorder = 2)

	#errorbar(YearsP, array(lineTP)+array(ExtentG) , yerr=prederror, color='r',fmt='',linestyle='',lw=0.4,capsize=0.5, zorder = 2)
	#if (np.isfinite(forecastVars[4])):

	#ax1.errorbar(yearsP, extentPredAbs , yerr=[1.96*x for x in perr], color='r',fmt='',linestyle='',lw=0.3,capsize=0.5, zorder = 2)

	forecastStr='%.2f' %(forecastVarsM[-1, 3])
	linearStr='%.2f' %(forecastVarsM[-1, 2])
	#observedStr='%.2f' %(extent[-1])

	ax1.annotate('Year: '+str(year)+'\nLinear trend forecast: '+linearStr+r' M km$^2$',
			xy=(0.8, 1.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')

	ax1.annotate('NASA GSFC forecast: '+forecastStr+r' M km$^2$',
			xy=(0.8, 0.9), xycoords='axes fraction', color='r', horizontalalignment='left', verticalalignment='top')

	#ax1.annotate('Year: '+str(year)+'\nObserved: '+observedStr+r' M km$^2$'+'\nTrend: '+linearStr+r' M km$^2$',
	#		xy=(1., 1.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')

	#ax1.annotate(,
	#		xy=(0.8, 0.85), xycoords='axes fraction', color='r', horizontalalignment='left', verticalalignment='top')

	#ax1.annotate('NASA GSFC forecast: '+forecastStr+r' M km$^2$'+'\nSkill:'+skill,
	#		xy=(0.8, 0.85), xycoords='axes fraction', color='r', horizontalalignment='left', verticalalignment='top')

	ax1.annotate(textStr, fontsize=5, 
			xy=(0.02, 0.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

	ax1.set_ylabel('Ice '+iceType+r' (Million km$^2$)')
	#ax1.set_xlabel('Years')
	ax1.set_xlim(1980, 2020)
	ax1.set_xticks(np.arange(1980, 2021, 10))
	ax1.set_xticks(np.arange(1980, 2021, 5), minor=True)
	#ax1.set_xticklabels([])
	ax1.set_ylim(minval, maxval-0.1)
	ax1.set_yticks(np.arange(minval, maxval, 1))

	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)

	ax1.yaxis.grid(which='major', linestyle='-', linewidth='0.1', color='k')
	#plt.tight_layout()
	#subplots_adjust(left=0.15, right=0.90, bottom=0.17, top=0.96, hspace=0)
	subplots_adjust(left=0.11, right=0.74, bottom=0.1, top=0.96, hspace=0)

	savefig(figPath+'/'+outStr+hemStr+'multi.png', dpi=300)
	savefig(figPath+'/'+outStr+hemStr+'multi.jpg', dpi=300)
	close(fig)


#-- run main program
if __name__ == '__main__':
	#main(2015, 6, 9)
	for y in range(2018, 2018+1, 1):
		main(y, 7, 9, iceType='extent', hemStr='N', minval=2, maxval=8, region=0, startYear=1979)


