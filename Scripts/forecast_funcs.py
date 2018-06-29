"""  forecast_funcs.py
	Various functions used by the forecast repo
	Author:
		Alek Petty

	Update history:
		04/20/2018: Version 1
"""

import matplotlib
matplotlib.use("AGG")
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from pylab import *
import numpy.ma as ma
from glob import glob
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
#from netCDF4 import Dataset


# rawdatapath='../Data/'
# datapath='../DataOutput/'

# skilldatapath='../DataOutput/SkillVals/'
# linedatapath='../DataOutput/TimeSeries/'
# weightdatapath='../DataOutput/Weights/'
# extdatapath = rawdatapath+'/IceExtent/'

def rms(var):
	"""calculate th root mean square of a given list """
	return sum([x**2 for x in var])/size(var)
                                                              


def get_detrended_yr(yearsTr, yearT, var_yearsT, var_yrT, num_years_req):
	"""Detrend a 2D array using linear regression

       Mask based on valid number of years in each grid cell.
    """
	var_yearsDT=ma.masked_all((var_yearsT.shape))
	var_yrDT=ma.masked_all((var_yrT.shape))

	# Loop over each dimension
	for i in range(var_yearsT.shape[1]):
		for j in range(var_yearsT.shape[2]):
			mask=~var_yearsT[:, i, j].mask
			var_yearsT_ma = var_yearsT[:, i, j][mask]	
				
			if (len(var_yearsT_ma)>num_years_req):
				trendT, interceptT, r_valsT, probT, stderrT = stats.linregress(yearsTr[mask],var_yearsT_ma)
				lineT = (trendT*yearsTr) + interceptT
				#print var_yearsT[:, i, j].shape, lineT.shape, yearsTr
				var_yearsDT[:, i, j]=var_yearsT[:, i, j]-lineT
				
				# Calculate the detrended var (linear trend persistence) fo the given forecast year
				lineT_yr=interceptT + (trendT*(yearT))
				var_yrDT[i, j]=var_yrT[i, j]-lineT_yr

	return var_yearsDT, var_yrDT

def plotForecastOneYear(figPath, years, extent, year, forecastVars, outVarStr, iceType, minval=0, maxval=10):
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
	
	im3 = plot(years[-1], extent[-1], marker='o', markersize=2, color='k')
	im3 = plot(year, forecastVars[2], marker='x', markersize=2, color='k')
	im3 = plot(year, forecastVars[3], marker='o', markersize=2, color='r')
	#errorbar(YearsP, array(lineTP)+array(ExtentG) , yerr=prederror, color='r',fmt='',linestyle='',lw=0.4,capsize=0.5, zorder = 2)
	#if (np.isfinite(forecastVars[4])):
	
	ax1.errorbar(year, forecastVars[3] , yerr=forecastVars[6], color='r',fmt='',linestyle='',lw=0.6,capsize=0.5, zorder = 2)
	#ax1.errorbar(yearsP, extentPredAbs , yerr=[1.96*x for x in perr], color='r',fmt='',linestyle='',lw=0.3,capsize=0.5, zorder = 2)

	forecastStr='%.2f' %(forecastVars[3])
	linearStr='%.2f' %(forecastVars[2])
	observedStr='%.2f' %(extent[-1])

	ax1.annotate('Year: '+str(year)+'\nObserved: '+observedStr+r' M km$^2$'+'\nTrend: '+linearStr+r' M km$^2$',
			xy=(1., 1.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')

	#ax1.annotate('Year: '+str(year)+'\nObserved: '+observedStr+r' M km$^2$',
 	#	xy=(1., 1.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top')

	ax1.annotate('\nForecast: '+forecastStr+r' M km$^2$',
 		xy=(1., 0.85), xycoords='axes fraction', color='r', horizontalalignment='left', verticalalignment='top')

	#ax1.annotate('June forecasts of September sea ice / @alekpetty / alekpetty.com', fontsize=5, 
 	#	xy=(0.02, 0.02), xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

	ax1.set_ylabel(iceType+r' (Million km$^2$)')
	#ax1.set_xlabel('Years')
	ax1.set_xlim(1980, 2020)
	ax1.set_xticks(np.arange(1980, 2021, 10))
	ax1.set_xticks(np.arange(1980, 2021, 5), minor=True)
	#ax1.set_xticklabels([])
	ax1.set_ylim(minval, maxval)

	ax1.spines['right'].set_visible(False)
	ax1.spines['top'].set_visible(False)

	#plt.tight_layout()
	subplots_adjust(left=0.13, right=0.75, bottom=0.1, top=0.96, hspace=0)

	savefig(figPath+'/forecast'+outVarStr+'.png', dpi=300)
	close(fig)


def get_psnlatslons(data_path):
	""" Get Arctic polar stereographic lon/lats from the NSIDC"""
	mask_latf = open(data_path+'/OTHER/psn25lats_v3.dat', 'rb')
	mask_lonf = open(data_path+'/OTHER/psn25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	return lats_mask, lons_mask

def get_psslatslons(data_path):
	""" Get Antarctic polar stereographic lon/lats from the NSIDC"""
	mask_latf = open(data_path+'/OTHER/pss25lats_v3.dat', 'rb')
	mask_lonf = open(data_path+'/OTHER/pss25lons_v3.dat', 'rb')
	lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [332, 316])
	lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [332, 316])

	return lats_mask, lons_mask

def get_month_concSN(datapath, year, month, alg=0, pole='AA',  mask=1, maxConc=0, lowerConc=0):
	""" Get monthly ice concentration from the NSIDC"""
	if (alg==0):
		team = 'NASA_TEAM'
		team_s = 'nt'
		header = 300
		datatype='uint8'
		scale_factor=250.
	if (alg==1):
		team = 'BOOTSTRAP'
		team_s = 'bt'
		header = 0
		datatype='<i2'
		scale_factor=1000.

	if (pole=='AA'):
		poleStr='ANTARCTIC'
		rows=332
		cols=316

	if (pole=='A'):
		poleStr='ARCTIC'
		rows=448
		cols=304

	month_str = '%02d' % (month+1)
	year_str=str(year)
	files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/monthly/'+team_s+'_'+year_str+month_str+'*.bin')
	fd = open(files[-1], 'r')
	data = fromfile(file=fd, dtype=datatype)
	data = data[header:]
	#FIRST 300 FILES ARE HEADER INFO
	ice_conc = reshape(data, [rows, cols])
	#divide by 250 to express in concentration
	ice_conc = ma.masked_where(ice_conc>250., ice_conc)
	ice_conc = ice_conc/scale_factor
	#GREATER THAN 250 is mask/land etc
	if (mask==1):
		ice_conc = ma.masked_where(ice_conc>1., ice_conc)
	
	if (maxConc==1):
		ice_conc = ma.where(ice_conc>1.,0, ice_conc)

	if (lowerConc==1):
		ice_conc = ma.where(ice_conc<0.15,0, ice_conc)


	return ice_conc


def get_month_concSN_NRT(datapath, year, month, alg=0, pole='A',  mask=1, maxConc=0, lowerConc=0, monthMean=0):
	""" Get near real-time (NRT) monthly ice concentration from the NSIDC, generated from daily data"""
	if (alg==0):
		team = 'NASA_TEAM'
		team_s = 'nt'
		header = 300
		datatype='uint8'
		scale_factor=250.
	if (alg==1):
		team = 'BOOTSTRAP'
		team_s = 'NH'
		header = 0
		datatype='<i2'
		scale_factor=1000.
	
	if (pole=='A'):
		poleStr='ARCTIC'
		rows=448
		cols=304
	if (pole=='AA'):
		poleStr='ANTARCTIC'
		rows=332
		cols=316

	month_str = '%02d' % (month+1)
	year_str=str(year)
	files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+'*')
	
	print ('Num conc files:', size(files), 'in month:'+month_str)
	ice_conc = ma.masked_all((size(files), rows, cols))
	
	for x in range(size(files)):
		fd = open(files[x], 'r')
		data = fromfile(file=fd, dtype=datatype)
		data = data[header:]
		#FIRST 300 FILES ARE HEADER INFO
		ice_conc[x] = reshape(data, [rows, cols])
		
	#divide by 250 to express in concentration
	ice_conc = ice_conc/scale_factor
	#GREATER THAN 250 is mask/land etc
	
	if (mask==1):
		ice_conc = ma.masked_where(ice_conc>1., ice_conc)
	
	if (maxConc==1):
		ice_conc = ma.where(ice_conc>1.,0, ice_conc)

	if (lowerConc==1):
		ice_conc = ma.where(ice_conc<0.15,0, ice_conc)

	if (monthMean==1):
		ice_conc=ma.mean(ice_conc, axis=0)

	return ice_conc
	



def plot_forecast(outVarStr):
	"""Plot forecast data """
    
	fig = figure(figsize=(3.5,2.2))
	ax1=subplot(2, 1, 1)
	im1 = plot(years, extent, 'k')
	#im2 = plot(Years[start_year_pred-start_year:], lineT[start_year_pred-start_year:]+ExtentG, 'r')
	im3 = plot(yearsP, extentPredAbs, 'r')
	#errorbar(YearsP, array(lineTP)+array(ExtentG) , yerr=prederror, color='r',fmt='',linestyle='',lw=0.4,capsize=0.5, zorder = 2)
	ax1.errorbar(yearsP, extentPredAbs , yerr=perr, color='r',fmt='',linestyle='',lw=0.6,capsize=0.5, zorder = 2)
	ax1.errorbar(yearsP, extentPredAbs , yerr=[1.96*x for x in perr], color='r',fmt='',linestyle='',lw=0.3,capsize=0.5, zorder = 2)

	ax1.set_ylabel(r'Extent (M km$^2$)')
	ax1.set_xlim(1978, 2017)
	ax1.set_xticks(np.arange(1980, 2018, 5))
	ax1.set_xticklabels([])
	#ylim(3, 9)

	ax2=subplot(2, 1, 2)
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	im21 = plot(yearsP[0:-1], extentObsDt, 'k')
	im3 = plot(yearsP, extentPredDt, 'r')
	ax2.errorbar(yearsP, extentPredDt , yerr=perr, color='r',fmt='',linestyle='',lw=0.6,capsize=0.5, zorder = 2)
	ax2.errorbar(yearsP, extentPredDt , yerr=[1.96*x for x in perr], color='r',fmt='',linestyle='',lw=0.3,capsize=0.5, zorder = 2)


	ax2.set_ylabel(r'Extent anomaly (M km$^2$)', rotation=270, labelpad=10)
	ax2.set_xlabel('Years')
	ax2.set_yticks([-2, -1, 0, 1, 2])
	ax2.set_xlim(1978, 2017)
	ax2.set_xticks(np.arange(1980, 2018, 5))
	ax2.axhline(0, linewidth=0.5,linestyle='--', color='k')
	ax2.annotate(r'$\sigma_{ferr}$='+errorFore+r' M km$^2$'+', S:'+skill, 
		xy=(0.03, 0.04), xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom')

	subplots_adjust(left=0.1, right=0.90, bottom=0.17, top=0.96, hspace=0)

	savefig(figpath+'/forecast'+str(startYear)+str(endYear)+'M'+str(fmonth)+outVarStr+'.pdf', dpi=300)
	close(fig)



def get_correlation_coeffs(var_yearsT, ExtentDT, num_years_req):
	""" Calculate the correlation coeficients between the detrended forecast
		variable and detrended ice extent
    """
	r_valsDT=np.zeros((var_yearsT.shape[1], var_yearsT.shape[2]))
	for i in range(var_yearsT.shape[1]):
		for j in range(var_yearsT.shape[2]):
			mask=~var_yearsT[:, i, j].mask
			var_yearsT_ma = var_yearsT[:, i, j][mask]
			if (len(var_yearsT_ma)>num_years_req):
				trendDT, interceptDT, r_valsDT[i, j], probDT, stderrDT = stats.linregress(ExtentDT[mask],var_yearsT_ma)
	return r_valsDT

def GetWeightedPredVar(deriveddatapath, yearsTr, yearT, extentDTT, predvarYrsT, predvar_yrT, varT, fmonth, pmonth, startYear, numYearsReq, region, hemStr, icetype, normalize=0, rneg=0, rpos=1, absr=1, weight=1, outWeights=0):
	""" Get forecast data and weight using historical correlation if selected
    """
	
	if (hemStr=='S'):
		savedatapath=deriveddatapath+'/Antarctic/'
	elif (hemStr=='N'):
		savedatapath=deriveddatapath+'/Arctic/'
	if (region=='A'):
		savedatapath=deriveddatapath+'/Alaska/'

		
	# Get detrended 2D forecast data
	predvarYrsDT, predvarYrDT = get_detrended_yr(yearsTr, yearT, predvarYrsT, predvar_yrT, numYearsReq)

	# Correlate detrended time series
	rvalsDT = get_correlation_coeffs(predvarYrsDT, extentDTT, numYearsReq)
	
	if (rneg==1):
		# Set positive R-vals to zero (assumed to be unphysical)
		rvalsDT[where(rvalsDT>0)]=0
	if (rpos==1):
		# Set negative R-vals to zero (assumed to be unphysical)
		rvalsDT[where(rvalsDT<0)]=0
	if (absr==1):
		# Use absolute values of correlation coefficeint
		rvalsDT=abs(rvalsDT)
	if (weight==0):
		print ('No weighting applied!')
		rvalsDT=np.ones((rvalsDT.shape))

	if (outWeights==1):
		rvalsDT.dump(savedatapath+'rvalsDT'+varT+icetype+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(yearT)+'.txt')
		predvarYrDT.dump(savedatapath+'predvarYrDT'+varT+icetype+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(yearT)+'.txt')
	
	# Calculated weighted forcast data
	weightedPredvar=[]
	for x in range(predvarYrsDT.shape[0]):
		weightedPredvar.append(ma.mean(rvalsDT*predvarYrsDT[x]))
	
	weightedPredvarYr = ma.mean(rvalsDT*predvarYrDT)
	
	if (normalize==1):
		# Normalize data (doesn't change single var forecasting, may be important for multivar)
		weightedPredvarN=(weightedPredvar-min(weightedPredvar))/(max(weightedPredvar)-min(weightedPredvar))
		weightedPredvarYrN=(weightedPredvarYr-min(weightedPredvar))/(max(weightedPredvar)-min(weightedPredvar))
		return rvalsDT, predvarYrDT, weightedPredvarN, weightedPredvarYrN
	else:
		return rvalsDT, predvarYrDT, weightedPredvar, weightedPredvarYr


def get_varDT(Years, Extent):
	""" Detrend linear time series  """
	trendT, interceptT, r_valsT, probT, stderrT = stats.linregress(Years,Extent)
	lineT = (trendT*Years) + interceptT
	ExtentDT=Extent-lineT
	return ExtentDT, lineT


def get_ice_extentN(extdatapath, Month, start_year, end_year, icetype='extent', version='', hemStr='N'):
	""" Get Arctic sea ice extent

	Data downlaoded from the NSIDC Arctic Sea Ice Index.

	Can also get ice area if icetype set to 'area', 
	   but beware of variable pole hole contaminating Arctic data

	"""
	Month_str = '%02d' %Month
	extent_data_path=extdatapath+'SeaIceIndex/'+hemStr+'_'+Month_str+'_extent_'+version+'.csv'
	ice_extent_data=pd.read_csv(extent_data_path,names=['year', 'extent', 'area'],skiprows=1, usecols=[0, 4, 5])
	#ice_extent_data = np.loadtxt(extent_data_path, delimiter=',',skiprows=1)
	Extent = ice_extent_data[icetype]
	Year = ice_extent_data['year']
	
	#print 'Y:', Year
	#Years=array(Year[start_year-1979:end_year-1979+1])
	Years=array(Year[(Year>=start_year)&(Year<=end_year)])
	Extents=array(Extent[(Year>=start_year)&(Year<=end_year)])

	Years=Years[where(Extents>0)]
	Extents=Extents[where(Extents>0)]

	#Extents=ma.masked_where(Extents<0, Extents)
	#Extent=array(Extent[start_year-1979:end_year-1979+1])

	return Years, Extents

def get_region_maskAA(datapath, mplot, xypts_return=0):
	header = 300

	datatype='uint8'
	file_mask = datapath+'/OTHER/region_s.msk'
	#1   non-region oceans
	#2 Weddell Sea
	#3 Indian Ocean
	#4 Pacific Ocean
	#5 Soss Sea
	#6 Amundsen/BHausen Sea
	#11 Land
	#12 Coast


	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	print (region_mask.shape)
	region_mask = reshape(region_mask[header:], [332, 316])


	if (xypts_return==1):
		mask_latf = open(datapath+'/OTHER/pss25lats_v3.dat', 'rb')
		mask_lonf = open(datapath+'/OTHER/pss25lons_v3.dat', 'rb')
		lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [332, 316])
		lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [332, 316])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask

def getIceExtentAreaPetty(dataOutPath, month, start_year, end_year, icetype='extent', alg=0, extraStr=''):
	""" Get Arctic sea ice extent using Petty/NSIDC method

	Data downloaded from the NSIDC and extent caluclated using the ASI

	Can also get ice area if icetype set to 'area', 
	   but beware of variable pole hole contaminating Arctic data

	"""

	mstr = '%02d' %month

	if (icetype=='area'):
		typeStr='Area'
	else:
		typeStr='Ext'
		
	Extent = loadtxt(dataOutPath+'ice'+typeStr+'Months'+str(start_year)+str(2016)+'-'+mstr+'Alg-'+str(alg)+extraStr)[:]
	#extent = extentT[start_year-1979:end_year-1979+1]
	Year=np.arange(start_year, 2016+1, 1)

	Years=array(Year[(Year>=start_year)&(Year<=end_year)])
	Extents=array(Extent[(Year>=start_year)&(Year<=end_year)])

	Years=Years[where(Extents>0)]
	Extents=Extents[where(Extents>0)]
	

	return Years, Extents

def CalcForecastMultiVar(rawdatapath, deriveddatapath, yearF, startYear, predvarYrs, fmonth, pmonth=9, region=0, anomObs=1 , outWeights=0, 
	icetype='extent', numYearsReq=5, weight=1, hemStr='N', siiVersion='v2.1'):
	""" The primary sea ice forecast function. 

	NB: This should probably be converted to a class at some point.

	"""

	# Initially set the year we are predicting to be the same as the year of the forecast data being used.
	# Note that this changes if we want to start using say December data and forecast the following spring
	yearP=0
	yearP=yearF
	# Thus we initialise the forecast year to be the same as the initial prediction year
	startYearF=0
	startYearF=startYear

	# However, if the month we are predicting is lower than the forecast month (i.e. December=12, January=1) 
	# then switch the initial forecast year back one.
	if (pmonth<fmonth):
		startYearF=startYear-1
		yearF=yearF-1
		print ('pmonth<fmonth:', yearF, yearP)
			


	if (region==0): 
	# use the NSIDC Arctic/Antarctic Sea Ice Index
		print (startYear, yearP-1, pmonth)

		# Get the NSIDC sea ice index data
		yearsPr, extentTr = get_ice_extentN(rawdatapath, pmonth, startYear, yearP-1, 
			icetype=icetype, version=siiVersion, hemStr=hemStr)

		# De-trend the extent data
		extentDTr, lineTr=get_varDT(yearsPr, extentTr)

		if (anomObs==1):
		# If we have observed sea ice extent data for the given forecast year to check forecast skill.
			print ('anomObs')
			# Get the NSIDC sea ice index data
			years2, extent2 = get_ice_extentN(rawdatapath, pmonth, startYear, yearP, 
				icetype=icetype, version=siiVersion,  hemStr=hemStr)

			# If the last year of the extent file doesn't match the prediction year dont compare with observed.
			if (years2[-1]!=yearP):
				anomObs=0

			# Get current year ice extent/area
			extentyr=extent2[-1]

	elif (region==-1): 
	# region ==-1 is Alek's pan-Arctic sea ice indicies (which includes filling the pole hole for SIA)
		print (startYear, yearP-1, pmonth)
		yearsPr, extentTr = getIceExtentAreaPetty(deriveddatapath, pmonth, startYear, yearP-1, icetype=icetype, alg=0)
		print (extentTr)

		extentDTr, lineTr=get_varDT(yearsPr, extentTr)
		
		if (anomObs==1):
		# If we have observed sea ice extent data for the given forecast year to check forecast skill.
			print ('anomObs')
			
			# Get the Petty sea ice index data up to forecast year
			years2, extent2 = getIceExtentAreaPetty(deriveddatapath, pmonth, startYear, yearP, icetype=icetype, alg=0)

			# If the last year of the extent file doesn't match the prediciton year dont compare with observed.
			if (years2[-1]!=yearP):
				print ('year check:', yearsPr[-1], yearP)
				anomObs=0
			# Get current year ice extent/area
			extentyr=extent2[-1]
	else: 
		
		# Get regional sea ice indices we generate
		if (hemStr=='N'):
			poleStr='A'
		elif (hemStr=='S'):
			poleStr='AA'

		extentALL=loadtxt(deriveddatapath+'/Extent/'+'ice_'+icetype+'_M'+str(pmonth)+'R'+str(region)+'_'+str(startYearF)+'2017'+poleStr)
		

		#get years and extent for years preceeding the given forecast year
		yearsPr=np.arange(startYear, yearP, 1)
		extentTr=extentALL[0:yearP-startYear]

		# De-trend the extent data
		extentDTr, lineTr=get_varDT(yearsPr, extentTr)
		
		if (anomObs==1):
			# Get current year ice extent/area
			extentyr=extentALL[yearP-startYear]
	
	# Need to get an array filled with ones to act as the intercept
	predVarsTYr=[1]
	predVars=np.ones((size(yearsPr)))

	# Needed for melt pond forecast
	if (fmonth>=6):
		# June
		pdate='56'
	else:
		# May
		pdate='31'

	# Get forecast years
	yearsFr=np.arange(startYearF, yearF, 1)

	print ('Training years', yearsPr)
	print ('Predicted year', yearP)
	print ('Forecast year', yearF)
	print ('Training start year', startYear)
	print ('Training start year', startYearF)
	print ('Forecast years', yearsFr)

	
	for varT in predvarYrs:
		#print 'Var:', varT
		if (varT in ['sst','conc','melt','melt_nan', 'pmas']):

			# Get the gridded forecast data for training
			VarYearsTr = get_gridvar(deriveddatapath, varT, fmonth, yearsFr, hemStr)

			# Get the gridded prediction data
			VarYear = get_gridvar(deriveddatapath, varT, fmonth, array(yearF), hemStr)
			
			# Weight the gridded forecast data with historical sea ice extent
			rvalsDT, unweightedpredVarT, predVarT, predVarTYr = GetWeightedPredVar(deriveddatapath, yearsFr, yearF, extentDTr, VarYearsTr, VarYear,varT, fmonth, pmonth, startYearF,numYearsReq, region, hemStr, icetype, normalize=0, outWeights=outWeights, weight=weight)
		
		# will be an array of 1 (intercept) and a number
		predVarsTYr.append(predVarTYr)
		# will be an array of 1s (intercepts) and a series of numbers
		predVars=np.column_stack((predVars, array(predVarT)))
	

	# Use SM to generate the regression model. Could have just used linregress (tested, gave same results, but this was just a bit neater)
	model=sm.OLS(extentDTr, predVars)
	fit=model.fit()

	# Forecast detrended sea ice extent!
	extentForrDT = fit.predict(predVarsTYr)[0]
	# Prediction uncertainty estimate
	prstd, iv_l, iv_u = wls_prediction_std(fit, exog=predVarsTYr)

	# Calculate ice extent assuming inear trend persistnce
	extTrendP=(lineTr[-1]+(lineTr[-1]-lineTr[-2]))

	extentForrAbs = extentForrDT+extTrendP
	
	#print ('detrended extent forecast :', extentForrDT, 'Linear trend extent:',extTrendP, 'Absolute extent forecast:',extentForrAbs)
	
	#print ('Did we produce an observed anomaly:', anomObs)

	if (anomObs==1):
		extentObsDT=extentyr-extTrendP
		anom=extentyr-extentForrAbs
		return  extentyr, extentObsDT, extTrendP, extentForrAbs, extentForrDT, anom, prstd[0]
	else:
		
		return  np.nan, np.nan, extTrendP, extentForrAbs, extentForrDT, np.nan, prstd[0]



def get_pmask(year, month):
	#remove half a degree as gridding around the pole hole edge
	if (year<1987):
		pmask=84.4
	elif((year==1987)&(month<=7)):
		pmask=84.4
	elif ((year==1987)&(month>7)):
		pmask=86.7
	elif ((year>1987)&(year<2008)):
		pmask=87.2
	else:
		pmask=89.2
	
	return pmask

def get_pmas_month(m, rawdatapath, year, month=4):

	fd = open(rawdatapath+'/PIOMAS/heff/heff.H'+str(year), 'rb')
	dataP = fromfile(file=fd, dtype='f')
	dataP = reshape(dataP, [12, 120*360])
	thickness=dataP[month]
	gridP = loadtxt(rawdatapath+'/PIOMAS/grid.dat.txt')

	lonsP = gridP[0:4320, :].flatten()
	latsP = gridP[4320:, :].flatten()
	xptsP,yptsP = m(lonsP, latsP)

	thickness=ma.masked_where(thickness<0.01, thickness)

	return xptsP, yptsP, thickness
	
def get_region_mask_sect(datapath, mplot, xypts_return=0):
	# Get NSIDC region masks
	datatype='uint8'
	file_mask = datapath+'/OTHER/sect_fixed_n.msk'
	# 1   non-region oceans
	# ;           = 2   Sea of Okhotsk and Japan
	# ;           = 3   Bering Sea
	# ;           = 4   Hudson Bay
	# ;           = 5   Gulf of St. Lawrence
	# ;           = 6   Baffin Bay/Davis Strait/Labrador Sea
	# ;           = 7   Greenland Sea
	# ;           = 8   Barents Seas
	# ;           = 9   Kara
	# ;           =10   Laptev
	# ;           =11   E. Siberian
	# ;           =12   Chukchi
	# ;           =13   Beaufort
	# ;           =14   Canadian Archipelago
	# ;           =15   Arctic Ocean
	# ;           =20   Land
	# ;           =21   Coast
	fd = open(file_mask, 'rb')
	region_mask = fromfile(file=fd, dtype=datatype)
	region_mask = reshape(region_mask, [448, 304])

	#mask_latf = open('/Volumes/GRAID_NASA/NOAA/DATA/ICE_CONC/BOOTSTRAP/psn25lats_v3.dat', 'rb')
	#mask_lonf = open('/Volumes/GRAID_NASA/NOAA/DATA/ICE_CONC/BOOTSTRAP/psn25lons_v3.dat', 'rb')
	#lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	#lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	#xpts, ypts = mplot(lons_mask, lats_mask)
	if (xypts_return==1):
		mask_latf = open(datapath+'/OTHER/psn25lats_v3.dat', 'rb')
		mask_lonf = open(datapath+'/OTHER/psn25lons_v3.dat', 'rb')
		lats_mask = reshape(fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
		lons_mask = reshape(fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

		xpts, ypts = mplot(lons_mask, lats_mask)

		return region_mask, xpts, ypts
	else:
		return region_mask


def get_conc_gridded(dataoutpath, yearsT, month, hemStr, concVersion='v2'):
	""" Get gridded ice concentration data

	Data gridded using linear interpolation of NASA Team concentration data onto a 100 km grid.
	Used monthly data, then monthly means of the daily NRT data for 2015 onwards.


	"""

	if (hemStr=='N'):
		poleStr='A'
	elif (hemStr=='S'):
		poleStr='AA'


	xpts=load(dataoutpath+concVersion+'xpts100km'+poleStr)
	ypts=load(dataoutpath+concVersion+'ypts100km'+poleStr)

	if (size(yearsT)>1):
		conc_years=ma.masked_all((size(yearsT),xpts.shape[0], xpts.shape[1]))
		x=0
		for year in yearsT:
			conc_years[x] = load(dataoutpath+concVersion+'ice_conc100km'+str(month)+str(year)+poleStr+concVersion)
			x+=1
	else:
		conc_years = load(dataoutpath+concVersion+'ice_conc100km'+str(month)+str(yearsT)+poleStr+concVersion)

	return xpts, ypts, conc_years


def get_meltonset_gridded(dataoutpath, yearsT, freezemelt_str, hemStr):
	""" Get gridded melt onset data

	Data gridded using linear interpolation of NASA's GSFC melt onset data onto a 100 km grid.

	"""

	if (hemStr=='N'):
		poleStr='A'
	elif (hemStr=='S'):
		poleStr='AA'

	xpts=load(dataoutpath+'xpts100km'+poleStr)
	ypts=load(dataoutpath+'ypts100km'+poleStr)
	Melt_onset_years=ma.masked_all((size(yearsT),xpts.shape[0], xpts.shape[1]))
	x=0
	if (size(yearsT)>1):
		Melt_onset_years=ma.masked_all((size(yearsT),xpts.shape[0], xpts.shape[1]))
		x=0
		for year in yearsT:
			Melt_onset_years[x] = load(dataoutpath+freezemelt_str+'100km'+str(year)+poleStr)
			x+=1
	else:
		Melt_onset_years = load(dataoutpath+freezemelt_str+'100km'+str(yearsT)+poleStr)

	return xpts, ypts, Melt_onset_years


def get_gridvar(griddatapath, fvar, fmonth, yearsT, hemStr, concVersion=''):
	""" Select which gridded forecast dataset to use in forecast

	NB pond data left out for now.

	"""
	# SUBTRACT ONE FROM FORECAST MONTH TO START MONTH INDEX AT 0.
	fmonth=fmonth-1

	if (fvar=='conc'):
		if (hemStr=='N'):
			dataoutpath=griddatapath+'IceConcA/'
		elif (hemStr=='S'):
			dataoutpath=griddatapath+'IceConcAA/'

		xpts, ypts, VarYears=get_conc_gridded(dataoutpath, yearsT, fmonth, hemStr, concVersion=concVersion)
		#rneg=0
		#rpos=1
	if ((fvar=='melt')|(fvar=='melt_nan')):
		meltdays=[31, 59, 90, 120, 151, 181, 212, 243]
		
		meltday=meltdays[fmonth]
		dataoutpath=datapath+'/MeltOnset/'
		xpts, ypts, VarYears=get_meltonset_gridded(dataoutpath, yearsT, fvar, hemStr)
		# Express melt onset relative to the given forecast date (end of the forecast month)
		VarYears=meltday-VarYears
		VarYears[where(VarYears<0)]=0
		#reverse to make consistent with concentration - i.e. low vals lead to low ice extent
		VarYears=-VarYears
	return VarYears


