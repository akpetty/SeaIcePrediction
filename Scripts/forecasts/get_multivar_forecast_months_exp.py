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

import argparse
import sys
import forecast_funcs as ff
from pylab import *

figpath='../Figures/'

#rawdatapath='../../../../DATA/'
rawdatapath = '/Volumes/PETTY_PASSPORT2/DATA/' 
extdatapath = '../DataOutput/Extent/'

# Options for this Script
plotSkill=1 # for plotting the skill
outSkill=1 # far saving the skill values
outLine=1 # for saving the forecast time series
outWeights=1 #for saving the weightings


startYear=1980 # Start of forecast training
endYear=2018
startYearPred=1985 # Start of forecasts

yearsP=arange(startYearPred, endYear+1, 1)
numYearsReq=5 #number of years required in a grid cell for it to count towards the training data

# Defaults if no arguments added after the call.
fmonth=11 #6=June, 9=Sep #  Forecast month
pmonth=2 #9=SEP # Predicted month
fvars=['conc']
weight=1 # Spatially weighting the data
hemStr='S' # Hemipshere (N or S)
iceType='area' # ice type being forecast (extent or area)


region=0



if (hemStr=='S'):
	skilldatapath='../DataOutput/Antarctic/SkillVals/'
	linedatapath='../DataOutput/Antarctic/TimeSeries/'
	datapath='../DataOutput/Antarctic/'
elif (hemStr=='N'):
	skilldatapath='../DataOutput/Arctic/SkillVals/'
	linedatapath='../DataOutput/Arctic/TimeSeries/'
	datapath='../DataOutput/Arctic/'

if (region=='A'):
	skilldatapath='../DataOutput/Alaska/SkillVals/'
	linedatapath='../DataOutput/Alaska/TimeSeries/'


#0 implies pan-Arctic or Antarctic
#2 Weddell Sea
#3 Indian Ocean
#4 Pacific Ocean
#5 Ross Sea
#6 Amundsen/BHausen Sea
#A Alaskan


#parser = argparse.ArgumentParser()
#parser.add_argument('--fmonth', type=int)
#parser.add_argument('--pmonth', type=int)
#parser.add_argument('--fvars', nargs='+')
#parser.add_argument('--weight', type=int)
#parser.add_argument('--numVars', type=int)
#args=parser.parse_args()
#fvars=args.fvars
#weight=args.weight
#fmonth=args.fmonth
#pmonth=args.pmonth

varStrsOut=''.join(fvars)

print 'fmonth:', fmonth, 'pmonth:', pmonth
print fvars
print 'Weighted:', weight


# Load empty lists (not very Pythony)
extentObsDt=[] # Detrended (observed) ice extent from linear trend persistence (LTP)
extentPredDt=[] # Detrended (forecast) ice extent from LTP
extentPredAbs=[] # Aboslute forecast of sea ice extent (add LTP)
anoms=[] # Forecast anomaly from observed
perr=[] # Estimated forecast error (1 SD)

# Generate a new forecast each year
for year in xrange(startYearPred, endYear+1, 1):
	if (year<endYear):
		vals=ff.CalcForecastMultiVar(rawdatapath, extdatapath, year, startYear, fvars, fmonth, pmonth=pmonth,
			region=region, anomObs=1, numYearsReq=numYearsReq, weight=weight, outWeights=outWeights, 
			icetype=iceType, hemStr=hemStr, siiVersion='v3.0')
		print 'Forecast year:', year, fvars, 'DT SIE forecast:', vals[1], 'DT SIE Obs :', vals[0],'SIE forecast:', vals[2] 
		
		extentObsDt.append(vals[0])
		extentPredDt.append(vals[1])
		extentPredAbs.append(vals[2])
		anoms.append(vals[3])
		perr.append(vals[4])

	elif(year>=endYear):	
		vals=ff.CalcForecastMultiVar(rawdatapath, extdatapath, year, startYear, fvars, fmonth, pmonth=pmonth,
			region=region, anomObs=0, numYearsReq=numYearsReq, weight=weight, outWeights=outWeights, 
			icetype=iceType, hemStr=hemStr, siiVersion='v2.1')
		print 'Forecast year:', year, fvars, 'DT SIE forecast:', vals[1], 'DT SIE Obs :', vals[0],'SIE forecast:', vals[2] 
		
		extentObsDt.append(np.nan)
		extentPredDt.append(vals[1])
		extentPredAbs.append(vals[2])
		anoms.append(np.nan)
		perr.append(vals[4])


rmsF= ff.rms(anoms)
rmsP= ff.rms(extentObsDt)
errorFore = '%.2f' %sqrt(rmsF)
errorExt = '%.2f' %sqrt(rmsP)
skill = '%.2f' %(1 - (rmsF/rmsP))
if (endYear<2017):
	skill2 = '%.2f' %(1 - (ff.rms(anoms[-9:]))/(ff.rms(extentObsDt[-9:])))
elif (endYear>=2017):
	skill2 = '%.2f' %(1 - (ff.rms(anoms[-10:-1]))/(ff.rms(extentObsDt[-10:-1])))
#skill3 = '%.2f' %(1 - (rms(anoms[-16:-8]))/(rms(extentObsDt[-16:-8])))

print str(startYearPred)+'-'+str(endYear)+' skill:',skill, '2008-'+str(endYear)+' skill:', skill2

if (outSkill==1):
	savetxt(skilldatapath+'Skill_'+varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)
		+'R'+str(region)+str(startYearPred)+str(endYear)+'W'+str(weight)+iceType+'.txt', 
		array([skill, errorFore, errorExt, skill2])[None], 
		header='skill error_forr error_obs skill2', fmt='%s')

if (outLine==1):
	array(extentPredDt).dump(linedatapath+iceType+'PredDt'
		   +varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYearPred)+str(endYear)
		   +'W'+str(weight)+'.txt')
	array(extentObsDt).dump(linedatapath+iceType+'ObsDt'
		   +varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYearPred)+str(endYear)
		   +'W'+str(weight)+'.txt')
	array(extentPredAbs).dump(linedatapath+iceType+'PredAbs'
		   +varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYearPred)+str(endYear)
		   +'W'+str(weight)+'.txt')

	array(anoms).dump(linedatapath+iceType+'anoms'
		   +varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYearPred)+str(endYear)
		   +'W'+str(weight)+'.txt')
	array(perr).dump(linedatapath+iceType+'perr'
		   +varStrsOut+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYearPred)+str(endYear)
		   +'W'+str(weight)+'.txt')



if (region==0):
	# GET extent AND DETREND
	years, extent = ff.get_ice_extentN(extdatapath, pmonth, startYear, endYear, icetype=iceType, version='v2.1',  hemStr=hemStr)
else: 
	if (hemStr=='N'):
		poleStr='A'
	elif (hemStr=='S'):
		poleStr='AA'
# If generating a regional forecast
	extentT=loadtxt(extdatapath+'ice_'+iceType+'_M'+str(pmonth)+'R'+str(region)+'_19792016'+poleStr)
	#get years and extent for years preceeding the given forecast year
	extent=extentT[0:endYear-startYear+1]

	if (endYear==2017):
		years=np.arange(startYear, endYear, 1)
		
	else:
		years=np.arange(startYear, endYear+1, 1)


if (plotSkill==1):
	ff.plot_forecast(years, extent, yearsP, extentPredAbs, extentObsDt, extentPredDt, perr, errorFore, errorExt, skill2, varStrsOut, fmonth, pmonth, weight, figpath, hemStr, iceType, region)


