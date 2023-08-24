############################################################## 
# Date: 10/05/16
# Name: calc_ice_extent_regionsAA.py
# Author: Alek Petty
# Description: Script to calculate Antarctic ice exent for various regions.

import matplotlib
matplotlib.use("AGG")
import sys
sys.path.append("../")

import numpy as np
from pylab import *
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap, shiftgrid
from glob import glob
from netCDF4 import Dataset
import forecast_funcs as ff

dataPath = '/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Data/'
outPath='/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/DataOutput/Extent/'
figPath='/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Figures/'

m = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )


lats, lons = ff.get_psnlatslons(dataPath)
areaF=reshape(fromfile(file=open(dataPath+'/OTHER/psn25area_v3.dat', 'rb'), dtype='<i4')/1000., [448, 304])/1e6
region_mask, xpts, ypts = ff.get_region_mask_sect(dataPath, m, xypts_return=1)



start_year=1979
end_year=2018
month = 9 # 9 = Sept
alg=0 #0=NASA TEAM, 1=BOOTSTRAP

ice_ext_mean=[]
ice_area_mean=[]

regions=[12, 13]


print ('Month:', month)
for year in range(start_year, end_year+1):
	print (year)

	if (year>2016):
		ice_conc = ff.get_month_concSN_NRT(dataPath, year, month-1, alg=alg, pole='A', lowerConc=1, maxConc=1, mask=1, monthMean=1)
		#ice_conc=ma.masked_where(ice_conc<=0.15, ice_conc)
	else:
		ice_conc = ff.get_month_concSN(dataPath, year, month-1, alg=alg, pole='A', lowerConc=1, maxConc=1, mask=1)
				
	ice_conc = ice_conc.filled(0)


	#if (year>2015):
		# Subtract 1 to get month index starting from zero.
	#	ice_conc = ff.get_month_concSN_NRT(dataPath, year, month-1, alg=alg, pole='A')
	#else:
	#	ice_conc = ff.get_month_concSN(dataPath, year, month-1, alg=alg, pole='A')
				
	#ice_conc = ice_conc.filled(0)
	
	#ice_conc = where((ice_conc <=0.15), 0, ice_conc)

	ice_conc = where((region_mask==regions[0])|(region_mask==regions[1])
		, ice_conc, 0)

	ice_area = ice_conc*areaF
	ice_ext = where((ice_conc >=0.15), 1, 0)*areaF
	

	ice_ext_mean.append(sum(ice_ext))
	ice_area_mean.append(sum(ice_area))

savetxt(outPath+'ice_extent_M'+str(month)+'RA'
+'_'+str(start_year)+str(end_year)+'A', ice_ext_mean)

savetxt(outPath+'ice_area_M'+str(month)+'RA'
+'_'+str(start_year)+str(end_year)+'A', ice_area_mean)


Years=np.arange(start_year, end_year+1, 1)

fig = figure(figsize=(5,2.2))
im1 = plot(Years, ice_ext_mean, 'k')
ylabel('Extent [M km2]')
savefig(figPath+'ice_extent_M'+str(month)+'Alaska.png')


fig = figure(figsize=(5,2.2))
im1 = plot(Years, ice_area_mean, 'k')
ylabel('Area [M km2]')
savefig(figPath+'ice_area_M'+str(month)+'Alaska.png')


