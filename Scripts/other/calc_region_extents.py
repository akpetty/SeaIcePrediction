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

dataPath = '/Users/aapetty/GitHub/SeaIcePrediction/Data/'
outPath='/Users/aapetty/GitHub/SeaIcePrediction/DataOutput/Extent/'
figPath='/Users/aapetty/GitHub/SeaIcePrediction/Figures/'
cdrdatapath='/Volumes/PETTY_E/DATA/ICE_CONC/CDR/monthly/'

def get_region_mask_v2(anc_data_path, m, xypts_return=0):
	""" Read in NSIDC Arctic Ocean mask and transofrm to given projection
	"""

	header = 0
	datatype='uint8'
	file_mask = anc_data_path+'RegionMask_NH_ps25km_304x448_2021.bin'

	#8 - Arctic Ocean
	#9 - Canadian Archipelago
	#10 - Gulf of St Lawrence
	#11 - Land

	fd = open(file_mask, 'rb')
	region_mask = np.fromfile(file=fd, dtype=datatype)
	region_mask = np.reshape(region_mask[header:], [448, 304])

	mask_latf = open(anc_data_path+'/psn25lats_v3.dat', 'rb')
	mask_lonf = open(anc_data_path+'/psn25lons_v3.dat', 'rb')
	lats_mask = np.reshape(np.fromfile(file=mask_latf, dtype='<i4')/100000., [448, 304])
	lons_mask = np.reshape(np.fromfile(file=mask_lonf, dtype='<i4')/100000., [448, 304])

	if (xypts_return==1):
		# switched lat and lon from basemap
		xpts, ypts = m(lons_mask, lats_mask)

		return region_mask, xpts, ypts, lons_mask, lats_mask
	else:
		return region_mask, lons_mask, lats_mask

def get_month_conc_cdr(fileT, variable='cdr_seaice_conc_monthly', hem='nh', lowerConc=True, maxConc=True, mask=True):

	"""
	Grab the CDR ice concentration

	The CDR variable only exists after July 1987 so use Bootstrap before then. 

	I need to try and preserve pole hole as masked not nan so later routines will linearly interpolate over it.

	"""

	f = Dataset(fileT, 'r')
	
	conc = (f.variables[variable][0])
	if np.size(conc[conc<1])<10000:
		print('Issue with grabbing valid NT CDR data')
		return

	f.close()

	print('max value of raw CDR', np.amin(conc), np.amax(conc), np.nanmin(conc), np.nanmax(conc))

	if maxConc:
		conc = np.where(conc>1.,np.nan, conc)

	if lowerConc:
		conc = np.where(conc<0.15,0, conc)

	if mask:
		conc = ma.masked_where(conc>1., conc)
	
	print('max value of filtered CDR', np.amin(conc), np.amax(conc), np.nanmin(conc), np.nanmax(conc))

	return conc

m = Basemap(projection='npstere',boundinglat=66,lon_0=0, resolution='l'  )


region_mask, lonsR, latsR = get_region_mask_v2(dataPath+'OTHER/', m, xypts_return=0)
areaF=reshape(fromfile(file=open(dataPath+'OTHER/psn25area_v3.dat', 'rb'), dtype='<i4')/1000., [448, 304])/1e6


start_year=1979
end_year=2021
month = 9 # 9 = Sept
poleStr='A'

mstr = '%02d' % (month)

if poleStr=='A':
	hem='nh'
else:
	hen='sh'

ice_ext_mean=[]
ice_area_mean=[]

regions=[9, 12]
region_str='Canadian_v2'

print ('Month:', month)
for year in range(start_year, end_year+1):
	print (year)

	if year > 2020:
		try:
			# Try final data first
			fileT=glob(cdrdatapath+'/seaice_conc_monthly_'+hem+'_'+str(year)+mstr+'*.nc')[0]
		except:
			fileT=glob(cdrdatapath+'nrt/seaice_conc_monthly_icdr_'+hem+'_'+str(year)+mstr+'*.nc')[0]
	else:
		fileT=glob(cdrdatapath+'/seaice_conc_monthly_'+hem+'_'+str(year)+mstr+'*.nc')[0]

	ice_conc = get_month_conc_cdr(fileT, hem=hem, lowerConc=True, maxConc=True, mask=False)

	mask = np.isin(region_mask, regions)
	ice_conc = where(mask, ice_conc, 0)
	#ice_conc = where((region_mask==regions[0])|(region_mask==regions[1]|(region_mask==regions[2]), ice_conc, 0)

	ice_area = ice_conc*areaF
	ice_ext = where((ice_conc >=0.15), 1, 0)*areaF

	ice_ext_mean.append(sum(ice_ext))
	ice_area_mean.append(sum(ice_area))

savetxt(outPath+'ice_extent_M'+str(month)
+'_'+str(start_year)+str(end_year)+region_str, ice_ext_mean)

savetxt(outPath+'ice_area_M'+str(month)
+'_'+str(start_year)+str(end_year)+region_str, ice_area_mean)


Years=np.arange(start_year, end_year+1, 1)

fig = figure(figsize=(5,2.2))
ax=gca()
im1 = plot(Years, ice_ext_mean, 'k')
ax.set_ylim([0, 1.5])
ylabel('Extent [M km2]')
savefig(figPath+'ice_extent_M'+str(month)+region_str+'.png')


fig = figure(figsize=(5,2.2))
im1 = plot(Years, ice_area_mean, 'k')
ylabel('Area [M km2]')
savefig(figPath+'ice_area_M'+str(month)+region_str+'.png')


