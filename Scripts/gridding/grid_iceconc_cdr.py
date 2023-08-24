"""  grid_iceconcA.py
	Script to grid sea ice concentration data onto our given domain
	Run with e.g. python grid_iceconcA.py 
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
#from scipy.io import netcdf
import numpy.ma as ma
from matplotlib import rc
from glob import glob
import matplotlib.patches as patches
from scipy.interpolate import griddata
import pandas as pd
from netCDF4 import Dataset
#from scipy import stats
import sys
sys.path.append('../')
import os

import forecast_funcs as ff



def plot_conc(figpath, m , xpts, ypts, conc_year, year, month, grid_str, outStr='A'):
	textwidth=4.
	fig = figure(figsize=(textwidth,textwidth))
	subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99)

	#ax1=subplot(1, 3, 1)
	minval=0
	maxval=1
	#ADD GRIDSIZE=NUMBER KWARG TO HEXBIN IF YOU WANT TO CHANGE SIZE OF THE BINS
	im1 = m.pcolormesh(xpts , ypts, conc_year, cmap=cm.Blues_r, vmin=minval, vmax=maxval,shading='flat', zorder=2)
	#im2 = m.contour(xpts , ypts, ma.mean(Pressure, axis=0),levels=[990, 1000, 1100],colors='k', zorder=4)
	m.drawcoastlines(linewidth=0.5, zorder=5)
	m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=3)
	m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=3)
	#m.plot(xptsR, yptsR, '--', linewidth = 2, color='k', zorder=5)

	#ADD COLORBAR TO MAP
	#bbox_args = dict(fc="white")
	#ax1.annotate('.                   \n             \n        ', xy=(0.02, 0.98), bbox=bbox_args,xycoords='axes fraction', horizontalalignment='left', verticalalignment='top', zorder=10)

	#ax1.annotate(files[x][-8:-4]+'-'+files[x][-4:-2]+'-'+files[x][-2:], xy=(0.98, 0.98), bbox=bbox_args,xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', zorder=10)
	label_str='Conc'
	#ax1.annotate('AIRS temp anomaly from 2003-2014 mean', xy=(0.02, 0.02), bbox=bbox_args,xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', zorder=10)
	cax = fig.add_axes([0.04, 0.9, 0.2, 0.035])
	cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
	cbar.set_label(label_str, labelpad=1)
	cbar.set_ticks(np.arange(minval, maxval+1, 1))
	cbar.solids.set_rasterized(True)
	#SHIFT COLOR SPACE SO OFF WHITE COLOR IS AT 0 m
	cbar.set_clim(minval-0.5, maxval)
	savefig(figpath+'/conc'+str(year)+str(month)+grid_str+outStr+'.png', dpi=300)
	close(fig)

def get_month_conc_cdr(fileT, variable='nsidc_nt_seaice_conc_monthly', hem='nh', lowerConc=True, maxConc=True, mask=True):

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

def main(year, month, poleStr='A', outputGrid=1):
	
	datapath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/Data/'
	#cdrdatapath='/Volumes/PETTY_E/DATA/ICE_CONC/CDR/monthly/'
	cdrdatapath='/Users/aapetty/DATA/IceConc/CDR/monthly/'

	if poleStr=='A':
		hem='nh'
		hem_long='north'
		lats, lons = ff.get_psnlatslons(datapath)
		m = Basemap(projection='npstere',boundinglat=65,lon_0=0, resolution='l'  )
		dataoutpath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/DataOutput/IceConcA/CDR/'
		figpath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/Figures/Arctic/IceConc/CDR/'
	
	else:
		hem='sh'
		hem_long='south'
		lats, lons = ff.get_psslatslons(datapath)
		m = Basemap(projection='spstere',boundinglat=-55,lon_0=180, resolution='l'  )
		dataoutpath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/DataOutput/IceConcAA/CDR/'
		figpath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/Figures/Antarctic/IceConc/CDR/'

	if not os.path.exists(datapath):
		os.makedirs(datapath)

	#if not os.path.exists(cdrdatapath):
	#	os.makedirs(cdrdatapath)

	if not os.path.exists(dataoutpath):
		os.makedirs(dataoutpath)

	if not os.path.exists(figpath):
		os.makedirs(figpath)

	dx_res = 100000.
	nx = int((m.xmax-m.xmin)/dx_res)+1; ny = int((m.ymax-m.ymin)/dx_res)+1
	grid_str=str(int(dx_res/1000))+'km'
	lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)	

	print(poleStr)

	if (outputGrid==1):
		xptsG.dump(dataoutpath+'xpts'+grid_str+poleStr)
		yptsG.dump(dataoutpath+'ypts'+grid_str+poleStr)

	mstr = '%02d' % (month+1)
	
	xpts, ypts =m(lons, lats)
	#f = Dataset(datapath+'/OTHER/NIC_valid_ice_mask.N25km.01.1972-2007.nc', 'r')
	#ice_flag = f.variables['valid_ice_flag'][:]
	#region_mask = ff.get_region_mask_sect(datapath, m, xypts_return=0)

	try:
		# Try final data first

		fileT=glob(cdrdatapath+'v4/'+hem_long+'/seaice_conc_monthly_'+hem+'_'+str(year)+mstr+'*.nc')[0]
	except:
		print('No final data so getting NRT')
		fileT=glob(cdrdatapath+'nrt/'+hem_long+'/seaice_conc_monthly_icdr_'+hem+'_'+str(year)+mstr+'*.nc')[0]
	
	ice_conc = get_month_conc_cdr(fileT, hem=hem, lowerConc=True, maxConc=True, mask=False)

	print('Ice conc shape:', ice_conc.shape)		
	#ice_conc = ice_conc.filled(0)
	
	#ice_conc = ma.masked_where(ice_conc>1., ice_conc)
	#ice_conc = ma.where(ice_conc>1.,0, ice_conc)
	#ice_conc = ma.where(ice_conc<0.15,0, ice_conc)
	#ice_conc = np.where((ice_flag >=1.5), 0, ice_conc)

	# get mean conc around the pole hole (time varying) and apply to the pole hole. 
	# do i need to do this with CDR? 
	#pmask=ff.get_pmask(year, month)
	#concHole=ma.mean(ice_conc[(lats>pmask-0.5) & (lats<pmask)])
	#ice_conc = where((lats >=pmask-0.5), concHole, ice_conc)
	
	#ice_conc[where(region_mask>18)]=0

	ice_concG = griddata((xpts.flatten(), ypts.flatten()),ice_conc.flatten(), (xptsG, yptsG), method='linear')
	ice_conc_ma=ma.masked_where(np.isnan(ice_concG), ice_concG)

	print(ice_conc)
	print(ice_concG)
	print(ice_conc_ma)

	print('max value of gridded/masked CDR data', np.amin(ice_conc_ma), np.amax(ice_conc_ma), np.nanmin(ice_conc_ma), np.nanmax(ice_conc_ma))

	plot_conc(figpath, m, xptsG, yptsG, ice_conc_ma, year, month, grid_str, outStr=hem+'cdr_nt')

	ice_conc_ma.dump(dataoutpath+'ice_conc'+grid_str+str(month)+str(year)+poleStr+'cdr_nt')





startYear=2023
endYear=2023

startMonth=5 #3=April, 7=August

endMonth=5
poleStr='A'
#-- run main program
if __name__ == '__main__':
	for y in range(startYear, endYear+1, 1):
		for m in range(startMonth, endMonth+1):
			print (y, m)
			main(y, m, poleStr=poleStr)





