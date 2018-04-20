############################################################## 
# Date: 01/01/17
# Name: grid_iceconc.py
# Author: Alek Petty
# Description: Grid the ice concentration data onto a coarser polar stereographic grid 

import matplotlib
matplotlib.use("AGG")
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from pylab import *
from scipy.io import netcdf
import numpy.ma as ma
from matplotlib import rc
from glob import glob
import matplotlib.patches as patches
from scipy.interpolate import griddata
import pandas as pd
from scipy import stats
import forecast_funcs as ff

from netCDF4 import Dataset

def plot_conc(xpts, ypts, conc_year, year, month, poleStr='A'):
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
	savefig(figpath+'/conc'+str(year)+str(month)+grid_str+poleStr+'v2.png', dpi=300)
	close(fig)
	
m = Basemap(projection='spstere',boundinglat=-52,lon_0=180, resolution='l'  )

datapath='../../Data/'
#datapath = '../../../../DATA/'
dataoutpath='../../DataOutput/IceConcAA/'
figpath='../../Figures/Antarctic/IceConc/'

start_year=2016
end_year=2016

startMonth=11 #3=April, 7=August
endMonth=11

poleStr='AA'

dx_res = 100000.
nx = int((m.xmax-m.xmin)/dx_res)+1; ny = int((m.ymax-m.ymin)/dx_res)+1
grid_str=str(int(dx_res/1000))+'km'
lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)
xptsG.dump(dataoutpath+'xpts'+grid_str+poleStr)
yptsG.dump(dataoutpath+'ypts'+grid_str+poleStr)

lats, lons = ff.get_psslatslons(datapath)
xpts, ypts =m(lons, lats)


alg=0 #0=Nasa team

#SSMIS 89.18 -	January 2008 to present
#SSM/I Pole Hole Mask	87.2	July 1987 through December 2007
#SMMR Pole Hole Mask - 84.5	November 1978 through June 1987

for month in xrange(startMonth, endMonth+1):
	print 'Month:', month
	for year in xrange(start_year, end_year+1):
		print year
		if (year>2015):
			
			ice_conc = ff.get_month_concSN_NRT(datapath, year, month, alg=alg, pole=poleStr, monthMean=1)
			

		else:
			ice_conc = ff.get_month_concSN(datapath, year, month, alg=alg, pole=poleStr)

		ice_conc=ma.masked_where(ice_conc<=0.15, ice_conc)
					
		ice_conc = ice_conc.filled(0)
		

		
		ice_concG = griddata((xpts.flatten(), ypts.flatten()),ice_conc.flatten(), (xptsG, yptsG), method='linear')
		ice_conc_ma=ma.masked_where(np.isnan(ice_concG), ice_concG)
		#ice_conc_ma=ma.masked_where((latsG>pmask), ice_conc_ma)

		plot_conc(xptsG, yptsG, ice_conc_ma, year, month, poleStr='A')

		ice_conc_ma.dump(dataoutpath+'ice_conc'+grid_str+str(month)+str(year)+poleStr)



