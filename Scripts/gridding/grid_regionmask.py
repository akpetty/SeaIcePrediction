""" grid_regionmask.py
	Script to grid and plot the region mask
	Run with e.g. python grid_regionmask.py 
	Author:
		Alek Petty

	Update history:
		06/28/2018: Version 1
"""
import matplotlib
matplotlib.use("AGG")
import sys
from pylab import *
sys.path.append('../')
import forecast_funcs as ff
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata


dataPath = '../../Data/'
dataOutPath= '../../DataOutput/Regions/'
figPath='../../Figures/'

# BASEMAP INSTANCE
m = Basemap(projection='npstere',boundinglat=65,lon_0=0, resolution='l'  )
dx_res = 100000.
nx = int((m.xmax-m.xmin)/dx_res)+1; ny = int((m.ymax-m.ymin)/dx_res)+1
grid_str=str(int(dx_res/1000))+'km'
lonsG, latsG, xptsG, yptsG = m.makegrid(nx, ny, returnxy=True)


region_mask, xpts, ypts = ff.get_region_mask_sect(dataPath, m, xypts_return=1)
region_mask=ma.masked_where(region_mask==1.5, region_mask)
region_mask=ma.masked_where(region_mask>19.5, region_mask)

# Grid data to 100 km using nearest neighbor
region_maskG = griddata((xpts.flatten(), ypts.flatten()),region_mask.flatten(), (xptsG, yptsG), method='nearest')
region_maskG.dump(dataOutPath+'regionMaskA'+grid_str)


rcParams['axes.labelsize'] = 9
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9
rcParams['font.size']=9
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
textwidth=4.
fig = figure(figsize=(textwidth,textwidth))
subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99)

#ax1=subplot(1, 3, 1)
minval=0
maxval=20
#ADD GRIDSIZE=NUMBER KWARG TO HEXBIN IF YOU WANT TO CHANGE SIZE OF THE BINS
im1 = m.pcolormesh(xpts , ypts, region_mask, cmap=cm.viridis, vmin=minval, vmax=maxval,shading='gouraud', zorder=2)
#im2 = m.contour(xpts , ypts, ma.mean(Pressure, axis=0),levels=[990, 1000, 1100],colors='k', zorder=4)
m.drawcoastlines(linewidth=0.5, zorder=5)
m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=3)
m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=3)

#ax1.annotate(files[x][-8:-4]+'-'+files[x][-4:-2]+'-'+files[x][-2:], xy=(0.98, 0.98), bbox=bbox_args,xycoords='axes fraction', horizontalalignment='right', verticalalignment='top', zorder=10)
label_str='Region'
#ax1.annotate('AIRS temp anomaly from 2003-2014 mean', xy=(0.02, 0.02), bbox=bbox_args,xycoords='axes fraction', horizontalalignment='left', verticalalignment='bottom', zorder=10)
cax = fig.add_axes([0.02, 0.88, 0.25, 0.035])
cbar = colorbar(im1,cax=cax, orientation='horizontal', extend='both', use_gridspec=True)
cbar.set_label(label_str, labelpad=1)
cbar.set_ticks(np.arange(minval, maxval+1, 4))
cbar.solids.set_rasterized(True)
#SHIFT COLOR SPACE SO OFF WHITE COLOR IS AT 0 m
#cbar.set_clim(minval, maxval)
savefig(figPath+'/regionsArctic.png', dpi=300)
#savefig(figpath+'/regions_'+str(region)+'.png', dpi=300)
close(fig)