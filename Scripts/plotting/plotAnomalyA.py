############################################################## 
# Date: 01/01/17
# Name: plot_iceconcAA.py
# Author: Alek Petty

import matplotlib
matplotlib.use("AGG")

from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
from pylab import *
from matplotlib import rc

rcParams['xtick.major.size'] = 2
rcParams['ytick.major.size'] = 2
rcParams['axes.linewidth'] = .3
rcParams['lines.linewidth'] = .3
rcParams['patch.linewidth'] = .3

rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize']=8
rcParams['ytick.labelsize']=8
rcParams['legend.fontsize']=8
rcParams['font.size']=8
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

m = Basemap(projection='npstere',boundinglat=65,lon_0=0, resolution='l'  )


figpath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/Figures/forecasts/Arctic/Anomalies/'
dataoutpath='/Users/aapetty/GitHub/akpetty/SeaIcePrediction/DataOutput/'

startYear=1979
yearT=2023
fmonth=6 #6=june
pmonth=9 #6=june


hemStr='N'
if (hemStr=='S'):
	weightsdataoutpath=dataoutpath+'forecasts/Antarctic/'
	poleStr='AA'
elif (hemStr=='N'):
	weightsdataoutpath=dataoutpath+'forecasts/Arctic/'
	poleStr='A'

#xpts, ypts, Melt_onset_years =pfuncs.get_meltonset(m, data_path, 'melt_onset', -0.5, 0, 0, str(start_year), str(end_year))

xpts=load(dataoutpath+'IceConcA/CDR/xpts100km'+poleStr, allow_pickle=True)
ypts=load(dataoutpath+'IceConcA/CDR/ypts100km'+poleStr, allow_pickle=True)

varStr='conc'
icetype='extent'
region=0

if ((region=='A')):
	regionOut='Alaska'
	weightsdataoutpath=dataoutpath+'/Alaska/'
else:

	regionOut='Arctic'


r_valsDT=load(weightsdataoutpath+'rvalsDT'+varStr+icetype+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(yearT)+'.txt', allow_pickle=True)
predvarDT=load(weightsdataoutpath+'predvarYrDT'+varStr+icetype+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(yearT)+'.txt', allow_pickle=True)


minval=0.5
maxval=-0.5


fig = figure(figsize=(3,3*0.8))
ax=gca()

im1 = m.contourf(xpts , ypts, r_valsDT,levels=[0.25, 1.], colors='none', hatches=['///'], zorder=4)
	
im2 = m.pcolormesh(xpts, ypts, predvarDT,vmin=minval, vmax=maxval, cmap=cm.RdBu_r, shading='interp', zorder=2)
#im1 = m.pcolormesh(xpts , ypts, rvals, cmap=cm.cubehelix, vmin=minval, vmax=maxval,shading='flat', zorder=2)

#m.fillcontinents(color='w',lake_color='0.9', zorder=2)
m.drawcoastlines(linewidth=0.25, zorder=5)
m.drawparallels(np.arange(90,-90,-10), linewidth = 0.25, zorder=3)
m.drawmeridians(np.arange(-180.,180.,30.), linewidth = 0.25, zorder=3)

label_str=r'$\Delta$ A'
cax = fig.add_axes([0.81, 0.05, 0.03, 0.45])
cbar = colorbar(im2,cax=cax, orientation='vertical', extend='both', use_gridspec=True)
cbar.set_label(label_str, labelpad=4, rotation=0)
cbar.set_ticks(np.linspace(minval, maxval, 3))
cbar.solids.set_rasterized(True)

#ax.annotate(varStr, xy=(0.5, 1.01),xycoords='axes fraction', horizontalalignment='center', verticalalignment='bottom', zorder=10)
ax.annotate(regionOut+'\n'+str(yearT)+'\nFM:'+str(fmonth)+' PM:'+str(pmonth)+'\n\nFvar:'+varStr+'\nPvar:'+icetype, xy=(1.005, 0.98), xycoords='axes fraction', horizontalalignment='left', verticalalignment='top', rotation=0, zorder=10)

subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.8)
savefig(figpath+'anomaly'+varStr+icetype+'fm'+str(fmonth)+'pm'+str(pmonth)+'R'+str(region)+str(startYear)+str(yearT)+poleStr+'.png', dpi=300)
close(fig)




