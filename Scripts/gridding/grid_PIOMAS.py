############################################################## 
# Date: 01/01/16
# Name: grid_PIOMAS.py
# Author: Alek Petty
# Description: Grid the PIOMAS data onto the forecast grid

import matplotlib
matplotlib.use("AGG")
import pred_funcs as pfuncs
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
import pred_funcs as pfuncs

m = Basemap(projection='npstere',boundinglat=65,lon_0=0, resolution='l'  )

rawdatapath='../../../DATA/'
pmasdatapath=rawdatapath+'/PIOMAS/heff_txt/'
dataoutpath='./Data_output/PMAS_OUT/'
meltoutpath='./Data_output/MELT_OUT/'
figpath='./Figures/'

grid_str='100km'

xpts100 =load(meltoutpath+'xpts'+grid_str)
ypts100 =load(meltoutpath+'ypts'+grid_str)

xpts100.dump(dataoutpath+'xpts'+grid_str)
ypts100.dump(dataoutpath+'ypts'+grid_str)

start_year=1979
end_year=2015
month=4 #June

for year in xrange(start_year, end_year+1, 1):
	print year
	xptsP, yptsP, thickness=pfuncs.get_pmas_month(m, rawdatapath, year,month=month)

	thickness_year = griddata((xptsP, yptsP),thickness, (xpts100, ypts100), method='linear')
	#thickness_year_ma=ma.masked_where(thickness_year<0.01, thickness_year)

	thickness_year.dump(dataoutpath+'pmas'+grid_str+str(year)+str(month))





