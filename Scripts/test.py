from glob import glob
from pylab import *

alg=0
pole='A'
datapath='../Data/'
year=2020
month=5

if (alg==0):
	team = 'NASA_TEAM'
	team_s = 'nt'
	header = 300
	datatype='uint8'
	scale_factor=250.

if (pole=='A'):
	poleStr='ARCTIC'
	rows=448
	cols=304


month_str = '%02d' % (month+1)
year_str=str(year)
files = glob(datapath+'/ICE_CONC/'+team+'/'+poleStr+'/NRT/*'+str(year)+month_str+'*')

fd = open(files[0], 'r')
data = fromfile(file=fd, dtype=datatype)
data = data[header:]

ice_conc[x] = reshape(data, [rows, cols])
		