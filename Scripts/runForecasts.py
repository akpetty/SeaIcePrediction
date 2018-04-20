
import subprocess
import shlex

month=2 # here month starts from 1
year=2018
sensor=f17
subprocess.call('./getData/wgetData.sh month year sensor %s %s %s' % (month, year, sensor) )


import grid_iceconcA 

grid_iceconcA.main(year, month-1)





