"""  getData.py
	Script to download ice conc data
	Run with e.g. python runForecasts.py 
	Author:
		Alek Petty

	Update history:
		04/20/2018: Version 1
"""

import subprocess
import shlex

month=9 # here month starts from 1
monthStr='%02d' %(month)
year=2018
yearStr=str(year)
sensor='f18'
hemisphere='n'



# ADD A CHECK TO SEE IF DATA EXISTS FIRST!

subprocess.check_call(['./wgetDataNRT.sh', monthStr, yearStr, sensor, hemisphere])