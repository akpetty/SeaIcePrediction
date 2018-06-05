"""  runForecasts.py
	Script to run the various forecast scripts
	Run with e.g. python runForecasts.py 
	Author:
		Alek Petty

	Update history:
		04/20/2018: Version 1
"""

import subprocess
import shlex

month=5 # here month starts from 1
monthStr='%02d' %(month)
year=2017
yearStr=str(year)
sensor='f18'
hemisphere='s'



# ADD A CHECK TO SEE IF DATA EXISTS FIRST!

subprocess.check_call(['./wgetDataNRT.sh', monthStr, yearStr, sensor, hemisphere])