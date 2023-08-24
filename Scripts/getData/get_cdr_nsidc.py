#!/usr/bin/env python
# NSIDC
# Sample script to download all the files within one directory on the FTP server
#
# Requires Python3 and the ftplib and os libraries

from ftplib import FTP
import os

### The following 3 variables can be changed ###
# 1. Set the directory you would like to download the files to
destdir='path to directory'

# 2. Set the path to the FTP directory that contains the data you wish to download.
# This example is for the daily northern hemisphere data from the Sea Ice Index
# https://nsidc.org/data/g02135
directory = '/DATASETS/NOAA/G02202/north/monthly'

# 3. Set the password which will be your email address
password = 'alekpetty@gmail.com'

############################################
### Don't need to change this code below ###
############################################
# FTP server
ftpdir = 'sidads.colorado.edu'

#Connect and log in to the FTP
print('Logging in')
ftp = FTP(ftpdir)
ftp.login('anonymous',password)

# Change to the directory where the files are on the FTP
print('Changing to '+ directory)
ftp.cwd(directory)

# Get a list of the files in the FTP directory
files = ftp.nlst()
files = files[2:]
print(files)

#Change to the destination directory on own computer where you want to save the files
os.chdir(destdir)

#Download all the files within the FTP directory
for file in files:
    print('Downloading...' + file)
    ftp.retrbinary('RETR ' + file, open(file, 'wb').write)

#Close the FTP connection
ftp.quit()