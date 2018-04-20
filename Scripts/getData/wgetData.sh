#!/bin/bash
# https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
# https://nsidc.org/support/faq/what-options-are-available-bulk-downloading-data-https-earthdata-login-enabled
#> cd ~
#> touch .netrc
#> echo "machine urs.earthdata.nasa.gov login uid_goes_here password password_goes_here" > .netrc
#> chmod 0600 .netrc


# MONTHLY FINAL NASA TEAM DATA
#wget -r --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -nH --cut-dirs=3 --directory-prefix=../../Data/Temp/  --reject "index.html*" -np -e robots=off -r https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0051_gsfc_nasateam_seaice/final-gsfc/north/monthly/



# DAILY NRT NASA TEAM DATA
#wget -r --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies 
#--keep-session-cookies --no-check-certificate --auth-no-challenge=on 
#-nH --cut-dirs=3 --directory-prefix=../../Data/Temp/  --reject "index.html*" 
#-np -e robots=off -r https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0081_nrt_nasateam_seaice/north/

month=3
monthStr=$(printf "%02d" $month)
year=2018
yearStr=$(printf "%04d" $year)
sensor='f18'
echo $monthStr
echo $yearStr

# merge files along time dimension
for d in {01..32}
do
dayStr=$(printf "%02d" $d)	

FILE=https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0081_nrt_nasateam_seaice/north/nt_"$yearStr$monthStr$dayStr"_"$sensor"_nrt_n.bin
echo $FILE
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
--keep-session-cookies --no-check-certificate --auth-no-challenge=on \
-nH --cut-dirs=5 --directory-prefix=../../Data/Temp/  --reject "index.html*" \
-np -e robots=off -r $FILE

done

cp -r ../../Data/Temp/* ../../Data/ICE_CONC/NASA_TEAM/ARCTIC/NRT/









#ncrcat $FOLDER"b.e11.B20TRC5CNBDRD.f09_g16."$foo".cice.h.hi_nh.192001-200512.nc" $FOLDER"b.e11.BRCP85C5CNBDRD.f09_g16."$foo".cice.h.hi_nh.208101-210012.nc" $FOLDER"mergedtime."$foo".nc"



#import wget
#wget('https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0081_nrt_nasateam_seaice/north/nt_20180414_f18_nrt_n.bin', '../../Data/Temp/file.bin')
