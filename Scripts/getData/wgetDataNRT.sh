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

#month=1
#monthStr=$(printf "%02d" $month)
#year=2018
#yearStr=$(printf "%04d" $year)
#sensor='f18'
#hemisphere='s'

echo 'Grabbing  data'

echo $1
echo $2
echo $3
echo $4

monthStr=$1
yearStr=$2
sensor=$3
hemisphere=$4

if [ "$hemisphere" = "n" ]
then
    folder1='north'
    folder2='ARCTIC'
else
	folder1='south'
    folder2='ANTARCTIC'
fi

echo $folder

# merge files along time dimension
for d in {01..32}
do
dayStr=$(printf "%02d" $d)	

FILE=https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0081_nrt_nasateam_seaice/"$folder1"/nt_"$yearStr$monthStr$dayStr"_"$sensor"_nrt_"$hemisphere".bin
echo $FILE
wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies \
--keep-session-cookies --no-check-certificate --auth-no-challenge=on \
-nH --cut-dirs=5 --directory-prefix=/Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Data/Temp/  --reject "index.html*" \
-np -e robots=off -r $FILE

done

# COPY TO CORRECT FOLDER
cp -r /Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Data/Temp/ /Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Data/ICE_CONC/NASA_TEAM/"$folder2"/NRT/
# REMOVE FROM TEMP FOLDER
rm -r /Users/aapetty/GitRepos/GitHub/SeaIcePrediction/Data/Temp/*








#ncrcat $FOLDER"b.e11.B20TRC5CNBDRD.f09_g16."$foo".cice.h.hi_nh.192001-200512.nc" $FOLDER"b.e11.BRCP85C5CNBDRD.f09_g16."$foo".cice.h.hi_nh.208101-210012.nc" $FOLDER"mergedtime."$foo".nc"



#import wget
#wget('https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0081_nrt_nasateam_seaice/north/nt_20180414_f18_nrt_n.bin', '../../Data/Temp/file.bin')
