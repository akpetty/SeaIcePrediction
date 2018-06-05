# Sea ice prediction

## Contact
Alek Petty - [alekpetty.com](http://www.alekpetty.com) - [@alekpetty](https://twitter.com/alekpetty)

## Introduction

Included in this repository is Python code to produce seasonal forecasts of Arctic and Antarctic sea ice. Results from our summer Arctic sea ice forecasts were published in JGR Earth's Future (see references below). In that paper we demonstrated 'skillful' seasonal forecasts of summer (September) Arctic sea ice extent. We define 'skillful' as a forecast that performs better than simply continuing the linear, downward, trend in Arctic sea ice extent.

The model can forecast the monthly mean sea ice area or extent across either hemisphere for any month or year (since 1980) using gridded spatial data of the sea ice state (primairly sea ice concetration, but also melt onset and melt pond coverage). The underlying forecast model utilizes a simple linear regression framework, detrended gridded input data, and a simple correlation/weighting scheme. Weighting the de-trended input data was shown to be an important process for improving forecast skill for longer (several month) lead times.

The original github repo for those scripts can be found here: https://github.com/akpetty/ArcticSeaIcePrediction2017. This is an updated repository that features a few upgrades to make the forecast scripts easier to use, and more flexible (e.g. for forecasting Antarctic sea ice, different months of the year etc.).

I have included some information below that should get you up and running. I recently ported the scripts over to work in Python 3 (3.6) and all seems to be the same as it was in Python 2.7.


### Desired model improvements.

* Investigate more sophisticated prediction models. This model currently uses a simple linear regression framwework at its core. More sophisticated machine learning methods are likely to offer improvements in skill, especially at predicting extreme deviations from the long-term trend, which linear regression models struggle with. The Jupyter Notebook should provide an easy framework to test different models for a given forecast.

* Explore forecasts using other input datasets. I have tested thickness data from the ice-ocean model PIOMAS (not shown) which did not show particularly improved forecast skill compared to using ice concentration, however this was only for the forecasts of summer sea ice extent. I'm also keen to explore sea surface temperature data (e.g. NOAA's OISST data) as a further input variable.  

* Improve the near real-time scripting to more easily generate a forecast and produce a plot of a forecast for a given month/year.


## Getting Started

### Conda installation

I recommend using conda to install your Python environment. You can build this from the environment I used and exported:
```
conda env create -n py36 -f environment.yml
```

Alternatively you can create your own Python environment with the following command:
```
conda create -n py36 python=3.6 anaconda
```
and the following packages:
```
conda install matplotlib
conda install basemap
conda install pandas
conda install statsmodels
conda install netCDF4
conda install wget
conda install imageio
``` 

The conda environment can be activated with ```source activate py36``` where py36 is the environment name.

To make sure the Jupyter Notebooks work, you may need to instal nb_conda as there can be issues with choosing the appropriate kernal

```
conda install nb_conda
``` 

### Jupyter notebooks

I have generated some Jupyter Notebooks to provide a walk through of the gridding and forecast methodology. 

+ ```Scripts/gridding/GriddingDemo.ipynb```: A Jupyter notebook explaining how we grid sea ice concentration data to our own grid.
+ ```Scripts/forecasts/PetysPredictionNotebook.ipynb```: A Jupyter notebook explaining the sea ice forecast methodology.

You can run the Jupyter Notebooks by running the following command in the relevant directory 

```
jupyter notebook
``` 

### Python scripts

The Python scripts found in the ```/Scripts/``` directory by executing
```
python script.py
``` 

Individual descriptions should be included at the top of each script.

The script ```/Scripts/runForecasts.py``` includes commands to run the various scripts needed to go from raw ice concentration data to a given sea ice forecast. This is a work in progress. Alternatively you can just run the individual scripts needed as follows:

* The shell script ```./getData/wgetData.sh month year sensor``` grabs sea ice concentration data from the NSIDC and stores it in the ```/Data``` folder.

* The Python script ```/Scripts/gridding/grid_iceconcA.py``` grids this raw sea ice concentration data to a consistent (100 km) polar stereographic grid

* The Python script ```/Scripts/forecasts/calcForecastYears.py``` uses the gridded data to generate a sea ice forecast for a given year/month etc. This also includes a main loop to run forecasts over a number of years and months.

* The Python script ```/Scripts/plotting/plotForecasts.py``` produces a plot of the sea ice forecast for the given year/month etc.

### Data

I have included the archived and near real-time ice concentration data in ```/Data/IceConc/```. These were obtained from the following, publically available, data repositories:

Sea ice concentration data (final): http://nsidc.org/data/nsidc-0051    
Sea ice concentration data (near real-time): https://nsidc.org/data/nsidc-0081   

In the paper below we also used data regarding melt onset and melt pond coverage:

Melt onset data: http://neptune.gsfc.nasa.gov/csb/index.php?section=54   
(note that the melt onset data are not made avilable each year near real-time, so contact me if required).

Simulated melt pond data were provided by CPOM-Reading, with the detrended forecast data included in this repo.


### References

Petty, A. A., D. Schroder, J. C. Stroeve, T. Markus, J. Miller, N. T. Kurtz, D. L. Feltham, D. Flocco (2017), Skillful spring forecasts of September Arctic sea-ice extent using passive microwave sea ice observations, Earth’s Future, 4 , doi:10.1002/2016EF000495.

Petty, A. A., J. C. Stroeve, P. R. Holland, L. N. Boisvert, A. C. Bliss, N. Kimura, W. N. Meier (2018), The Arctic sea ice cover of 2016: A year of record-low highs and higher-than-expected lows, The Cryosphere, 12, 433-452, doi:10.5194/tc-12-433-2018.

More information can be found on my website www.alekpetty.com

Contact me if you any any questions!

Alek
