# Sea ice prediction

## Contact
Alek Petty - [alekpetty.com](http://www.alekpetty.com) - [@alekpetty](https://twitter.com/alekpetty)

Python scripts used to produce tseasonal forecasts of Arctic and Antarctic sea ice. Results from our Arctic sea ice forecasts analysis were published in JGR Earth's Future (see references below). The original github repo for those scripts can be found here: https://github.com/akpetty/ArcticSeaIcePrediction2017. This is an updated repository that features a few upgrades to make the forecast scripts easier to use, and more felxible (e.g. for forecasting Antarctic sea ice, different months of the year etc.).

More information to be added shortly.

## Getting Started

### Conda installation

I recommend using conda to install your Python environment. You can build this from the environment I used and exported:
```
conda env create -f environment.yml
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
+ ```Scripts/forecasts/PetysPredictionNotebook.ipynb```: A Jupyter notebook explaining how we produce a sea ice forecast.

You can run the Jupyter Notebooks by running the following command in the relevant directory 

```
jupyter notebook
``` 

### Python scripts

The Python scripts found in the Scripts directory by executing
```
python script.py
``` 

Individual descriptions should be included at the top of each script.

Add info about running and using Jupyter Notebooks.

Add info about how to generate forecasts using the runFOrecasts.py script.


### Data

I include the archived and near real-time ice concentration data in ```/Data/IceConc/```. These were obtained from the following, pubclically available, data repositories:

Sea ice concentration data (final): http://nsidc.org/data/nsidc-0051    
Sea ice concentration data (near real-time): https://nsidc.org/data/nsidc-0081   

In the paper below we also used data regarding melt onset and melt pond coverage:

Melt onset data: http://neptune.gsfc.nasa.gov/csb/index.php?section=54   
(note that the melt onset data are not made avilable each year near real-time, so contact me if required).

Simulated melt pond data were provided by CPOM-Reading, with the detrended forecast data included in this repo.

### References

Petty, A. A., D. Schroder, J. C. Stroeve, T. Markus, J. Miller, N. T. Kurtz, D. L. Feltham, D. Flocco (2017), Skillful spring forecasts of September Arctic sea-ice extent using passive microwave sea ice observations, Earth’s Future, 4 , doi:10.1002/2016EF000495.

Petty, A. A., J. C. Stroeve, P. R. Holland, L. N. Boisvert, A. C. Bliss, N. Kimura, W. N. Meier (2018), The Arctic sea ice cover of 2016: A year of record-low highs and higher-than-expected lows, The Cryosphere, 12, 433-452, doi:10.5194/tc-12-433-2018.

More informwation can be found on my website www.alekpetty.com

Contact me if you any any questions!

Alek
