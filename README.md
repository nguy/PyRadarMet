PyRadarMet
=============== 
Python Fundamental Calculations in Radar Meteorology package notes

Originally Created:   5 February 2014

## Author
Nick Guy - nick.guy@uwyo.edu

Special thanks to Timothy Lang and Kai Muehlbauer for the insights and contributions.

## Package Details
attenuation.py – Routines to calculate coefficients useful in attenuation calculations

conversion.py – Routines to convert linear to log reflectivity and vice versa

doppler.py – Routines to calculate a number of fundamental Doppler radar characteristics
              including unambiguous range and velocity, 
              “Doppler dilemma” equation, dual PRF Vmax,…

geometry.py – Routines to calculate such characteristics as effective radius, 
               half-power radius, ray height, sample volumes, range corrections, 
               beam blockage fractions,…

system.py – Routines to calculate such characteristics as wavelength, frequency, 
             pulse length, radar constant, effective antenna area, thermal noise, …

variables.py – Routines to calculate such characteristics as CDR, LDR, ZDR, ZDP,…  
                This is really basic at the moment as I haven’t even attempted
                dual-pol calcs as of yet.
                
zdrcal.py – Routines to calculate ZDR offset of a dual-polarimetric radar.  This is called
              by an executable called cal_zdr.
              
BeamBlock - A class that allows the calculation of the geometric beam blocking.

## News
This module has be ported over to an [R package](http://cran.r-project.org/web/packages/radar/) by Jose Gama.

Also, PyRadarMet is undergoing a port into the [wradlib](http://wradlib.bitbucket.org/) python package.

Most fundamental modules are (or have been) ported currently into both packages.

Special packages (e.g. BeamBlock, ZDR cal scripts) are currently only available in PyRadarMet.

## Installation
There are two options for install.  The setup script now downloads the data listed below
automatically.  

To download and install the data within the package directory:

```python
python setup.py install
```

This is the most robust way, as the data is always a relative path stored within the python 'egg'.

To download the data and install it to a user-specified location an environmental variable must be set.
The variable is 'GTOPO_DATA' and can be set on unix systems (in a bash environment) as such:
```python
export GTOPO_DATA=/some/direcoty/for/the/data
```

It is recommended that this be done in a shell profile file like .bashrc or .bash_profile 
Then run the install command:
```python
python setup.py install
```

The install took about 150 seconds on a university connection when installed in an external directory.

The install took about 250 seconds when installed as part of the package.

##Data
The data folder holds Digital Elevation Model data used in the BeamBlock routine.
The data is the GTOPO30, a 30-arc second (~ 1 km) elevation data set.  
Data available from the U.S. Geological Survey.
Downloaded from the [long term archive](https://lta.cr.usgs.gov/),

 using the [EarthExplorer](http://earthexplorer.usgs.gov/) tool.

## Dependencies

Developed on the Anaconda distribution (1.9.1 & Python 2.7.7), tested to:
Anaconda 2.1.0 and Python 2.7.8
[Anaconda](https://store.continuum.io/cshop/anaconda/)

It uses a typical scientific python stack:
[Numpy](http://www.scipy.org)
[Scipy](http://www.scipy.org)
[matplotlib](http://matplotlib.org)

## Notes
This software was originally an attempt to help me get used to programming in Python environment.  
I have tested more since it was originally developed and have modified accordingly.

Please feel free to contact me with questions or suggestions.  
Do note this is a side project and it may take me some time to respond.

This is open-source software, with no warranties extended.

Nick Guy (nick.guy@noaa.gov)
