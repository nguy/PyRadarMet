PyRadarMet
=============== 

Python Fundamental Calculations in Radar Meteorology package notes

Created:   5 February 2014	Nick Guy (NRC; NOAA/NSSL)
Updated:  13 February 2014      NG - Added functions to system and variables
Updated:   8 April 2014         NG - Added extensive docstring formatting, correct Vmax
                                      equation.
          19 June 2014          NG - Added a script to list raw Sigmet characteristics


Further Details::
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
              
BeamBlock - A class that allows the calculation of the geometric beam blocking .

## Installation
The setup.py should now work.  The package has been put in a more standard format as of 21 Nov 2014.
Install should be able to be done by the following

e.g.
```python
python setup.py install
```

The install takes some time because of all of the DEM files in the data directory.

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
This is somewhat immature software, as it was originally an attempt to help me get  
 used to programming in Python environment.  
I have tested more since it was originally developed and have modified accordingly.

Please feel free to contact me with questions or suggestions.  
Do note this is a side project and it may take me some time to respond.

This is open-source software, with no warranties extended.

Nick Guy (nick.guy@noaa.gov)
