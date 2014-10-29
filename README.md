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

## Installation
As of now, I have used these as a standalone package that was not installed.
To access, I just added the path where the folder was unpacked to my 
PYTHONPATH environmental variable path (in my .bashrc file)

e.g.
```python
export PYTHONPATH=/Users/nickguy/programs/python/pythonlib
```

There is a setup.py file provided, this should work to install to a common python installation directory.

## Notes
This is somewhat immature software, as it was originally an attempt to help me get  
 used to programming in Python environment.  
I have tested more since it was originally developed and have modified accordingly.

Please feel free to contact me with questions or suggestions.  
Do note this is a side project and it may take me some time to respond.

This is open-source software, with no warranties extended.

Nick Guy (nick.guy@noaa.gov)
