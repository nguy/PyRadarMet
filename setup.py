"""PyRadarMet: Python Fundamental Calculations in Radar Meteorology 

PyRadarMet is a toolkit that contains a variety of utilities that can be used
 to calculate fundamental radar meteorology values:

* Convert between linear reflectivity(Z) and log reflectivity (dBZ)
* Calculate coefficients for attenuation calculations
* Calculate fundamental radar system characteristics
* Calculate Doppler radar characteristics
* Calculate geometrical characterisistics of radar
* Calculate variables from radar output
* Calculate ZDR bias of radar system

"""

from numpy.distutils.core import setup, Extension
import os

#- Pull the header into a variable 
doclines = __doc__.split("\n")

#- Set variables for setup
packages = ['pyradarmet']
package_dirs={'pyradarmet'}
datafiles = glob.glob(os.path.join(pathout,'*'))
datafiles = [os.path.join('data',os.path.basename(f)) for f in datafiles]
package_data = {'pyradarmet':datafiles}

#- Run setup
setup (name = 'pyradarmet',
       version = '0.1.1',
       author = 'Nick Guy',
       author_email = 'nick.guy@noaa.gov'
       packages = packages,
       package_dir = package_dirs,
       package_data = package_data,
       url = 'https://github.com/nguy/PyRadarMet',
       license='LICENSE.txt',
       description = doclines[0],
       long_description = """
A toolkit that contains a variety of utilities that can be used
 to calculate fundamental radar meteorology values:

* Convert between linear reflectivity(Z) and log reflectivity (dBZ)
* Calculate coefficients for attenuation calculations
* Calculate fundamental radar system characteristics
* Calculate Doppler radar characteristics
* Calculate geometrical characterisistics of radar
* Calculate variables from radar output
* Calculate ZDR bias of radar system""",
       install_requires = ['Numpy >=1.7.2','pyart'],
       )
