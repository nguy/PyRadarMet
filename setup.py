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

import os
import sys

from setuptools import setup, find_packages
from distutils.sysconfig import get_python_lib

#- Pull the header into a variable 
doclines = __doc__.split("\n")

VERSION = '0.1.1'

DEM_DATADIR = os.sep.join([os.path.dirname(__file__), 'data'])

#- Set variables for setup
PACKAGES = ['pyradarmet']

package_dirs={'pyradarmet'}
#datafiles = glob.glob(os.path.join(pathout,'*'))
#datafiles = [os.path.join('data',os.path.basename(f)) for f in datafiles]
#package_data = {'pyradarmet':datafiles}

#- Run setup
setup(
      name='pyradarmet',
      version=VERSION,
      url='https://github.com/nguy/PyRadarMet',
      author='Nick Guy',
      author_email='nick.guy@noaa.gov',
      description=doclines[0],
      license='LICENSE.txt',
      packages=PACKAGES,
      package_data={'pyradarmet': ['data/*']},
      include_package_data=True,
      classifiers=["""
        Development Status :: 3 - Alpha,
        Programming Language :: Python",
        Topic :: Scientific/Engineering
        Topic :: Scientific/Engineering :: Atmospheric Science
        Operating System :: Unix
        Operating System :: POSIX :: Linux
        Operating System :: MacOS
        """],
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
      install_requires = ['Numpy >=1.7.2','pyart >= 1.0.0'],
      )
