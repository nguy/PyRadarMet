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
import urllib2
import shutil
#from shutil import copyfile
import tarfile
import time
from setuptools import setup

#- Pull the header into a variable 
doclines = __doc__.split("\n")

VERSION = '0.1.2'

#- Set variables for setup
PACKNAME = 'pyradarmet'
PACKAGES = [PACKNAME]
URL = 'https://github.com/nguy/PyRadarMet'
AUTHOR = 'Nick Guy'
EMAIL = "nick.guy@uwyo.edu"
LICENSE = 'LICENSE.txt'

#DEMDIR = PACKNAME + os.sep.join([os.path.dirname(__file__), 'data'])

package_dirs={'':PACKNAME}
#datafiles = glob.glob(os.path.join(pathout,'*'))
#datafiles = [os.path.join('data',os.path.basename(f)) for f in datafiles]
#package_data = {'pyradarmet':datafiles}

# Check for Pyart package - Avoid conflict with PyPI

try:
    import pyart
except ImportError as e:
    print("Could not find PyArt Install\n\
           Please go to https://github.com/ARM-DOE/pyart for PyArt install")
    sys.exit()

# Read the text file with location of where to install
#demfiletxt = open('pyradarmet/demfile_locator.txt', 'r')
#DEMInstall = demfiletxt.readline()
#demfiletxt.close()

# Check for trailing directory "/"
#if DEMInstall[-1] is not '/':
#        DEMInstall = DEMInstall + '/'
########

gtopourl="ftp://edcftp.cr.usgs.gov/data/gtopo30/global/"

try:
  basefolder = os.environ["GTOPO_DATA"]
  print "Topo data will be installed to %s" % os.environ["GTOPO_DATA"]
except KeyError:
  basefolder = "./pyradarmet/"
  print "GTOPO_DATA environmental variable not set, "\
  "using default path in package directory"

archfolder=os.path.join(basefolder, 'gz')
datafolder=os.path.join(basefolder, 'data')

names =["antarcps", "e060n40", "e100n40", "e140n40", "w020n40", "w060n90",
        "w100n90", "w140n90", "w180s10", "e020n40", "e060n90", "e100n90",
        "e140n90", "w020n90", "w060s10", "w100s10", "w140s10", "w180s60",
        "e020n90", "e060s10", "e100s10", "e140s10", "w020s10", "w060s60",
        "w120s60", "w180n40", "e020s10", "e060s60", "e120s60", "w000s60",
        "w060n40", "w100n40", "w140n40", "w180n90"]

for folder in [archfolder, datafolder]:
  print folder
  if not os.path.exists(folder):
    os.makedirs(folder)

start_time = time.time()

for filename in names:
  archfilebase = filename+".tar.gz"
  archfile = os.path.join(archfolder, archfilebase)
  hdrfilebase = filename+".hdr"
  demfilebase = filename+".dem"

  if not os.path.exists(archfile):
    print("downloading: {0}".format(archfilebase))
    ftpfile = urllib2.urlopen(gtopourl+archfilebase)
    localfile = open(archfile, "wb")
    shutil.copyfileobj(ftpfile, localfile)
    localfile.close()

  for exfilebase in [hdrfilebase, demfilebase]:
    exfile = os.path.join(datafolder,exfilebase)
    if not os.path.exists(exfile):
      print("extracting: {0} from {1}".format(exfilebase, archfilebase))
      tar = tarfile.open(archfile)  
      tarobj = tar.extractfile(exfilebase.upper())
      localfile = open(exfile, "wb")
      shutil.copyfileobj(tarobj, localfile)
      tarobj.close()
      localfile.close()
      tar.close()

download_time = time.time() - start_time
print "Total download time: %g seconds"%(download_time)
#print "Removing tar directory: %s" % archfolder
#shutil.rmtree(archfolder)
#########
#- Run setup
if datafolder == './pyradarmet/data/':
	setup(
		  name=PACKNAME,
		  version=VERSION,
		  url=URL,
		  author=AUTHOR,
		  author_email=EMAIL,
		  description=doclines[0],
		  license=LICENSE,
		  packages=PACKAGES,
		  package_data={PACKNAME: [datafolder+'*']},
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
		  install_requires = ['Numpy >=1.7.2'],
		  )
else:
    print "In EXTERNAL"

    setup(
		  name=PACKNAME,
		  version=VERSION,
		  url=URL,
		  author=AUTHOR,
		  author_email=EMAIL,
		  description=doclines[0],
		  license=LICENSE,
		  packages=PACKAGES,
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
		  install_requires = ['Numpy >=1.7.2'],
		  )
		  

install_time = time.time() - start_time
print "Total install time: %g seconds"%(install_time)
