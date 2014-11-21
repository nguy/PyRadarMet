# -*- coding: utf-8 -*-
"""
BeamBlock.py - Class for Beam Blockage calculations

Author::
Nick Guy - OU CIMMS/Univ of Miami

This program was ported from code written in IDL used at Colorado State University.
It is believed that the original code at CSU was written by Steven Nesbitt.  
Timothy Lang also contributed to the development of this program.

This program is not particularly fast as it reads in large DEM files.
"""
# Import required libraries
import struct

import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from glob import glob
import os

import pyradarmet as rmet
from .geometry import r_effective, ray_height
from .geometry import half_power_radius, range_correct, beam_block_frac
from scipy.io.idl import readsav
import scipy.ndimage as scim
from mpl_toolkits.basemap import Basemap

###################################

VERSION = '0.0.1'

##################################
# CONSTANTS #
#############
# Radius of earth (m)
RE = 6371000.

# Plotting setup
TITLEDICT = {'fontsize': 10}
PAN3_AX1 = [.08, .5, .43, .43]
PAN3_AX2 = [.58, .5, .43, .43]
PAN3_AX3 = [.10, .08, .80, .35]

# Set up the data directory for DEM files
DEM_DATADIR = os.sep.join([os.path.dirname(__file__), 'data'])

class BeamBlock(object):
    """
    A class instance to calculalate beam blockage of weather radar
    """

###################################################

    def __init__(self, dem_filename=None, verbose=False,
                 dem_hdr_file=None, upper_left_x=None, upper_left_y=None,
                 radar_lon=None, radar_lat=None, radar_alt=None,
                 radar_antenna_height=None,
                 elev_angle=None, beamwidth=None, refract_grad=None,
                 range_dist=None, azimuthal_angle=None,
                 num_rng_pts=None, num_az_pts=None,
                 lon_box_size=6., lat_box_size=6.):
        """
        If initialized with a filename (incl. path), will call
        read_bin_netcdf() to populate the class instance.
        
        If not, it simply instances the class but does not populate
        its attributes.
        
        Parameters::
        ----------
        dem_filename : str
            Full path and filename of Digital Elevation Model file to use.
            If not supplied, a file will be selected, along with the 
            appropriate dem_hdr_file, based on radar_lon and radar_lat
        dem_hdr_file : str
            Full path and filename of Digital Elevation Model header file to use.
            This file will find upper_left_x and upper_left_y. 
            If not set, upper_left_x and upper_left_y must be set by user.
        upper_left_x : float
            Coordinate of upper left X position of file [decimal degrees]
            Only required if dem_hdr_file not supplied
        upper_left_y : float
            Coordinate of upper left Y position of file [decimal degrees]
        verbose: boolean
            Set to True for text output. Useful for debugging.
        
        radar_lon : float
            Radar longitude [decimal degrees]
        radar_lat : float
            Radar latitude [decimal degrees]
        radar_alt : float
            Radar altitude [meters]
        radar_antenna_height : float
            Radar antenna height [meters]
        elev_angle : float
            Elevation angle for analysis [degrees]
        beamwidth : float
            Radar beam width [degrees]
        refract_grad : float
            Vertical gradient of refraction [1/km]
        range_dist : float
            Radial distance of radar beam [km]
        azimuthal angle : float
            Azimuthal angle [decimal degrees] to display beam propagation along
        num_rng_pts : int
            Number of points along the range direction (i.e. # gates)
        num_az_pts : int
            Number of points along the azimuthal direction (i.e. # angles)
        lon_box_size : float
            The size of the longitudinal box to subset, radar in middle
        lat_box_size : float
            The size of the latitudinal box to subset, radar in middle
        """
        self.dem_file = dem_filename
        self.dem_hdr_file = dem_hdr_file
        self.upper_left_x = upper_left_x
        self.upper_left_y = upper_left_y
        self.lon_box_size = lon_box_size
        self.lat_box_size = lat_box_size
        self.ralt = radar_alt
        
        # Set radar latitude and longitude
        if (radar_lon is None) or (radar_lat is None):
            print "Radar longitude, latitude, and altitude need to be defined!"
            exit()
        else:
            self.rlon, self.rlat= radar_lon, radar_lat
        
        # Set the height of the radar antenna above ground
        if radar_antenna_height is None:
            self.rheight = 0.
        else:
            self.rheight = radar_antenna_height
            
        # Set elevation angle (degrees) [e.g. 0.8, 1.3, 1.8, 2.5]
        if elev_angle is None:
            self.E = 0.8
        else:
            self.E = elev_angle

        # Set radar beam width (degrees)
        if beamwidth is None:
            self.BW = 0.91
        else:
            self.BW = beamwidth

        # Set the vertical refractivity gradient (1/km) [Standard Atm = -40]
        if refract_grad is None:
            self.dNdH = -40.
        else:
            self.dNdH = refract_grad

        # Set the radial distance of radar [range] in km
        if range_dist is None:
            self.range = 150.
        else:
            self.range = range_dist

        # Azimuthal direction to plot beam propagation
        if azimuthal_angle is None:
            self.Az_plot = 0
        else:
            self.Az_plot = azimuthal_angle

        # Set the number of points along the radial
        if num_rng_pts is None:
            self.nrng = 601
        else:
            self.nrng = num_rng_pts
        
        # Set the number of azimuthal points
        if num_az_pts is None:
            self.naz = 361
        else:
            self.naz = num_az_pts
        
        # Check the DEM file
        self._check_dem_file()
        
        # Read the DEM file
        self.read_dem()
        
        # Process inputs to calculate beam properties
        self.calc_beam_blockage()
            
    ##############

    def help(self):

        _method_header_printout('help')
        print 'To define a new MosaicTile(), use instance = MosaicTile().'
        print 'then read in a file to populate attributes.'
        print 'Available read methods:'
        print '    read_mosaic_netcdf(<FILE>):'
        print '    read_mosaic_binary(<FILE>):'
        print '    Read v1 MRMS mosaics (7/30/2013 & earlier)'
        print '    Also v2 MRMS mosaics (7/30/2013 & after)'
        print 'Can also use instance = MosaicTile(filepath+name).'
        print 'Set binary=True as keyword above if reading an MRMS binary file.'
        print 'Other available methods:'
        print 'diag(), get_comp(), plot_vert(), plot_horiz(), three_panel_plot()'
        print 'subsection(), write_mosaic_binary(), output_composite()'
        _method_footer_printout()
        
##########################
# Read data file methods #
##########################

    def read_dem(self):
        '''Read a DEM file if passed or search for correct file
        '''
        if self.dem_hdr_file is None:
            if (upper_left_x is None) or (upper_left_y is None):
                print "Must prived upper_left_x and upper_left_y if \
                       no DEM header file is supplied"
        self.read_gt30_dem()

    def read_gt30_dem(self):
        '''
        Read in a binary GTOPO30 DEM file
        The following code was adapted from an example at:
        http://stevendkay.wordpress.com/tag/dem/
        '''
        if isinstance(self.dem_hdr_file, str) != False:
            hdr = np.loadtxt(self.dem_hdr_file, dtype='S20')
            nrows = int(hdr[2][1])
            ncols = int(hdr[3][1])
            missing = hdr[9][1]
            upper_left_x = float(hdr[10][1])
            upper_left_y = float(hdr[11][1])
            xdim = float(hdr[12][1])
            ydim = float(hdr[13][1])
        else:
            # From the output_parameters.hdr file or product documentation
            upper_left_x = self.upper_left_x
            upper_left_y = self.upper_left_y
            nrows = 6000
            ncols = 4800
            missing = -9999.
            xdim = 0.00833333333333
            ydim = 0.00833333333333
            
        # Extract GTOPO30 BIL file
        fi = open(self.dem_file, "rb")
        databin = fi.read()
        fi.close()
            
        # Unpack binary data into a flat tuple z
        s = ">%dH" % (nrows * ncols,)
        z = struct.unpack(s, databin)
        
        topo = np.zeros((nrows, ncols))
        for r in range(0, nrows):
            for c in range(0, ncols):
                elevation = z[((ncols) * r) + c]
                if (elevation == 65535 or elevation < 0 or elevation > 20000):
                    # may not be needed depending on format, and the "magic number"
                    # value used for 'void' or missing data
                    elevation = 0.#np.nan
                topo[r][c] = float(elevation)

        lat = upper_left_y - ydim * np.arange(nrows)
        lon = upper_left_x + xdim * np.arange(ncols)
        
        # Find indices for a subset of DEM data
        londel = self.lon_box_size / 2.
        latdel = self.lat_box_size / 2.
        lonInd = np.where((lon >= (self.rlon - londel)) & \
                          (lon <= (self.rlon + londel)))
        latInd = np.where((lat >= (self.rlat - latdel)) & \
                          (lat <= (self.rlat + latdel)))
        
        lonSub = lon[np.min(lonInd):np.max(lonInd)]
        latSub = lat[np.min(latInd): np.max(latInd)]
        self.lon, self.lat = np.meshgrid(lonSub, latSub)
        self.topo = topo[np.min(latInd):np.max(latInd), np.min(lonInd): np.max(lonInd)]
                            
        self.topo = np.ma.masked_outside(self.topo, 0, 20000)
        self.topo = np.ma.masked_invalid(self.topo)
        
        # Find the elevation of the radar
        if self.ralt is None:
            radar_lon_index = self._get_array_index(self.rlon, self.lon[0, :])
            radar_lat_index = self._get_array_index(self.rlat, self.lat[:, 0])
        
            radar_alt_map = scim.map_coordinates(self.topo,
                                np.vstack((radar_lat_index, radar_lon_index)),
                                prefilter=False)
            self.ralt = self.rheight + radar_alt_map

###################################################
# Calculation methods #
###################################################

    def calc_beam_blockage(self):
        '''Calculate the properties of the beam for calculation'''
        # Initialize arrays
        self.inX = np.empty(self.nrng)
        self.inY = np.empty(self.nrng)

        self.rng_lon = np.empty([self.naz, self.nrng])
        self.rng_lat = np.empty([self.naz, self.nrng])
        self.terr = np.empty([self.naz, self.nrng])
        
        # Points along each ray to calculate beam propagation
        # e.g. From 0 to 150 km, with 0.25 km bin spacing
        self.rng = np.linspace(0, self.range*1000., self.nrng)
        
        # Calculate the effective radius
        self.Reff = r_effective()

        # Calculate the height of center of beam
        self.h = ray_height(self.rng, self.E, 
                                          self.ralt, R1=self.Reff)

        # Calculate the beam radius
        self.a = half_power_radius(self.rng, self.BW)

        # Calculate rng_gnd (actual distance along ground)
        self.rng_gnd = range_correct(self.rng, self.h, self.E)

        # Set up arrays for calculations
        self.phis = np.linspace(0, 360, self.naz)
        
        # Initialize the beam blockage arrays with 0
        self.PBB = np.zeros((self.naz, self.nrng))
        self.CBB = np.zeros((self.naz, self.nrng))

        for jj, az in enumerate(self.phis):
            # Calculate the Lons/Lats along the rng_gnd
            self.rng_lat[jj, :] = np.degrees(np.arcsin(np.sin(np.radians(self.rlat)) * \
                                            np.cos(self.rng_gnd / RE) +\
                                            np.cos(np.radians(self.rlat)) * \
                                            np.sin(self.rng_gnd / RE) *\
                                            np.cos(np.radians(az))\
                                            )\
                                            )
            self.rng_lon[jj, :] = self.rlon + \
                     np.degrees(np.arctan2(np.sin(np.radians(az)) * \
                                np.sin(self.rng_gnd / RE) *\
                                np.cos(np.radians(self.rlat)), \
                                np.cos(self.rng_gnd / RE) - \
                                np.sin(np.radians(self.rlat)) * \
                                np.sin(np.radians(self.rng_lat[jj, :]))\
                                )\
                                )
                                
            # Find the indices for interpolation
            for ii, junk in enumerate(self.inX):
                self.inX[ii] = self._get_array_index(self.rng_lon[jj, ii], self.lon[0, :])
                self.inY[ii] = self._get_array_index(self.rng_lat[jj, ii], self.lat[:, 0])

            
            # Interpolate terrain heights to the radar grid
            self.terr[jj, :] = scim.map_coordinates(self.topo, 
                                                    np.vstack((self.inY, self.inX)),
                                                    prefilter=False)

            # Calculate PBB along range
            self.PBB[jj, :] = beam_block_frac(self.terr[jj, :], 
                                                            self.h, self.a)

        self.PBB = np.ma.masked_invalid(self.PBB)

        for ii in range(self.nrng):
            if ii == 0:
                self.CBB[:, ii] = self.PBB[:, ii]
            else:
                self.CBB[:,ii] = np.fmax([self.PBB[:, ii]], [self.CBB[:, ii-1]])
                
        self.terr = np.ma.masked_less(self.terr, -1000.)

###############
# Get methods #
###############

    def _get_array_index(self, value, array):
        '''Calculate the exact index position within latitude array'''
        # Find the spacing
        dp = np.absolute(array[1] - array[0])
        # Calculate the relative position
        pos = np.absolute(value - array[0]) / dp
        return pos

#################
# Check methods #
#################
    
    def _check_dem_file(self):
        '''Check if dem file is passed and set if not'''
        if self.dem_file is None:
            if (self.rlon >= -180.) and (self.rlon < 180.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30GT30antarcps'
            elif (self.rlon >= 20.) and (self.rlon < 60.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30e020n40'
            elif (self.rlon >= 20.) and (self.rlon < 60.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30e020n90'
            elif (self.rlon >= 20.) and (self.rlon < 60.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30e020s10'
                
            elif (self.rlon >= 60.) and (self.rlon < 100.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30e060n40'
            elif (self.rlon >= 60.) and (self.rlon < 100.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30e060n90'
            elif (self.rlon >= 60.) and (self.rlon < 100.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30e060s10'
            elif (self.rlon >= 60.) and (self.rlon < 120.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30e060s60'
                
            elif (self.rlon >= 100.) and (self.rlon < 140.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30e100n40'
            elif (self.rlon >= 100.) and (self.rlon < 140.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30e100n90'
            elif (self.rlon >= 100.) and (self.rlon < 140.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30e100s10'
            elif (self.rlon >= 120.) and (self.rlon < 180.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30e120s60'
                
            elif (self.rlon >= 140.) and (self.rlon < 180.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30e140n40'
            elif (self.rlon >= 140.) and (self.rlon < 180.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30e140n90'
            elif (self.rlon >= 140.) and (self.rlon < 180.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30e140s10'
                
            elif (self.rlon >= 0.) and (self.rlon < 60.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30w000s60'
                
            elif (self.rlon >= -20.) and (self.rlon < 20.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30w020n40'
            elif (self.rlon >= -20.) and (self.rlon < 20.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30w020n90'
            elif (self.rlon >= -20.) and (self.rlon < 20.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30w020s10'
                
            elif (self.rlon >= -60.) and (self.rlon < -20.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30w060n40'
            elif (self.rlon >= -60.) and (self.rlon < -20.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30w060n90'
            elif (self.rlon >= -60.) and (self.rlon < -20.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30w060s10'
            elif (self.rlon >= -60.) and (self.rlon < 0.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30w060s60'
                
            elif (self.rlon >= -100.) and (self.rlon < -60.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30w100n40'
            elif (self.rlon >= -100.) and (self.rlon < -60.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30w100n90'
            elif (self.rlon >= -100.) and (self.rlon < -60.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30w100s10'
            elif (self.rlon >= -120.) and (self.rlon < -60.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30w120s60'
                
            elif (self.rlon >= -140.) and (self.rlon < -100.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30w140n40'
            elif (self.rlon >= -140.) and (self.rlon < -100.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30w140n90'
            elif (self.rlon >= -140.) and (self.rlon < -100.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30w140s10'
                
            elif (self.rlon >= -180.) and (self.rlon < -140.) and \
               (self.rlat >= -10.) and (self.rlat < 40.):
                dem_filename = 'gt30w180n40'
            elif (self.rlon >= -180.) and (self.rlon < -140.) and \
               (self.rlat >= 40.) and (self.rlat < 90.):
                dem_filename = 'gt30w180n90'
            elif (self.rlon >= -180.) and (self.rlon < -140.) and \
               (self.rlat >= -60.) and (self.rlat < -10.):
                dem_filename = 'gt30w180s10'
            elif (self.rlon >= -180.) and (self.rlon < -120.) and \
               (self.rlat >= -90.) and (self.rlat < -60.):
                dem_filename = 'gt30w180s60'
                
            self.dem_file = DEM_DATADIR + "/" + dem_filename + ".dem"
            self.dem_hdr_file = DEM_DATADIR + "/" + dem_filename + ".hdr"
    
################
# Plot methods #
################

    def plot_3panel(self, beam_blockage='complete',
                    lon_plot_range=2., lat_plot_range=2.,
                    terrain_cmap=None, beam_blockage_cmap=None,
                    elev_min=None, elev_max=None,
                    range_rings=None,
                    saveFig=False):
        '''
        Create a 3-panel plot characterization of beam blockage
        
        Parameters::
        ----------
        beam_blockage : str
            'complete' or 'partial' chooses which field to plot
            Defaults to complete.
        lon_plot_range : float
            Horizontal coverage plots will show +/- lon_plot_range
        lat_plot_range : float
            Horizontal coverage plots will show +/- lat_plot_range
        terrain_cmap : str
            Matplotlib colormap to use for elevation map
        beam_blockage_cmap : str
            Matplotlib colormap to use for beam blockage map
        elev_min : float
            Minumum elevation to display on beam propagation and terrain height [meters]
        elev_max : float
            Maximum elevation to display on beam propagation and terrain height [meters]
        range_rings : float
            A list with location of range rings - e.g. [50., 100.]
        saveFig : boolean
            True to save figure as output file, False to show
        '''
        
        # Set the radar min and max longitude to plot
        self.minlat, self.maxlat = self.rlat - lat_plot_range, self.rlat + lat_plot_range
        self.minlon, self.maxlon = self.rlon - lon_plot_range, self.rlon + lon_plot_range
        if terrain_cmap is None:
            terrain_cmap = 'BrBG_r'
        if beam_blockage_cmap is None:
            beam_blockage_cmap = 'PuRd'
        self.terr_cmap, self.bb_cmap = terrain_cmap, beam_blockage_cmap
        
        fig = plt.figure(figsize=(10,8))
        ax1 = fig.add_axes(PAN3_AX1)
        ax2 = fig.add_axes(PAN3_AX2)
        ax3 = fig.add_axes(PAN3_AX3)
        
        self.draw_terrain_height_map(fig, ax1)
        if beam_blockage == 'partial':
            self.draw_bb_map(fig, ax2, BB=self.PBB, range_rings=range_rings)
        elif beam_blockage == 'complete':
            self.draw_bb_map(fig, ax2, range_rings=range_rings)
        self.draw_beam_terrain_profile(fig, ax3, ymin=elev_min, ymax=elev_max)
        
        if saveFig:
            plt.savefig('BBmap.png', format='png')
        else:
            plt.show()
    
    def draw_terrain_height_map(self, fig, ax):
        '''Draw the terrain heights'''
        bm1 = Basemap(projection='cea', resolution='l', area_thresh = 10000.,
                    llcrnrlon=self.minlon, urcrnrlon=self.maxlon,
                    llcrnrlat=self.minlat, urcrnrlat=self.maxlat,
                    ax=ax)
        ax.set_title('Terrain within 150 km of Radar (km)', fontdict=TITLEDICT)
        bm1.drawmeridians(np.arange(self.minlon, self.maxlon, 1.), labels=[1,0,0,1])           
        bm1.drawparallels(np.arange(self.minlat, self.maxlat, 1.), labels=[1,0,0,1])             
        bm1.drawcountries()
        bm1.drawcoastlines()
        bm1.drawrivers()
        
        xbm1, ybm1 = bm1(self.lon, self.lat)
        Htmap = bm1.pcolormesh(xbm1, ybm1, self.topo/1000., vmin=.05, vmax=3.,
                               cmap=self.terr_cmap)
        bm1.plot(self.rlon, self.rlat, 'rD', latlon=True)
        #plot_range_ring(50., bm=bm1, color='w')
        fig.colorbar(Htmap, ax=ax)
        
    def draw_bb_map(self, fig, ax, BB=None, range_rings=None):
        '''Draw the Beam Blockage'''
        if BB is None:
            BB = self.CBB
            
        bm2 = Basemap(projection='cea', resolution='l', area_thresh = 10000.,
                    llcrnrlon=self.minlon, urcrnrlon=self.maxlon,
                    llcrnrlat=self.minlat, urcrnrlat=self.maxlat,
                    ax=ax)
        ax.set_title('Beam-blockage fraction', fontdict=TITLEDICT)
        bm2.drawmeridians(np.arange(self.minlon, self.maxlon, 1.), labels=[1,0,0,1])           
        bm2.drawparallels(np.arange(self.minlat, self.maxlat, 1.), labels=[1,0,0,1])             
        bm2.drawcountries()
        bm2.drawcoastlines()
        bm2.drawrivers()

        xbm2, ybm2 = bm2(self.rng_lon, self.rng_lat)
        BBmap = bm2.pcolormesh(xbm2, ybm2, BB, vmin=0., vmax=1., cmap=self.bb_cmap)
        if range_rings is not None:
            for nn in range(len(range_rings)):
                self.plot_range_ring(range_rings[nn], bm=bm2)
        fig.colorbar(BBmap, ax=ax)
         
    def draw_beam_terrain_profile(self, fig, ax, ymin=None, ymax=None):
        '''Draw the Beam height along with terrain and PBB'''
        if ymin is None:
            ymin = 0.
        if ymax is None:
            ymax = 5000.
        
        bc, = ax.plot(self.rng_gnd / 1000., self.h / 1000., '-b', 
                     linewidth=3, label='Beam Center')
        b3db, = ax.plot(self.rng_gnd / 1000., (self.h + self.a) / 1000., ':b', 
                       linewidth=1.5, label='3 dB Beam width')
        ax.plot(self.rng_gnd / 1000., (self.h - self.a) / 1000., ':b')
        tf = ax.fill_between(self.rng_gnd / 1000., 0., 
                        self.terr[self.Az_plot, :] / 1000., 
                        color='0.75')
        ax.set_xlim(0., self.range)
        ax.set_ylim(ymin / 1000., ymax / 1000.)
        ax.set_title(r'dN/dh = %d km$^{-1}$; Elevation angle = %g degrees at %d azimuth'% \
                      (self.dNdH, self.E, self.Az_plot), fontdict=TITLEDICT)
        ax.set_xlabel('Range (km)')
        ax.set_ylabel('Height (km)')

        axb = ax.twinx()
        bbf, = axb.plot(self.rng_gnd / 1000., self.CBB[self.Az_plot, :], '-k', 
                       label='BBF')
        axb.set_ylabel('Beam-blockage fraction')
        axb.set_ylim(0., 1.)
        axb.set_xlim(0., self.range)
        
        ax.legend((bc, b3db, bbf), ('Beam Center', '3 dB Beam width', 'BBF'),
                                   loc='upper left', fontsize=10)

    def plot_range_ring(self, range_ring_location_km, bm=None,
                        color='k', ls='-'):
        """
        Plot a single range ring.
        Parameters::
        ----------
        range_ring_location_km : float
            Location of range ring in km.
        npts: int
            Number of points in the ring, higher for better resolution.
        ax : Axis
            Axis to plot on. None will use the current axis.
        """
        npts = 100
        bm.tissot(self.rlon, self.rlat, 
                  np.degrees(range_ring_location_km * 1000. / RE), npts,
                  fill=False, color='black', linestyle='dashed')
        
        