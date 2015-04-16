"""
pyradarmet.zdrcal
=========================

Calculate the differential reflectivity (Zdr) offset from polarimetric radars.

Adapted by Nick Guy from original Fortran code written by David Jorgensen, 
  implemented more robust processing to reduce noise in data.  Thanks to Joseph
  Hardin regarding processing and posting a notebook online that got this started.
  http://www.josephhardinee.com/blog/?p=35


.. autosummary::
    :toctree: generated/

    calculate_zdr_offset


"""
import numpy as np
import matplotlib.pylab as plt

import pyart
#from .conversion import dBZ2Z, Z2dBZ
#======================================================================
def calculate_zdr_offset(radar, debug=False, 
                          remove_first_n_gates=None, sample_min=100,
                          rhv_min=None, sig_min=0.5,Htmax=None, dbz_min=None,
                          zdr_field=None, refl_field=None, phidp_field=None,
                          rhv_field=None, kdp_field=None):
    """
    Differential reflectivity calibration from 'birdbath' scans.

    Parameters
    ----------
    radar : Radar
        Radar object (from PyArt) to use for attenuation calculations.  Must have
        copol_coeff, norm_coherent_power, proc_dp_phase_shift,
        reflectivity_horizontal fields.
    debug : bool
        True to print debugging information, False supressed this printing.

    Returns
    -------
    zdr_bias : dict
        Field dictionary containing the differential reflectivity statistics.

    Other Parameters
    ----------------
    remove_first_n_gates : float
        Number of gates at the beginning of each ray to to remove from the
        calculation.
    sample_min : int
        The minimum number of samples required to perform calculation.
    rhv_min : float
        Minimum copol_coeff value to consider valid.
    sig_min : float
        Minimum standard deviation to consider valid and include in calculations.
    dbz_min : float
        Minimum reflectivity to consider valid and include in calculations.
    Htmax : float
        Maximum height [km] to use in bias calculations.
    zdr_field, refl_field, rhv_field, phidp_field, kdp_field : str
        Field names within the radar object which represent the horizonal
        reflectivity, normal coherent power, the copolar coefficient, and the
        differential phase shift. A value of None for any of these parameters
        will use the default field name as defined in the Py-ART
        configuration file.

    Usage
    ----------
    
    """
    # Pull out some information to send back in zdr_bias dictionary   
    Rlat = radar.latitude['data'] # Radar latitude
    Rlon = radar.longitude['data'] # Radar longitude
    RAlt = radar.altitude['data'] # Radar altitude
    RadName = radar.metadata['instrument_name'] # Radar name
    SitName = radar.metadata['source'] # Site/project name
    GenName = radar.metadata['institution'] # Facilitiy name
    time = radar.time # Time at middle of each ray
    range = radar.range # Range along ray
    azimuth = radar.azimuth # Azimuth of scan
    
    # Create 2D arrays for masking data later
#    Ht2D, az2D = np.meshgrid(range['data']/1000.,azimuth['data'])
#    print radar.fields.keys()
    # parse the field parameters
    if zdr_field is None:
        zdr_field = pyart.config.get_field_name('differential_reflectivity')
        if zdr_field not in radar.fields:
           zdr_field = pyart.config.get_field_name('ZDR')
    if refl_field is None:
        refl_field = pyart.config.get_field_name('reflectivity')
        if refl_field not in radar.fields:
           refl_field = pyart.config.get_field_name('DBZ')
    if rhv_field is None:
        rhv_field = pyart.config.get_field_name('cross_correlation_ratio')
        if rhv_field not in radar.fields:
           rhv_field = pyart.config.get_field_name('RHOHV')
    if phidp_field is None:
        # use corrrected_differential_phae or unfolded_differential_phase
        # fields if they are available, if not use differential_phase field
        phidp_field = pyart.config.get_field_name('corrected_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = pyart.config.get_field_name('unfolded_differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = pyart.config.get_field_name('differential_phase')
        if phidp_field not in radar.fields:
            phidp_field = pyart.config.get_field_name('PHIDP')
    if kdp_field is None:
        kdp_field = pyart.config.get_field_name('specific_differential_phase')

    # Extract data fields
    ZDR = radar.fields[zdr_field]['data']
    dBZ = radar.fields[refl_field]['data']
    Rhohv = radar.fields[rhv_field]['data']
    PhiDP = radar.fields[phidp_field]['data']
    KDP = radar.fields[kdp_field]['data']
        
    # Extract parameters from radar
#    nsweeps = int(radar.nsweeps)
    nGates = int(radar.ngates) # Number of gates
    nrays = int(radar.nrays) # Number of rays

    # Apply mask across variables for missing data
    ZDR, dBZ, Rhv, PhiDP, KDP = mask_variables(ZDR, dBZ, Rhohv, PhiDP, KDP)

    # Find corresponding axes for range (height) and azimuth        
    axAz, axHt = get_ray_range_dims(ZDR, nrays)
    
    # Apply mask for chosen gates nearest to radar and above height threshold
    ZDR, dBZ, Rhv, PhiDP, KDP = mask_gates(ZDR, dBZ, Rhv, PhiDP,KDP,
                                    range['data']/1000., azimuth['data'], axAz,
                                    Htmax=Htmax, remove_first_n_gates=remove_first_n_gates)
                                    
    # Apply mask for reflectivity below a threshold (proxy for low SNR)
    ZDR, dBZ, Rhv, PhiDP, KDP = mask_refl(ZDR, dBZ, Rhv, PhiDP, KDP, mindBZ=dbz_min)

    # Find individual ray properties (along each ray at azimuthal step)
    DR_RaySum, DR_RayAvg, DR_RayStd, NGood_Ray = calc_stats(ZDR, axAz)
    
    # Find individual range properties (along a constant distance [height] from radar)
    DR_HtSum, DR_HtAvg, DR_HtStd, NGood_Ht = calc_stats(ZDR, axHt)
    RH_HtSum, RH_HtAvg, RH_HtStd, NGood_RH_Ht = calc_stats(Rhv, axHt)
    dBZ_HtSum, dBZ_HtAvg, dBZ_HtStd, NGood_dBZ_Ht = calc_stats(dBZ, axHt)
    KDP_HtSum, KDP_HtAvg, KDP_HtStd, NGood_KDP_Ht = calc_stats(KDP, axHt)

    # Apply mask for data with range averaged Std Dev and azimuthal averaged Std Dev
    # greater than threshold
    ZDR, dBZ, Rhv, PhiDP, KDP = mask_sigma(ZDR, dBZ, Rhv, PhiDP, KDP, DR_HtStd, DR_RayStd,
                                    axAz, minSig=sig_min)
    
    # Find volume properites    
    DR_Avg = np.ma.mean(ZDR)
    DR_Std = np.ma.std(ZDR, ddof=1)
    NGood = np.ma.count(ZDR)
    
    if (NGood > sample_min):
        # Find statistics for filtered (masked) variables
        DR_HtSum_filt, DR_HtAvg_filt, DR_HtStd_filt, NGood_Ht_filt = calc_stats(ZDR, axHt)
        RH_HtSum_filt, RH_HtAvg_filt, RH_HtStd_filt, NGood_RH_Ht_filt = calc_stats(Rhv, axHt)
        dBZ_HtSum_filt, dBZ_HtAvg_filt, dBZ_HtStd_filt, NGood_dBZ_Ht_filt = calc_stats(dBZ, axHt)
        KDP_HtSum_filt, KDP_HtAvg_filt, KDP_HtStd_filt, NGood_KDP_Ht_filt = calc_stats(KDP, axHt)
        PhiDP_HtSum_filt, PhiDP_HtAvg_filt, PhiDP_HtStd_filt, NGood_PhiDP_Ht_filt = calc_stats(PhiDP, axHt)
    
        # Create an output dictionary
        zdr_bias = create_dict(DR_Avg, DR_Std, NGood, 
                    DR_RaySum, DR_RayAvg, DR_RayStd, NGood_Ray,
                    DR_HtAvg, DR_HtStd, NGood_Ht,
                    RH_HtAvg, dBZ_HtAvg, KDP_HtAvg,
                    ZDR, DR_HtAvg_filt, DR_HtStd_filt,
                    RH_HtAvg_filt, dBZ_HtAvg_filt, KDP_HtAvg_filt, PhiDP_HtAvg_filt,
                    RadName, GenName, Rlat, Rlon, RAlt, time, range, azimuth)
    
    else:
        zdr_bias = None
    # Send back the dictionary containing data.
    return zdr_bias
#======================================================================
def mask_variables(ZDR, dBZ, RhoHV, PhiDP, KDP, rhv_min=0.8):
    # Combine the masks for all variables to ensure no missing data points
    maskZDR = np.ma.getmask(ZDR)
    maskZ = np.ma.getmask(dBZ)
    maskRhv = np.ma.getmask(RhoHV)
    maskPdp = np.ma.getmask(PhiDP)
    is_cor = RhoHV < rhv_min # Mask values where Correl Coeff < threshold 
    
    mskt1 = np.logical_and(maskZ, maskZDR) # Combine refl and ZDR masks
    mskt2 = np.logical_and(mskt1, maskRhv) # Combine with missing Rho HV mask
    mskt3 = np.logical_and(mskt2, is_cor) # Combine with min threshold Rho HV
    mask = np.logical_and(mskt3, maskPdp) # Combine with missing Phidp mask
    
    ZDR = np.ma.masked_where(mask, ZDR)
    dBZ = np.ma.masked_where(mask, dBZ)
    Rhv = np.ma.masked_where(mask, RhoHV)
    PhiDP = np.ma.masked_where(mask, PhiDP)
    KDP = np.ma.masked_where(mask, KDP)
    
    return ZDR, dBZ, Rhv, PhiDP, KDP
#======================================================================
def get_ray_range_dims(ZDR, nrays):
    if nrays == ZDR.shape[0]:
        axAz = 1
        axHt = 0
    elif nrays == ZDR.shape[1]:
        axAz = 0
        axHt = 1
    else:
        print "Makes sure you are using polar coordinate data!"
        return
    return axAz, axHt
#======================================================================
def mask_gates(ZDR, dBZ, Rhv, PhiDP, KDP, 
               Ht, Az, axAz, Htmax=10., remove_first_n_gates=None):
    # Create 2D arrays for masking data later
    Ht2D, az2D = np.meshgrid(Ht, Az)
    
    # Check to see what gate to start with user can set or default to 0
    if remove_first_n_gates is None:
        GateBeg=1
    elif remove_first_n_gates == 0:
        GateBeg=1
    else:
        GateBeg=remove_first_n_gates # because of 0 based indexing

    # Check which axis to perform calculations over
    # Mask out the lower gates based upon remove_first_n_gates keyword    
    if axAz == 1:
            # Mask lower gates
        ZDR[:, 0:GateBeg-1] = np.ma.masked
        dBZ[:, 0:GateBeg-1] = np.ma.masked
        Rhv[:, 0:GateBeg-1] = np.ma.masked
        PhiDP[:, 0:GateBeg-1] = np.ma.masked
        KDP[:, 0:GateBeg-1] = np.ma.masked
        
        # Mask upper gates
        ZDR = np.ma.masked_where(Ht2D > Htmax, ZDR)
        dBZ = np.ma.masked_where(Ht2D > Htmax, dBZ)
        Rhv = np.ma.masked_where(Ht2D > Htmax, Rhv)
        PhiDP = np.ma.masked_where(Ht2D > Htmax, PhiDP)
        KDP = np.ma.masked_where(Ht2D > Htmax, KDP)
    elif axAz == 0:
        ZDR[0:GateBeg-1, :] = np.ma.masked
        dBZ[0:GateBeg-1, :] = np.ma.masked
        Rhv[0:GateBeg-1, :] = np.ma.masked
        PhiDP[0:GateBeg-1, :] = np.ma.masked
        KDP[0:GateBeg-1, :] = np.ma.masked
        
        # Mask upper gates
        ZDR = np.ma.masked_where(Ht2D.transpose() > Htmax, ZDR)
        dBZ = np.ma.masked_where(Ht2D.transpose() > Htmax, dBZ)
        Rhv = np.ma.masked_where(Ht2D.transpose() > Htmax, Rhv)
        PhiDP = np.ma.masked_where(Ht2D.transpose() > Htmax, PhiDP)
        KDP = np.ma.masked_where(Ht2D.transpose() > Htmax, KDP)
        
    return ZDR, dBZ, Rhv, PhiDP, KDP
#======================================================================
def mask_sigma(ZDR, dBZ, Rhv, PhiDP, KDP, HtStd, AzStd, axAz, minSig=0.5):

    # Create 2D arrays for masking data later
    HtStd2D, RayStd2D = np.meshgrid(HtStd, AzStd)

    # Check which axis to perform calculations over
    # Mask out the lower gates based upon remove_first_n_gates keyword    
    if axAz == 1:
        # Mask along the range (height)
        ZDR = np.ma.masked_where(HtStd2D > minSig, ZDR)
        dBZ = np.ma.masked_where(HtStd2D > minSig, dBZ)
        Rhv = np.ma.masked_where(HtStd2D > minSig, Rhv)
        PhiDP = np.ma.masked_where(HtStd2D > minSig, PhiDP)
        KDP = np.ma.masked_where(HtStd2D > minSig, KDP)
        
        # Mask along the Azimuths
        ZDR = np.ma.masked_where(RayStd2D > minSig, ZDR)
        dBZ = np.ma.masked_where(RayStd2D > minSig, dBZ)
        Rhv = np.ma.masked_where(RayStd2D > minSig, Rhv)
        PhiDP = np.ma.masked_where(RayStd2D > minSig, PhiDP)
        KDP = np.ma.masked_where(RayStd2D > minSig, KDP)
    elif axAz == 0:
        # Mask along the range (height)
        ZDR = np.ma.masked_where(HtStd2D.transpose() > minSig, ZDR)
        dBZ = np.ma.masked_where(HtStd2D.transpose() > minSig, dBZ)
        Rhv = np.ma.masked_where(HtStd2D.transpose() > minSig, Rhv)
        PhiDP = np.ma.masked_where(HtStd2D.transpose() > minSig, PhiDP)
        KDP = np.ma.masked_where(HtStd2D.transpose() > minSig, KDP)
        
        # Mask along the Azimuths
        ZDR = np.ma.masked_where(RayStd2D.transpose() > minSig, ZDR)
        dBZ = np.ma.masked_where(RayStd2D.transpose() > minSig, dBZ)
        Rhv = np.ma.masked_where(RayStd2D.transpose() > minSig, Rhv)
        PhiDP = np.ma.masked_where(RayStd2D.transpose() > minSig, PhiDP)
        KDP = np.ma.masked_where(RayStd2D.transpose() > minSig, KDP)
        
    return ZDR, dBZ, Rhv, PhiDP, KDP
#======================================================================
def mask_refl(ZDR, dBZ, Rhv, PhiDP, KDP, mindBZ=-40.):
    # Mask out the points with reflectivity < mindBZ
    ZDR = np.ma.masked_where(dBZ < mindBZ, ZDR)
    dBZ = np.ma.masked_where(dBZ < mindBZ, dBZ)
    Rhv = np.ma.masked_where(dBZ < mindBZ, Rhv)
    PhiDP = np.ma.masked_where(dBZ < mindBZ, PhiDP)
    KDP = np.ma.masked_where(dBZ < mindBZ, KDP)
        
    return ZDR, dBZ, Rhv, PhiDP, KDP
#======================================================================
def calc_stats(Var, Ax):
    Sum = np.ma.cumsum(Var, axis=Ax)
    Avg = np.ma.mean(Var, axis=Ax)
    Std = np.ma.std(Var, axis=Ax)
    nValid = np.ma.count(Var, axis=Ax)
    
    return Sum, Avg, Std, nValid
#======================================================================
def create_dict(DR_Avg, DR_Std, NGood, 
                DR_RaySum, DR_RayAvg, DR_RayStd, NGood_Ray,
                DR_HtAvg, DR_HtStd, NGood_Ht,
                RH_HtAvg, dBZ_HtAvg, KDP_HtAvg,
                ZDR, DR_HtAvg_filt, DR_HtStd_filt,
                RH_HtAvg_filt, dBZ_HtAvg_filt, KDP_HtAvg_filt, PhiDP_HtAvg_filt,
                RadName, GenName, Rlat, Rlon, RAlt, time, range, azimuth):
                
    zdr_bias = {'volume_average': DR_Avg, 
                'volume_standard_deviation': DR_Std,
                'volume_number_good': NGood,
                'ray_cumulative_sum': DR_RaySum, 
                'ray_average': DR_RayAvg,
                'ray_standard_deviation': DR_RayStd, 
                'ray_number_good': NGood_Ray,
                'height_average': DR_HtAvg,
                'height_standard_deviation': DR_HtStd, 
                'height_number_good': NGood_Ht, 
                'rhv_height_average': RH_HtAvg, 
                'dbz_height_average': dBZ_HtAvg,
                'kdp_height_average': KDP_HtAvg,
                'filtered_zdr': ZDR, 
                'height_average_filt': DR_HtAvg_filt,
                'height_standard_deviation_filt': DR_HtStd_filt,
                'rhv_height_average_filt': RH_HtAvg_filt, 
                'dbz_height_average_filt': dBZ_HtAvg_filt,
                'kdp_height_average_filt': KDP_HtAvg_filt,
                'phidp_height_average_filt': PhiDP_HtAvg_filt,
                'radar_name': RadName, 
                'facility_name': GenName,
                'radar_latitude':Rlat, 
                'radar_longitude':Rlon, 
                'radar_altitude':RAlt,
                'time': time, 
                'range': range, 
                'azimuth':azimuth}
    
    return zdr_bias
#======================================================================
def plot_zdr_rhv(X1, X2, Y, xlims1=None, xlims2=None, ymax=None,
               label=True, labelTx=' ', labelpos=None):
    """Create a vertical plot of differential reflectivity and copolar correlation
    coefficient.
    """
    fig, ax = plt.subplots()
    axes = [ax, ax.twiny()]
    
    if xlims1 is None:
      xlims1 = (-1.5, 1.5)
    if xlims2 is None:
      xlims2 = (.8, 1.)
    if ymax is None:
      ymax = 10.
      
    axes[0].plot(X1, Y, label='Differential Reflectivity', color='b')
    axes[1].plot(X2, Y, label='CoPolar Correlation Coefficient', color='r')
    axes[0].set_xlim(xlims1)
    axes[1].set_xlim(xlims2)
    axes[0].set_xlabel('Differential Reflectivity (dB)')
    axes[1].set_xlabel('CoPolar Correlation Coefficient')
    axes[0].tick_params(axis='x', colors='b')
    axes[1].tick_params(axis='x', colors='r')
    axes[0].xaxis.label.set_color('b')
    axes[1].xaxis.label.set_color('r')
    axes[0].set_ylim(0., ymax)
    axes[0].set_ylabel('Altitude (km)')
    axes[0].vlines(0, 0, ymax)
    axes[0].xaxis.grid(color='b', ls=':', lw=1)
    axes[1].xaxis.grid(color='r', ls=':', lw=1)
    axes[0].yaxis.grid(color='k', ls=':', lw=.5)
    
    handles1, labels1 = axes[0].get_legend_handles_labels()
    handles2, labels2 = axes[1].get_legend_handles_labels()
    plt.legend([handles1[0], handles2[0]],
            ['Differential Reflectivity','CoPolar Correlation Coefficient'],
            loc=6, fontsize =11)
    # Add label to lower left in plot unless location specified
    if label:
        if labelpos is None:
            labelpos = (.15, .15)
        plt.figtext(labelpos[0], labelpos[1], labelTx)

    return fig, ax, axes
#======================================================================
def oplot_zdr_rhv(X1, X2, Y, ax=None, alf=1.):
    """Add additional differential reflectivity and copolar correlation
    coefficient lines to existing graphic."""
    ax[0].plot(X1, Y, color='b', alpha=alf)
    ax[1].plot(X2, Y, color='r', alpha=alf)
#======================================================================
def plot_add_stats(instance=None, avg=None, sd=None, points=None,
                   good_samples=None, tot_samples=None, labelpos=None):
    """Add text to the figure that contains statistics."""
    if labelpos is None:
        labelpos=(.15, .7)
    if instance == 'single':
        plt.figtext(labelpos[0], labelpos[1],
            "Avg Zdr = %g\nStd Dev = %g\nPoints = %g\n"%(avg, sd, points))
    elif instance == 'multi':
        plt.figtext(labelpos[0], labelpos[1],
            "Avg Zdr = %g\nStd Dev = %g\nSamples = %g used of %g total\nPoints = %g\n"%
            (avg, sd, good_samples, tot_samples, points))
    else:
        print "Need to supply either single or multiple run statistics!!"
           
#======================================================================
def plot_3panel(X1, X2, X3, Y, 
               xlims1=None, xlims2=None, xlims3=None, ymax=None,
               X1lab='dBZ', X2lab='ZDR', X3lab=r'Kdp (dB km$^{-1}$)', ylab=None,
               label=True, labelTx=' ', labelpos=None):
    """Create a 3-panel vertical plot of reflectivity, differential reflectivity,
      and specific differential phase.
    """
    fig, (ax1, ax2, vax3) = plt.subplots(1, 3, sharey=True)
    
    if xlims1 is None:
      xlims1 = (-20, 30.)    
    if xlims2 is None:
      xlims2 = (-1.5, 1.5)
    if xlims3 is None:
      xlims3 = (.2, .4)
    if ymax is None:
      ymax = 10.
    ax1.plot(X1, Y, color='k')
    ax2.plot(X2, Y, color='k')
    ax3.plot(X3, Y, color='k')
    ax1.set_xlim(xlims1)
    ax2.set_xlim(xlims2)
    ax3.set_xlim(xlims3)
    ax1.set_xlabel(X1lab)
    ax2.set_xlabel(X2lab)
    ax3.set_xlabel(X3lab)
    
    # Set y axis characteristics
    ax1.set_ylim(0., ymax)
    if ylab is None:
        ax1.set_ylabel('Altitude (km)')
    else:
        ax1.set_ylabel(ylab)
    ax1.vlines(0,0,ymax)
    ax2.vlines(0,0,ymax)
    ax1.xaxis.grid(color='k', ls=':', lw=1)
    ax2.xaxis.grid(color='k', ls=':', lw=1)
    ax3.xaxis.grid(color='k', ls=':', lw=1)
    ax1.yaxis.grid(color='k', ls=':', lw=.5)
    ax2.yaxis.grid(color='k', ls=':', lw=.5)
    ax3.yaxis.grid(color='k', ls=':', lw=.5)
    # Add label to above plot unless location specified
    if label:
        if labelpos is None:
            labelpos = (.15,.91)
        plt.figtext(labelpos[0], labelpos[1], labelTx)

    return fig, ax1, ax2, ax3
#======================================================================
def oplot_3panel(X1, X2, X3, Y, ax1=None, ax2=None, ax3=None, alf=1.):
    """Add additional differential reflectivity and copolar correlation
    coefficient lines to existing graphic."""
    ax1.plot(X1, Y, color='k', alpha=alf)
    ax2.plot(X2, Y, color='k', alpha=alf)
    ax3.plot(X3, Y, color='k', alpha=alf)
#======================================================================