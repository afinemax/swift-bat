 #!/usr/bin/env/python
#=============================================================================#
# NAME:               swift_chime_pipeline                                    #
# INPUTS:         see -H                                                      #
# OUTPUTS:   output json, and diagnostic plots                                #
# PURPOSE:     searches for gamm-rays from CHIME/FRBs  with swift/bat         #
#                                                                             #
# LAST MODIFIED: April 2023 by Maxwell A. Fine                                #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2023 Maxwell A. Fine and CHIME collaboration                  #
#                                                                             #
# Contact:                                                                    #
#     Maxwell Fine: https://afinemax.github.io/afinemax1/                     #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
#=============================================================================#


# imports
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib as mpl
import os, sys, time, shutil
import astropy.io.fits as pyfits
import re
import subprocess # for running shell commands
import sys
from tqdm import tqdm
from swifttools.swift_too import Data
from swifttools.swift_too import Clock, ObsQuery, VisQuery
from datetime import datetime,timedelta
import json
from astropy.stats import poisson_conf_interval
import xspec # for spectrum models
from scipy.signal import boxcar
from scipy import signal
from swifttools.swift_too import Data
from swifttools.swift_too import ObsQuery
from os import listdir
from os.path import isfile, join
import json
from scipy import special
from matplotlib.gridspec import GridSpec
import argparse
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.gridspec as gridspec


def make_outdir(outdir):
    '''makes outdir dir for results dirs to go into'''
    os.makedirs(outdir, exist_ok = True)
'''
expects a .csv file with the followiing
swift_id, trig_time, chime_id
'''


def heasoft_rsp(evt_file, swift_id, cwd, ra, dec, stdout=None, stderr=None,
                    clobber = True):
    '''Runs a series of HEASoft commands to produce the .rsp file for the
    SWIFT/BAT target
    Parameters
    ----------
    evt_file: str
        name of SWIFT/BAT event file for target
    swift_id : str
        SWIFT observation ID of target
    cwd: str
        dir where HEASoft commands are executed
    ra: float
        ra of SWIFT target, in degrees
    dec: float
        dec of SWIFt targt, in degrees
    stdout:
        STDOUT for logging
    stderr:
        STDERR for logging
    clobber: True/False
        if True, overwrites files
    '''

    evt_file = 'sw' + str(swift_id) + evt_file
    #print('making .rsp file')


    # need to move sat file into the dir
    # i dont use the aux file
    # cwd dir is the events dir
    att = 'sw' + str(swift_id) +'sat.fits.gz'

    #aux = 'sw' + str(swift_id) + 'bevtr.fits'



    # make aux file
    # need the 'aux' file made by batmaskwtevt
    #  batmaskwtevt,

    subprocess.run(['cp ' + str(evt_file) + ' evt_file'],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)

    subprocess.run([' gzip -d -f ' + str(evt_file) ],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)

    subprocess.run(['cp ' + ' evt_file '+ str(evt_file) ],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)

    new_evt = evt_file[:-3]
    evt_file = new_evt


    subprocess.run(['batmaskwtevt detmask=frb.mask clobber = True auxfile = auxfile infile=' +str(evt_file) + ' attitude = ' +str(att)
                   + ' ra =' + str(ra) + ' dec =' +str(dec)],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)
   # print('made aux file!')



   #  I think I need to decrompress the evt file?
    # how to decompress?
    #

  #  print('uncompressed file name ', new_evt)


    # make pha file
    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = frb.pha energybins=CALDB:80 timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=pha detmask=frb.mask weighted=yes timedel=0.0'],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True, )
  #  print('made frb.pha file')

    # from the documentation
    # adds the BAT systematic error vector to a BAT spectral file.
    # The systematic error vector contains an estimate of the fractional
    # systematic error in each BAT spectral channel.
    # Normally this information comes from the BAT calibration database (CALDB).
   # print('error below?')
    subprocess.run(['batphasyserr frb.pha CALDB'],stdout=stdout , stderr=stderr, cwd=cwd , shell=True, )
   # print('updated frb.pha file!')

    # give read write acess for frb.pha
   # subprocess.run(['chmod u=rw,g=rw,o=rw ' + 'frb.pha'],stdout=stdout ,
   #                stderr=stderr, cwd=cwd , shell=True,)


    # produce a mask weights collum in the evt file
    # this makes it a targeted file
    # I thought we wanted to make it untargeted?
    # ask aaron

    subprocess.run(['batmaskwtevt detmask=frb.mask clobber = True auxfile = auxfile infile=' +str(evt_file) + ' attitude = ' +str(att)
                   + ' ra =' + str(ra) + ' dec =' +str(dec)],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)
   # print('produce mask weights check')

    # from documentation:
    # batupdatephakw updates the BAT ray tracing columns in a mask weighted spectrum.
    # This is important especially for gamma-ray bursts or other events where the spacecraft has slewed.
    # Without the update, the keywords will reflect the ray tracing status at the end of the event
    # data set, which will lead to an erroneous detector response matrix.

    subprocess.run(['batupdatephakw frb.pha auxfile clobber = True'],stdout=stdout ,
                   stderr=stderr, cwd=cwd , shell=True,)
   # print('batupdatephakw check')

    # finally make the matrix!
    subprocess.run(['batdrmgen frb.pha frb.rsp clobber= True'],stdout=stdout ,
                   stderr=stderr, cwd=cwd , shell=True,)


    return

# functions
xspec.XspecSettings.chatter = 0

def fit_spectrum(countlim_at_burst,cwd):

    fluence_limits = []

    # make inputs!
    model_strs = ["tbabs(bb)","tbabs(bb)","tbabs(cutoffpl)","tbabs(cutoffpl)","tbabs(po)","tbabs(po)"]
    model_pars = [[0.1,10.,1.],[10.0,10.,1.],[0.1,0.5,500.,1.],[10.0,0.5,500.,1.],[0.1,2.0,1.],[10.,2.0,1.]]
    mod_fluence_lims,mod_avg_e, plot_vals = [],[],[]


    for i,(model_str,model_par) in enumerate(zip(model_strs,model_pars)):
        model = xspec.Model(model_str)
        model.setPars(*model_par)
        fluence_lim, weighted_E = get_fluence_limit(countlim_at_burst, cwd, model)
        fluence_limits.append(fluence_lim)
    return fluence_limits, weighted_E


def get_fluence_limit(countlim_at_burst, cwd,model):

    fsp_file = cwd  + 'frb.rsp'

    with fits.open(fsp_file) as hdul:
        hdul.info()
        header_dict = hdul[0].header
        data = hdul[1].data
        ENERG_LO = hdul[1].data["ENERG_LO"]
        ENERG_HI = hdul[1].data['ENERG_HI']
        N_GRP = hdul[1].data['N_GRP']
        F_CHAN = hdul[1].data['F_CHAN']
        N_CHAN = hdul[1].data['N_CHAN']
        MATRIX = hdul[1].data['MATRIX']
    energy_lims = np.linspace(15,150,1000) #kev
    xspec.AllModels.setEnergies('15 150 1000') # limits of swift sky images

    energies = ENERG_LO[(ENERG_LO > energy_lims[0] ) & (ENERG_LO  < energy_lims[-1])]
    mod_energies = np.array(model.energies(0),dtype='float')[:-1]

    weighted_E = np.sum(mod_energies* model.values(0)) / np.sum(model.values(0)) # model weighted photon E

    resp = 1 # dummy variable, RSPs have allraedy been normalises to their respective area
    avg_area = np.sum(resp*model.values(0)) / np.sum(model.values(0))
    KeV_in_erg = 1.60218e-9
    fluence_lim = countlim_at_burst / (avg_area / weighted_E / KeV_in_erg)

    return fluence_lim,  weighted_E


def wrapper_fluence_limit(evt_file, swift_id, cwd, ra, dec, sky_image, w, bkg_image,
                        stdout=None, stderr=None, clobber=True):
    #'makes .rsp file', calculates fluence limit



    # make .rsp file
    # need to move .sat file!
    # cwd is the evt dir

    try:
        heasoft_rsp(evt_file, swift_id, cwd, ra, dec, stdout=None, stderr=None,
                    clobber = True)

    except Exception as e:
        error_msg += str(e) + ','
        #fluence_lim = None


    # determine coutlim at burst from sky image
    # assuming 1 pixel accuracy
    data = sky_image
    #countlim_at_burst = 0.03170389 # fake number
    sky = SkyCoord(ra=ra, dec=dec, unit=u.degree)
    x, y = w.world_to_pixel(sky)





    # check for valid values of x, y in pixel coords

    if np.isnan(x):
        return None, None
    if np.isnan(y):
        return None, None
    if x < 0:
        return None, None
    if y< 0 :
        return None, None

    evt_counts = data[int(x),int(y)] # this could still be <0
    bkg_counts = bkg_image[int(x),int(y)]




    if evt_counts > 0:
        upper_limit = evt_counts + np.sqrt(bkg_counts) # 1 sigma limit

    if evt_counts < 0:
        upper_limit = np.sqrt(bkg_counts) # 1 sigma limit

    if evt_counts<0:
        evt_counts = 0


    lower_limit, upper_limit = poisson_conf_interval(int(evt_counts), interval='kraft-burrows-nousek',
    confidence_level=0.95, background=bkg_counts)




#    lower_limit, upper_limit = poisson_conf_interval(evt_counts, interval='kraft-burrows-nousek', background=bkg_counts)
#    lower_limit, upper_limit = poisson_conf_interval(evt_counts)

    if np.isnan(upper_limit):
        upper_limit =0
    if upper_limit <=0:
        upper_limit =0
    countlim_at_burst = upper_limit
    fluence_limit, weighted_e = fit_spectrum(countlim_at_burst,cwd,)

    # fluence limit is a list of arrays!

    fluence_lim_list = [arr[0] for arr in fluence_limit]


    #print('fluence limit resultssss')
    #print(np.shape(fluence_limit), np.shape(upper_limit))
    #print(np.dtype(fluence_limit), np.dtype(upper_limit))
    #print(fluence_limit, upper_limit)

    return fluence_lim_list, upper_limit


def read_in_catalog(incatalog):
    data = np.genfromtxt(str(incatalog), dtype=str, delimiter=',')
    swift_ids = data[:,0]
    trig_time = data[:,1]
    chime_ids = data[:,2]
    ra = data[:,3]
    dec = data[:,4]

    # convert trig time into swift time
    # there must be a faster way of doing this
    #swift_time_arr = np.zeros(len(trig_time))
    #for i in range(len(trig_time)):


    return swift_ids, trig_time, chime_ids, ra ,dec


def create_heasoft_outdir(swift_id, cwd, stdout=None, stderr=None):
    '''Creates outdir for outputs of pipeline to be placed in

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    cwd : str
        Directory tree where command should be executed

    stdout:
        STDOUT for logging

    stderr: STDERR for logging
        '''

    subprocess.run(['mkdir -p ' + str(swift_id)], cwd= cwd, shell=True,
                        stdout=stdout, stderr=stderr)

    return

def move_ess_heasoft_files(swift_id, evt_file, outdir, datadir,
                            stdout=None, stderr = None):
    '''Moves necessary SWIFT/BAT files from SWIFT/BAT data directory into outdir
    for convience in generating heasoft outputs

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    evt_file: str
        name of SWIFT/BAT event file for target

    outdir: str
        out directory where results & working files are placed for target

    datadir: str
        directory which has the data directories for the SWIFT targets

    stdout:
        STDOUT for logging

    stderr: STDERR for logging

    '''

    swift_id = str(swift_id)
    hk_dir = datadir + '/' + swift_id+ '/bat/hk'

    cwd = datadir

    evt_file = 'sw' + swift_id + evt_file

    datadir = datadir + '/' + swift_id
    outdir = outdir + '/' + swift_id

    sat_file = 'sw' + swift_id + 'sat.fits.gz'

    #sw00096133003bevshpo_uf.evt.gz

    # datadir = ../data
    # outdir = ../results

    #/bigger dir/ # this file

    #bigger dir = ..


    subprocess.run(['cp ' + cwd + '/' +str(swift_id) + '/bat/event/' + str(evt_file)
                + ' ' + str(outdir) + '/' + str(evt_file) ],
                    shell=True, cwd = cwd, stdout = stdout, stderr =stderr)

    subprocess.run(['cp ' + cwd + '/' +str(swift_id) + '/auxil/' + str(sat_file)
                + ' ' +  str(swift_id) + '/bat/event/' ],
                    shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return



def heasoft_mask(swift_id, evt_file, outdir, stdout=None,
                    stderr=None, clobber = True):
    '''Procduces SWIFT/bat detector mask by running 'bathotpix
        and 'batbinevt'

        Parameters
        ----------
        swift_id : str
            SWIFT observation ID of target

        evt_file: str
            name of SWIFT/BAT event file for target

        outdir: str
            out directory where results & working files are placed for target
            additonaly, dir  where HEASoft commands are executed

        stdout:
            STDOUT for logging

        stderr:
            STDERR for logging

        clobber: True/False
            if True, overwrites files

        '''

    subprocess.run(['batbinevt ' + str(evt_file) + ' weighted=no outunits=counts'
                    + ' outtype=dpi energy=- clobber=' + str(clobber)
                    + ' outfile= frb.dpi.mask'], cwd = outdir, shell=True,
                    stdout = stdout, stderr =stderr)

    subprocess.run(['bathotpix detmask=sw' +  swift_id + 'bdecb.hk.gz outfile ='
                    +' frb.mask infile= frb.dpi.mask clobber =' + str(clobber) ],
                    shell=True, cwd=outdir, stdout = stdout, stderr =stderr)
    return


def wrapper_bgck_heasoft(swift_id, outdir, evt_file, tstart,
                         tstop, datadir, clobber = 'True',
                         stdout=None, stderr = None):
    '''Wrapper function: moves files, produces mask , and produces background dpi

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    outdir: str
        out directory where results & working files are placed for target
        additonaly, dir where HEASoft commands are executed

    evt_file: str
        name of SWIFT/BAT event file for target

    tstart: float
        start time of DPI, in seconds and in MET

    tstop: float
        end time of DPI, in seconds and in MET

    datadir: str
        directory which has the data directories for the SWIFT targets

    stdout:
        STDOUT for logging

    stderr:
        STDERR for logging

    clobber: True/False
        if True, overwrites files

    '''

    move_ess_heasoft_files(swift_id, evt_file, outdir, datadir,
                            stdout=None, stderr = None)

    cwd=outdir

    heasoft_mask(swift_id, evt_file, outdir, stdout=stdout,
                                stderr=None, clobber = clobber)

    heasoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    return



def get_swift_bat_data(swift_ids, datadir, overwrite=False):
    '''Downloads Swift data corresponding to the swift_ids and places
        data files in datadir

        Parameters
        ----------
        swift_ids: arr
                   swift target IDs to download

        datadir: str
                relative path to directory to store swift data

        overwrite: True / False
                   if true, overwrites swift data'''

    for i in range(len(swift_ids)):
        data = Data()
        data.obsid = str(swift_ids[i]) # query data for these observations,
                                       # note id is a str
        data.bat =    True # fetch data for bat

        data.submit() # checks for existence of data
        data.clobber = overwrite
        data.download(outdir=datadir)


def  check_valid_evt_file(outdir, evt_file, trig_time, ra, dec):
    ''' Checks that trig_time is within the time_start, time_stop values of the
    Swift/bat observation & that the taret is in the FOV

    Parameters
    ----------
    outdir: str
        out directory where results & working files are placed for target
        additonaly, dir where HEASoft commands are executed

    evt_file: str
        name of SWIFT/BAT event file for target

    time_window: float
        window centered at trig_time that is searched, in seconds

    trig_time: float
        trigger time of SWIFT target in MET, used for centering of the search
        window

    ra: float
        ra of SWIFT target, in degrees

    dec: float
        dec of SWIFt targt, in degrees

    '''
    #trig_time = trig_time
    with fits.open(outdir + '/' + evt_file) as hdul:
        #hdul.info()
        header_dict = hdul[0].header
        #print(header_dict)
        ra_pnt = header_dict['RA_PNT'] # degrees
        dec_pnt = header_dict['DEC_PNT'] # degrees
        data = hdul[1].data
        #ra = header_dict['RA_PNT']
        #dec = header_dict['DEC_PNT']
        time_start = data[0][0]
        time = hdul[1].data["Time"]
    #print(np.min(time))
    #print(np.max(time))
    distance= np.sqrt((ra_pnt-ra)**2 + (dec_pnt-dec)**2)
    if trig_time <= time[0]:
        print('bad time window, trig time too low')
        return False
    if trig_time >= time[-1]:
        print('bad time window trig time too high')
        return False
    if time[-1] - time[0] <= 190:
        print('bad time window, not full 200 seconds')
        return False

    if distance >= 40:
        # bat is ~60 deg by ~30deg so ~40 will almost always be in the FOV
        print('Target more than 40 degrees away from pointing')
        return False
    else:
        return True


def get_lightcurve(outdir, evt_file, time_resolution, energy_bans, lc_file,
                   mask_no_data=True, mask_noise_spikes=True,  mask_cosmic_ray=True):
    '''generates a lightcurve from a .evt swift/bat file'''

    # open lc file, read out data

    with fits.open(outdir + '/' + lc_file) as hdul:
        header_dict = hdul[0].header
    #hdul.info()

        data =hdul[1].data
        time = hdul[1].data["Time"]
        det_id = hdul[1].data['DET_ID']
        evts = hdul[1].data['EVENT_FLAGS']
        pha = hdul[1].data['PHA']
        #mask_w = hdul[1].data['MASK_WEIGHT']
        detx = hdul[1].data['DETX']
        dety = hdul[1].data['DETY']
        pi = hdul[1].data['PI']
        energy = hdul[1].data['ENERGY']

        #rate = np.transpose(rate)
        names = hdul[1].data.columns

   # test#

    times_old = np.copy(time)


    if mask_cosmic_ray == True:
        # mask cosmic ray hits
        # this is based on the nitrates paper

        # "light curve with bins of 50 µs is made using events
        # with energy > 50 keV. Any bin with > 40 counts
        # and > 10 times the average counts at ± 1 s is flagged
        # as a cosmic ray hit."
        cosmic_dt = 1e-5
        # only take > 50 KeV
        bins_cosmic=int(np.ptp(time[energy>40])/cosmic_dt)
        counts_cosmic, bins_cosmic = np.histogram(time[energy>40], bins=bins_cosmic,)
        time_cosmic = np.asarray([(bins_cosmic[i]+bins_cosmic[i+1])/2. for i in range(len(bins_cosmic)-1)])

        mask1 = [counts_cosmic>20] # too few counts still! whaa
        cosmic_snr = (counts_cosmic - np.mean(counts_cosmic) / np.std(counts_cosmic))[mask1]
        mask2 = [cosmic_snr>10]

        # loop through and chop out of time, and energy arr
        bad_times = time_cosmic[mask1][mask2]

        for i in range(len(bad_times)):
            t_bad_start = np.argmin(np.abs(time-bad_times[i] + 1e-4)) # idx values
            t_bad_end = np.argmin(np.abs(time-bad_times[i] - 1e-4))

            if t_bad_end == t_bad_start:
                t_bad_end +=3
            # cut out bad indexes

            time = np.concatenate((time[:t_bad_start], time[t_bad_end:] ))
            energy = np.concatenate((energy[:t_bad_start], energy[t_bad_end:] ))


    # now we make the LC
    bins=int(np.ptp(time)/time_resolution)
    total_counts, bins = np.histogram(time, bins=bins,)
    new_time = np.asarray([(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)])

    if mask_no_data == True:
        # cut out time bins (at time_resolution) with no counts
        # this is intended to cut out regions of no data
        # not reconmended to use at fine time resolutions
        mask_no_data_arr = np.nonzero(total_counts > 0)
        new_time = new_time[mask_no_data_arr]
        total_counts = total_counts[mask_no_data_arr]


    bans = energy_bans.split(',')

    rate = []
    bins=int(np.ptp(time)/time_resolution)
    for i in range(len(bans)):
        lower_e, upper_e = bans[i].split('-')
        lower_e, upper_e = float(lower_e), float(upper_e)

        good_index1 = np.where(energy >= float(lower_e))

        new_energy = energy[good_index1]
        bin_time1 = time[good_index1]

        good_index2 = np.where(new_energy <= float(upper_e))
        bin_time2 = bin_time1[good_index2]

        ban_e = new_energy[good_index2]

        # now we have bin_time2 which is the times of detections of photons in the ban


        rate_counts_i, bins = np.histogram(bin_time2, bins=bins)
        if mask_no_data == True:
            rate_counts_i = rate_counts_i[mask_no_data_arr]
        rate.append(rate_counts_i)


    return new_time, total_counts, np.transpose(np.array(rate)), time, energy,


def get_sky_image(swift_id, chime_id, ra, dec, outdir, datadir,
                  time_0, trig_time, time_window, evt_file, ra_err=0, dec_err=0, clobber='True'):
    '''Computes event dpi, and background dpi, produces sky image for target'''

    # produce paths to evt, hk, auxil

    datadir = '../data_lc' + '/' + swift_id

    # evt files
    #evt_file = datadir +'/bat/event/sw' + swift_id + evt_file
    # Get the absolute path to the "data_lc" directory
    evt_file = 'sw' + swift_id + evt_file #00014920013bevshpo_uf.evt.gz'
    datadir = os.path.abspath(datadir )
    evt_dir = os.path.abspath(datadir + '/bat/event/')

    #sw' +  swift_id + 'bdecb.hk.gz
   # obsid/auxil/swNNNNNNNNNNNsat.fits.gz
    #
    sat_file = os.path.abspath(datadir + '/auxil/sw' + swift_id + 'sat.fits.gz')

    # sw00014920013bdecb.hk.gz
    # hk file
  #  00014920013/bat/hk

    #qualilty map: obsid/bat/hk/swNNNNNNNNNNNbcbdq.hk.gz
    hk_file = os.path.abspath(datadir + '/bat/hk/sw' + swift_id + 'bdecb.hk.gz')

    f = open("output.txt", "w")

    stderr = f
    stdout = f

    # produce  dpi mask, and mask
    subprocess.run(['batbinevt ' + str(evt_file) + ' weighted=no outunits=counts'
                    + ' outtype=dpi energy=- clobber=' + str(clobber)
                    + ' outfile= frb.dpi.mask'], cwd = evt_dir, shell=True,
                    stdout = stdout, stderr =stderr)


    subprocess.run(['bathotpix detmask=' + str(hk_file) + ' outfile ='
                    +' frb.mask infile= frb.dpi.mask clobber =' + str(clobber) ],
                    shell=True, cwd=evt_dir, stdout = stdout, stderr =stderr)

    # make background dpi
    tstart = time_0
    tstop = trig_time  - (2*time_window)
    bgk_duration = tstop-tstart


    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = ' + 'bkg.frb.dpi' + ' energybins=14-150 timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=dpi detmask=frb.mask weighted=no ecol=ENERGY timedel=0.0 outunits=RATE '
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=evt_dir , shell=True,  )

    # make evt.dpi
    tstart =trig_time - time_window
    tstop = trig_time + time_window
    evt_duration = tstop-tstart



    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = ' + 'evt.frb.dpi' + ' energybins=14-150 timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=dpi detmask=frb.mask weighted=no ecol=ENERGY timedel=0.0 outunits=RATE '
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=f , stderr=f, cwd=evt_dir , shell=True,  )

    # make sky image
    subprocess.run(['batfftimage infile=evt.frb.dpi outfile=evt.frb.sky.img attitude=' + str(sat_file)
                    + ' bkgfile = bkg.frb.dpi detmask=frb.mask clobber=YES' ],
                       stdout=f , stderr=f, cwd=evt_dir , shell=True,  )


    # make noise image
    subprocess.run(['batcelldetect infile=evt.frb.sky.img outfile=batcelldetect.out snrthresh=6 bkgvarmap=bkg.frb.sky.img clobber=YES' ],
                       stdout=f , stderr=f, cwd=evt_dir , shell=True,  )



    #bkg.frb.dpi
    # make time
    '''
    subprocess.run(['batfftimage infile=evt.frb.dpi outfile=evt.frb.sky.img attitude=' + str(sat_file)
                    + ' detmask=frb.mask clobber=YES' ],
                       stdout=f , stderr=f, cwd=evt_dir , shell=True,  )

    subprocess.run(['batfftimage infile=bkg.frb.dpi outfile= bkg.frb.sky.img attitude=' + str(sat_file)
                    + ' detmask=frb.mask clobber=YES' ],
                       stdout=f , stderr=f, cwd=evt_dir , shell=True,  )
        '''
#    print(evit_dir)


    # move sky image into outdir

    # open sky image data and return image data + wcs
    sky_image = os.path.abspath(datadir + '/bat/event/evt.frb.sky.img')
    with fits.open(sky_image) as hdul:

        #header_dict = hdul[0].header
        data = hdul[0].data
        w = WCS(hdul[0].header,)#key='T')

    bkg_sky_image = os.path.abspath(datadir + '/bat/event/bkg.frb.sky.img')
    with fits.open(bkg_sky_image) as hdul:

        #header_dict = hdul[0].header
        bkg_data = hdul[0].data
        #w = WCS(hdul[0].header,)#key='T')

    return data, w, bkg_data #(bkg_data *evt_duration / bgk_duration)


def get_lightcurve(outdir, evt_file, time_resolution, energy_bans, lc_file,
                   mask_no_data=True, mask_noise_spikes=True,  mask_cosmic_ray=True):
    '''generates a lightcurve from a .evt swift/bat file'''

    # open lc file, read out data

    with fits.open(outdir + '/' + lc_file) as hdul:
        header_dict = hdul[0].header
    #hdul.info()

        data =hdul[1].data
        time = hdul[1].data["Time"]
        det_id = hdul[1].data['DET_ID']
        evts = hdul[1].data['EVENT_FLAGS']
        pha = hdul[1].data['PHA']
        #mask_w = hdul[1].data['MASK_WEIGHT']
        detx = hdul[1].data['DETX']
        dety = hdul[1].data['DETY']
        pi = hdul[1].data['PI']
        energy = hdul[1].data['ENERGY']

        #rate = np.transpose(rate)


        names = hdul[1].data.columns

    if np.ptp(time) > 250:
         return None

   # test#

    times_old = np.copy(time)








    if mask_cosmic_ray == True:
        # mask cosmic ray hits
        # this is based on the nitrates paper

        # "light curve with bins of 50 µs is made using events
        # with energy > 50 keV. Any bin with > 40 counts
        # and > 10 times the average counts at ± 1 s is flagged
        # as a cosmic ray hit."
        cosmic_dt = 1e-5
        # only take > 50 KeV
        bins_cosmic=int(np.ptp(time[energy>40])/cosmic_dt)
        counts_cosmic, bins_cosmic = np.histogram(time[energy>40], bins=bins_cosmic,)
        time_cosmic = np.asarray([(bins_cosmic[i]+bins_cosmic[i+1])/2. for i in range(len(bins_cosmic)-1)])

        mask1 = [counts_cosmic>20] # too few counts still! whaa
        cosmic_snr = (counts_cosmic - np.mean(counts_cosmic) / np.std(counts_cosmic))[mask1]
        mask2 = [cosmic_snr>10]

        # loop through and chop out of time, and energy arr
        bad_times = time_cosmic[mask1][mask2]

        for i in range(len(bad_times)):
            t_bad_start = np.argmin(np.abs(time-bad_times[i] + 1e-4)) # idx values
            t_bad_end = np.argmin(np.abs(time-bad_times[i] - 1e-4))

            if t_bad_end == t_bad_start:
                t_bad_end +=1
            # cut out bad indexes

            time = np.concatenate((time[:t_bad_start], time[t_bad_end:] ))
            energy = np.concatenate((energy[:t_bad_start], energy[t_bad_end:] ))


    if mask_noise_spikes == True:
        noise_dt = 0.016
        bins_low=int(np.ptp(time)/noise_dt)

        counts_low, bins_low = np.histogram(time[energy<15], bins=bins_low,)
        time_low = np.asarray([(bins_low[i]+bins_low[i+1])/2. for i in range(len(bins_low)-1)])


        #
        bins_low=int(np.ptp(time)/noise_dt)
         # I should use the same bins!
        bins_low=int(np.ptp(time)/noise_dt) #
        counts_high, bins_low = np.histogram(time[energy>50], bins=bins_low,)
        time_high = np.asarray([(bins_low[i]+bins_low[i+1])/2. for i in range(len(bins_low)-1)])

        low_snr = (counts_low - np.mean(counts_low)) / np.std(counts_low)
        high_snr = (counts_high - np.mean(counts_high)) / np.std(counts_high)

        mask1 = time_high[low_snr>20] # too few counts still! whaa
        mask2 = time_high[high_snr<2.5]
        # lets assume it is the same index number
        mask = np.intersect1d(mask1, mask2)
        for i in range(len(mask)):
            t_bad_start = np.argmin(np.abs(time-mask[i] + 1e-4)) # idx values
            t_bad_end = np.argmin(np.abs(time-mask[i] - 1e-4))

            if t_bad_end == t_bad_start:
                t_bad_end +=1

            time = np.concatenate((time[:t_bad_start], time[t_bad_end:] ))
            energy = np.concatenate((energy[:t_bad_start], energy[t_bad_end:] ))


    # now we make the LC

    bins=int(np.ptp(time)/time_resolution)
    total_counts, bins = np.histogram(time, bins=bins,)
    new_time = np.asarray([(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)])

    if mask_no_data == True:
        # cut out time bins (at time_resolution) with no counts
        # this is intended to cut out regions of no data
        # not reconmended to use at fine time resolutions
        mask_no_data_arr = np.nonzero(total_counts > 0)
        new_time = new_time[mask_no_data_arr]
        total_counts = total_counts[mask_no_data_arr]

    bans = energy_bans.split(',')

    rate = []
    bins=int(np.ptp(time)/time_resolution)

    for i in range(len(bans)):
        lower_e, upper_e = bans[i].split('-')
        lower_e, upper_e = float(lower_e), float(upper_e)

        good_index1 = np.where(energy >= float(lower_e))

        new_energy = energy[good_index1]
        bin_time1 = time[good_index1]

        good_index2 = np.where(new_energy <= float(upper_e))
        bin_time2 = bin_time1[good_index2]

        ban_e = new_energy[good_index2]

        # now we have bin_time2 which is the times of detections of photons in the ban


        rate_counts_i, bins = np.histogram(bin_time2, bins=bins)
        if mask_no_data == True:
            rate_counts_i = rate_counts_i[mask_no_data_arr]
        rate.append(rate_counts_i)


    return new_time, total_counts, np.transpose(np.array(rate)), time, energy,

# new new lc analysis
# new lc_analysis

# new new lc analysis
# new lc_analysis

def lc_analysis(swift_id, evt_file, trig_time, time_window,  energy_bans, outdir):

    # initalize results dict
    lc_analysis_dict = {}
    swift_id = str(swift_id)

    time_resolution = 1e-2
    lc_file = 'sw' + swift_id + evt_file

    outdir = outdir +'/' + swift_id

    # open lc file, read out data

    lc_return = get_lightcurve(outdir,
                                                                              lc_file,
                                                                              time_resolution,
                                                                              energy_bans, lc_file)

    if lc_return is None:
        return None
    time, totcounts, rate, old_time, energy  = lc_return

    time = np.asarray(time)

    # identify energy bans
    bans = energy_bans.split(',')




    # calculate mean and noise
    totcounts_mean = np.mean(totcounts)
    totcounts_noise = np.std(totcounts)
    totcounts_snr = (totcounts-totcounts_mean) / totcounts_noise

    lc_analysis_dict['totcounts_dict'] = {}
    lc_analysis_dict['totcounts_dict']['snr_arr'] = totcounts_snr

    # indexs to cut SNR curves to search window size
    start_i = np.where(np.abs(time - trig_time + time_window) == np.min(np.abs(time-trig_time + time_window)))[0][0]
    stop_i = np.where(np.abs(time - trig_time - time_window) == np.min(np.abs(time-trig_time - time_window)))[0][0]


    # also in rate domain
    rate_snr_dict = {}
    lc_analysis_dict['bans_dict'] = {}
    for i in range(len(bans)):
        # double check that bans is not in reverse order
        mean_i = np.mean(rate[:,i])
        noise_i = np.std(rate[:,i])
        rate_snr_dict[bans[i]] = (rate[:,i] - mean_i ) / noise_i
        lc_analysis_dict['bans_dict'][bans[i]] = {}


    # boxcar search over SNR with varying sized boxcars
    # we want to do a sliding box car search each energy ban + the total counts arr
    # first lets move into SNR space

    # -2, 0, 3
    time_res = np.logspace(-2, np.log10(2*time_window), 13) #, 0.01, 0.1, 1
    for i in tqdm(range(len(time_res)),leave=False,
                  desc='searching ' + str(len(time_res)) +' time windows'):


        time_size_per = int(np.abs((time[-1] - time[0]))) / len(time) # time size per index
        window_len = int(time_res[i] / time_size_per)
        window = np.ones(window_len,dtype='float32')/ np.sqrt(window_len)

        if len(window) == 0:
            window = [1]
        #print('stop-start', int(stop_i-start_i))
        #print('swift_id', swift_id)
        # bad one looks short



        # normalize by width of boxcar

        totcounts_convolve = signal.convolve(window, totcounts_snr[start_i: stop_i])




        totcounts_snr_max = np.max(totcounts_convolve) # peak from the boxcar search


        lc_analysis_dict['totcounts_dict'][time_res[i]] = totcounts_snr_max

        # now for each rate space SNR
        for j in range(len(bans)):


            bans_convolve_j = signal.convolve(window, rate_snr_dict[bans[j]][start_i: stop_i])




            bans_snr_max_j = np.max(bans_convolve_j) # peak from the boxcar search


            lc_analysis_dict['bans_dict'][bans[j]][time_res[i]] = bans_snr_max_j


    # returns dicts containg SNR, and max SNR for the boxcar search per time scale

    return lc_analysis_dict, rate_snr_dict, time, energy, old_time


def expected_snr(n_samples):

    snr = special.erfinv((-1/n_samples) + 1) * np.sqrt(2)

    # this sometimes returns negative values
    # we only care about the + ones


    return np.abs(snr) # return postive value only


def lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, trig_time,
                energy_bans, time_window, outdir, chime_id, swift_id, old_time, sky_image, w, ra, dec):
    '''produces plots of the SNR for the best time scale, and plots SNR vs time scale '''


    # determine best SNR
    # I need to filter out the low SNR energy channel prior to this
    # will be done in LC analysiss or based on the energy bands
    # plot

    # indexs to cut SNR curves to search window size
    #start_i = np.where(np.isclose(time - trig_time + time_window, 0, rtol=1e-05) == True)[0][0]
    #stop_i = np.where(np.isclose(time - trig_time - time_window, 0, rtol=1e-05) == True)[0][0]

    start_i = np.where(np.abs(time - trig_time + time_window) == np.min(np.abs(time-trig_time + time_window)))[0][0]
    stop_i = np.where(np.abs(time - trig_time - time_window) == np.min(np.abs(time-trig_time - time_window)))[0][0]



    # new plot
    # we want SNR vs time scale,
    # LC and zoomed in LC
    # if we can print off peak time scale and snr
    # maybe print swift id and chime id
    #

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(24, 9),#constrained_layout=True,
        gridspec_kw={'width_ratios': [1, 1, 2]})

    #fig = plt.figure(constrained_layout=True, figsize=(20, 30))

    #gs = GridSpec(3, 3, figure=fig)
    fs =15
    plt.yticks(fontsize=fs, color='k')
    plt.xticks(fontsize=fs, color='k')


    # top left
#    ax1 = fig.add_subplot(gs[0, 0])
    counts_hist, bins_hist, etc = axs[0,0].hist(energy[energy>0], bins=int(500), color='green')

    axs[0, 0].set_title('Energy Distrubution Total', fontsize=fs, color='k')
    axs[0, 0].set_xlabel('KeV', fontsize=fs, color='k')
    axs[0, 0].set_ylabel('N', fontsize =fs, color='k')
    axs[0, 0].set_xlim(0, 350) # look into xlim

    # middle

    e_idx_start = np.argmin(np.abs(old_time -trig_time +time_window))
    e_idx_end = np.argmin(np.abs(old_time -trig_time  - time_window))

    counts_hist, bins_hist, etc = axs[0,1].hist(energy[e_idx_start:e_idx_end][energy[e_idx_start:e_idx_end]>0], bins=int(500), color='green')

    axs[0, 1].set_title('Energy Distrubution Search Window', fontsize=fs, color='k')
    axs[0, 1].set_xlabel('KeV', fontsize=fs, color='k')
    axs[0, 1].set_ylabel('N', fontsize =fs, color='k')
    axs[0, 1].set_xlim(0, 350) # look into xlim
#    axs[0,0].set_ylim(0, np.max(counts_hist)*1.05)

    # bottom left
    #axs[1, 0].plot(x, -y, 'tab:green')
    #axs[1, 0].set_title('Axis [1, 0]')
    # lets loop through and grap all SNR vs time scale

    #ax2 = fig.add_subplot(gs[1, 0])
    time_res = []
    snr_val = []
    my_keys = list(lc_analysis_dict['totcounts_dict'].keys())


    for i in range(len(my_keys) -1):
        i += 1
        time_res.append(float(my_keys[i]))
        snr_val.append(lc_analysis_dict['totcounts_dict'][my_keys[i]])
    axs[1, 0].scatter(time_res, snr_val, label='Total Counts')

    peak_time_scale = time_res[np.argmax(snr_val)]
    peak_snr = np.max(snr_val)
    peak_snr_arr = [peak_snr]
    time_scales_arr = time_res

    out_snr = []



    # now loop through energy bans
    bans = list(lc_analysis_dict['bans_dict'].keys())
    for i in range(len(bans)):
        energy_ban = bans[i]
        time_res_keys = list(lc_analysis_dict['bans_dict'][energy_ban].keys())

        time_res_list = []
        snr_val_list = []
        for j in range(len(time_res_keys)):
            time_res_list.append(float(time_res_keys[j]))
            snr_val_list.append(lc_analysis_dict['bans_dict'][energy_ban][time_res_keys[j]]) # changed to j from i
            # if statment for getting max SNR, and peak time scale



        axs[1, 0].scatter(time_res_list, snr_val_list, label=energy_ban + ' KeV')
        out_snr.append(snr_val_list)

        if np.max(snr_val) > peak_snr:
            peak_snr = np.max(snr_val)
            peak_time_scale = time_res(np.argmax(snr_val))

    time_resolution = peak_time_scale
    out_time_res = time_res_list
    out_bans = bans
    # out_snr



    # I should just be able to rerun LC analysis?
    # now wait this should still be faster

    #plt.legend(facecolor='white', framealpha=1)
    axs[1, 0].set_xlabel('Time Scale (s)', fontsize=fs, color='k')
    axs[1, 0].set_ylabel('Peak SNR', fontsize =fs, color='k')
    #axs[1, 0].set_yticks(fontsize=fs, color='k')
    #axs[1, 0].set_xticks(fontsize=fs, color='k')
    axs[1, 0].set_title('Time Scale vs Peak SNR', fontsize= 15)
    axs[1, 0].set_xscale('log')

    timescales = np.logspace(-2,np.log10(2*6.3), 100)#1.8,1000)
    #timescales = np.linspace(1e-8, 6, 1000)

    snr = np.zeros(len(timescales))
    time_window = 6
    n_samples = (2 * time_window )/ timescales


    snr = expected_snr(n_samples)
    #plt.scatter( timescales,snr)
    #axs[1,0].set_xlim(.7e-1, np.log10(2*time_window))
    axs[1, 0].plot(timescales, snr, ':', color='k', label='Expected Noise')
    axs[1, 0].legend(fontsize=12, labelcolor='k',facecolor='white', ncol=2 )

    #plt.show()
    # right side
    #ax3 = fig.add_subplot(gs[0, 1])
    #ax4 = fig.add_subplot(gs[1, 1])

    #axs[0, 2]

    legend_time_res = 'Time Resolution: '+ str(1e-2) + ' s'

    snr_arr_tot = lc_analysis_dict['totcounts_dict']['snr_arr']

    #fig.suptitle()
    axs[0, 2].set_title('SNR of Lightcurve', fontsize=15, color='k')
    axs[1, 2].set_title('SNR of Search Window', fontsize=15, color='k')

    axs[0, 2].set_xlabel('Time (s)', fontsize=fs, color='k')
    axs[0, 2].set_ylabel('SNR', fontsize=fs, color='k')
    axs[0, 2].step(time - time[0],snr_arr_tot, color='purple')
    axs[0, 2].axvline(x=trig_time - time[0], linestyle='--', color='k',
                    label="Trigger Time")
    axs[0, 2].legend(fontsize=fs-2,title=legend_time_res)

    axs[1, 2].set_xlabel('Time (s)', fontsize=fs, color='k')
    axs[1, 2].set_ylabel('SNR', fontsize=fs, color='k')
    axs[1, 2].step((time - time[0])[start_i:stop_i],
    (snr_arr_tot)[start_i:stop_i], color='purple')
    axs[1, 2].axvline(x=trig_time - time[0], linestyle='--', color='k',
                    label="Trigger Time")
    axs[1, 2].legend(fontsize=fs-2,title=legend_time_res)


    # 1 1

    #axs(projection=w,) #figsize=(15, 15))

    print('made it here!')
    #axs[1,1].remove() # remove to add again with projection
    #ax[1,1] = fig.add_subplot(1, 1, projection=wcs)
    sky_loc = SkyCoord(ra * u.degree, dec* u.degree)
    x, y = w.world_to_pixel(sky_loc)
    is_nan = False

    if np.isnan(x):
        x=0
        is_nan = True
    if np.isnan(y):
        y=0
        is_nan = True


    if is_nan is False:

        #x_win, y_win = w.world_to_pixel(SkyCoord((ra) * u.degree, (dec)* u.degree))
        win = 30
        x_win = win
        y_win = win
        image_data = sky_image

        img = axs[1,1].imshow(image_data[int(y -y_win):int(y +y_win),int(x-x_win):int(x+x_win),] , cmap='Greys_r', origin='lower')

        #axs[1,1].set_colorbar()
        fig.colorbar(img)

        axs[1,1].set_ylabel('DEC', size=15)
        axs[1,1].set_xlabel('RA',size=15)
        axs[1,1].scatter(x_win,y_win, label='Target', color='r')
        axs[1,1].legend()

    #print('made it here!2.0')


    fig.suptitle("SWIFT ID " + swift_id + '     Peak SNR: '+ str(peak_snr)[0:4]
    + '\nCHIME ID ' + chime_id +'    Peak Timescale ' + str(time_resolution)[0:4] +'s', size=fs)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    fig.tight_layout()
    plt.savefig(outdir + '/' + 'swift_id_'+ swift_id +'_chime_id_' + chime_id+'_ diagnostic.pdf', bbox_inches='tight')
    #plt.show()
    plt.close()

    # returns peak snr, and peak time scale
    # want to also return an array of the peak snr, and time scales

    return peak_snr, time_resolution, out_time_res, out_bans, out_snr


def new_lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, trig_time,
                energy_bans, time_window, outdir, chime_id, swift_id, old_time, sky_image, w, ra, dec):
    '''produces plots of the SNR for the best time scale, and plots SNR vs time scale '''

    start_i = np.where(np.abs(time - trig_time + time_window) == np.min(np.abs(time-trig_time + time_window)))[0][0]
    stop_i = np.where(np.abs(time - trig_time - time_window) == np.min(np.abs(time-trig_time - time_window)))[0][0]



    fig = plt.figure(figsize=(24,9))
    gs = fig.add_gridspec(2, 3, width_ratios=[1,1,2])

    fs =15
    #plt.yticks(fontsize=fs, color='k')
    #plt.xticks(fontsize=fs, color='k')


    # top left
    ax1 = fig.add_subplot(gs[0,0])
    counts_hist, bins_hist, etc = ax1.hist(energy[energy>0], bins=int(500), color='green')
    ax1.set_title('Energy Distrubution Total', fontsize=fs, color='k')
    ax1.set_xlabel('KeV', fontsize=fs, color='k')
    ax1.set_ylabel('N', fontsize =fs, color='k')
    ax1.set_xlim(0, 350) # look into xlim


    # middle
    ax2 = fig.add_subplot(gs[0,1])
    e_idx_start = np.argmin(np.abs(old_time -trig_time +time_window))
    e_idx_end = np.argmin(np.abs(old_time -trig_time  - time_window))

    counts_hist, bins_hist, etc = ax2.hist(energy[e_idx_start:e_idx_end][energy[e_idx_start:e_idx_end]>0], bins=int(500), color='green')

    ax2.set_title('Energy Distrubution Search Window', fontsize=fs, color='k')
    ax2.set_xlabel('KeV', fontsize=fs, color='k')
    ax2.set_ylabel('N', fontsize =fs, color='k')
    ax2.set_xlim(0, 350) # look into xlim

    # bottom left
    ax3 = fig.add_subplot(gs[1,0])
    time_res = []
    snr_val = []
    my_keys = list(lc_analysis_dict['totcounts_dict'].keys())


    for i in range(len(my_keys) -1):
        i += 1
        time_res.append(float(my_keys[i]))
        snr_val.append(lc_analysis_dict['totcounts_dict'][my_keys[i]])
    ax3.scatter(time_res, snr_val, label='Total Counts')

    peak_time_scale = time_res[np.argmax(snr_val)]
    peak_snr = np.max(snr_val)
    peak_snr_arr = [peak_snr]
    time_scales_arr = time_res

    out_snr = []



    # now loop through energy bans
    bans = list(lc_analysis_dict['bans_dict'].keys())
    for i in range(len(bans)):
        energy_ban = bans[i]
        time_res_keys = list(lc_analysis_dict['bans_dict'][energy_ban].keys())

        time_res_list = []
        snr_val_list = []
        for j in range(len(time_res_keys)):
            time_res_list.append(float(time_res_keys[j]))
            snr_val_list.append(lc_analysis_dict['bans_dict'][energy_ban][time_res_keys[j]]) # changed to j from i

        ax3.scatter(time_res_list, snr_val_list, label=energy_ban + ' KeV')
        out_snr.append(snr_val_list)

        if np.max(snr_val) > peak_snr:
            peak_snr = np.max(snr_val)
            peak_time_scale = time_res(np.argmax(snr_val))

    time_resolution = peak_time_scale
    out_time_res = time_res_list
    out_bans = bans

    ax3.set_xlabel('Time Scale (s)', fontsize=fs, color='k')
    ax3.set_ylabel('Peak SNR', fontsize =fs, color='k')
    ax3.set_title('Time Scale vs Peak SNR', fontsize= 15)
    ax3.set_xscale('log')
    timescales = np.logspace(-2,np.log10(2*6.3), 100)#1.8,1000)
    snr = np.zeros(len(timescales))
    time_window = 6
    n_samples = (2 * time_window )/ timescales
    snr = expected_snr(n_samples)
    ax3.plot(timescales, snr, ':', color='k', label='Expected Noise')
    ax3.legend(fontsize=12, labelcolor='k',facecolor='white', ncol=2 )

    # light curves
    ax4 = fig.add_subplot(gs[0,2])
    ax5 = fig.add_subplot(gs[1,2])

    legend_time_res = 'Time Resolution: '+ str(1e-2) + ' s'

    snr_arr_tot = lc_analysis_dict['totcounts_dict']['snr_arr']

    #fig.suptitle()
    ax4.set_title('SNR of Lightcurve', fontsize=15, color='k')
    ax5.set_title('SNR of Search Window', fontsize=15, color='k')

    ax4.set_xlabel('Time (s)', fontsize=fs, color='k')
    ax4.set_ylabel('SNR', fontsize=fs, color='k')
    ax4.step(time - time[0],snr_arr_tot, color='purple')
    ax4.axvline(x=trig_time - time[0], linestyle='--', color='k',
                    label="Trigger Time")
    ax4.legend(fontsize=fs-2,title=legend_time_res, ) #color='k')

    ax5.set_xlabel('Time (s)', fontsize=fs, color='k')
    ax5.set_ylabel('SNR', fontsize=fs, color='k')
    ax5.step((time - time[0])[start_i:stop_i],
    (snr_arr_tot)[start_i:stop_i], color='purple')
    ax5.axvline(x=trig_time - time[0], linestyle='--', color='k',
                    label="Trigger Time")
    ax5.legend(fontsize=fs-2,title=legend_time_res, )#color='k')

    # sky image
    ax6 = fig.add_subplot(gs[1,1], projection=w)
    sky_loc = SkyCoord(ra * u.degree, dec* u.degree)
    x, y = w.world_to_pixel(sky_loc)
    is_nan = False


    if np.isnan(x):
        x=0
        is_nan = True
    if np.isnan(y):
        y=0
        is_nan = True


    if is_nan is False:

        win = 4
        x_win = win
        y_win = win
        image_data = sky_image

        sky_loc = SkyCoord(ra * u.degree, dec* u.degree)
        x, y = w.world_to_pixel(sky_loc)



        img = ax6.imshow(image_data[int(y -y_win):int(y +y_win+1),int(x-x_win):int(x+x_win+1),] , cmap='Greys_r', origin='lower')


        cbar = fig.colorbar(img)
        cbar.ax.tick_params(colors='white')
        ax6.tick_params(axis='x', colors='white')
        ax6.tick_params(axis='y', colors='white')


        ax6.set_ylabel('DEC', size=15, color='white')
        ax6.set_xlabel('RA',size=15, color='white')

        # limits
        print('ra, dec')
        print(ra, dec)
        print('pixel locations')
        print(x, y)
        print(x,)

        sky_loc_a = SkyCoord((ra +1/60)  * u.degree, (dec +1/60)* u.degree)
        x_a, y_a = w.world_to_pixel(sky_loc_a)
        print('x_a, y_a')
        print(x_a, y_a)

        sky_loc_b = SkyCoord((ra -1/60)  * u.degree, (dec -1/60)* u.degree)
        x_b, y_b = w.world_to_pixel(sky_loc_b)

        print('x_b, y_b')
        print(x_b, y_b)

        # if x_a < x_b:
        #     ax6.set_xlim(x_a, x_b)
        # else:
        #     ax6.set_xlim(x_b, x_a)
        #
        # if y_a < y_b:
        #     ax6.set_xlim(y_a, y_b)
        # else:
        #     ax6.set_xlim(y_b, y_a)


        #ax6.set_xlim(x-10, ra + 1.5)
        #ax6.set_ylim(y-10, dec-1.5)

        ax6.scatter(x_win,y_win, label='Target', color='r', marker='x', s=3)
        plt.setp(ax6.get_xticklabels(), visible=True)
        ax6.set_frame_on(True)
        ax6.legend()




    fig.suptitle("SWIFT ID " + swift_id + '     Peak SNR: '+ str(peak_snr)[0:4]
    + '\nCHIME ID ' + chime_id +'    Peak Timescale ' + str(time_resolution)[0:4] +'s', size=fs)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    fig.tight_layout()
    plt.savefig(outdir + '/' + 'swift_id_'+ swift_id +'_chime_id_' + chime_id+'_ diagnostic.pdf', bbox_inches='tight')
    #plt.show()
    plt.savefig('ast425_final_sky.png', bbox_inches='tight', dpi=300, transparent=True)
    plt.close()

    # returns peak snr, and peak time scale
    # want to also return an array of the peak snr, and time scales

    return peak_snr, time_resolution, out_time_res, out_bans, out_snr

def outcatalog(outcatalog_file, swift_id, chime_id, timescale, peak_snr, out_time_res, out_bans, out_snr):
    # Check if file is empty
    is_empty = os.stat(outcatalog_file).st_size == 0

    with open(outcatalog_file, 'a', encoding='utf-8') as f:
        if is_empty:
            f.write('[]')

        # Move the file pointer to the second-to-last byte
        f.seek(-2, os.SEEK_END)

        # If the last character is a closing square bracket,
        # move the pointer back one more byte so we don't overwrite it
        if f.read(1) == ']':
            f.seek(-1, os.SEEK_CUR)

        # Add the new JSON object to the end of the file
        json.dump({
            'swift_id': swift_id,
            'chime_id': chime_id,
            'timescale': timescale,
            'peak_snr': peak_snr,
            'out_time_res': out_time_res,
            'out_bans': out_bans,
            'out_snr': out_snr
        }, f)

        # Add a comma and closing square bracket to the end of the file
        f.write('],')
    return


def append_to_json_file(file_path, data):
    with open(file_path, 'a') as f:
        # Check if file is empty
        is_empty = f.tell() == 0

        if not is_empty:
            f.seek(-1, 2)
            f.truncate()

        if is_empty:
            f.write('[')

        if not is_empty:
            f.write(',')

        json.dump(data, f)
        f.write(']')


def main():
    #test

    """
    Start the function to perform LC search if called from the command line.
    """

    # Help string to be shown using the -h option
    descStr = """Inputs yada yada wip"""

    epilog_text=""" outputs yada yada"""

    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,epilog=epilog_text,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", dest="incatalog", default='frb_test_cat_feb_2023.csv', nargs=1,
                        help="Incatalog file to read in. Formated as swift_id, trig_time, chime,id" )
    parser.add_argument("-e", dest="energy_bans",
                        default ='0-15,15-25,25-50,50-100,100-150,200-250,250-350',
                        help="Energy bans to look in. Formated as '0-15,15-25'")
    parser.add_argument("-f", dest="evtfile", action="store_false", default='bevshpo_uf.evt.gz',
                        help="SWIFT/BAT evt file ending to look into.")
    parser.add_argument("-d", dest="swift_data_dir",  default='../data_lc',
                        help="Directory to Check for & Download SWIFT/BAT data. ")
    parser.add_argument("-o", dest="out_dir", default='../lc_pipeline_res',
                        help="Directory to place results into.")
    parser.add_argument("-t", dest="search_time_window", type=float, default=3,
                        help="Plus and minus this value to search around the trigger")
    parser.add_argument("-r", dest="outcatalog_file", type=str, default='test_outcatalog_batcelldetect_1.2.json',
                        help="name of outputed results catalog.")
    parser.add_argument("-l", dest="download_data", type=bool, default='False',
                        help="If True downloads swift/bat data, if false does not")
    args = parser.parse_args()


    incatalog = args.incatalog
    energy_bans = args.energy_bans
    evt_file = args.evtfile
    datadir = args.swift_data_dir
    outdir = args.out_dir
    time_window = args.search_time_window
    outcatalog_file = args.outcatalog_file
    download_data = args.download_data

    swift_ids, trig_time, chime_ids, ra, dec = read_in_catalog(incatalog)

    # rewrite outcatalog so it is empty
    with open(str(outdir + '/' + outcatalog_file), 'w') as file:
            file.writelines('')

    # make outdir
    make_outdir(outdir)
    download_data = False

    if download_data == True:
        while  download_data is True:
            try:
                get_swift_bat_data(swift_ids, datadir, overwrite=False)
                download_data = False
            except:
                 download_data = True


    result_dict = {}
    for i in tqdm(range(len(swift_ids)),leave=False, desc='searching ' + str(len(swift_ids)) +' Targets'):


        error_msg = ''
        try:
            move_ess_heasoft_files(swift_ids[i], evt_file, outdir, datadir,
                                                            stdout=None, stderr = None)

            cc = Clock()
                                        #datetime(2016,12,31,23,59,59)
            cc = Clock(swifttime=str(trig_time[i]))
            swift_time = cc.met
                                        #print(swift_time)

            my_trig_time = swift_time

            analysis_result = lc_analysis(swift_ids[i],
                                    evt_file, my_trig_time, time_window,  energy_bans, outdir)

            if analysis_result is None:
                error_msg += 'lc_analysis result is None ' + ','

                continue

            lc_analysis_dict, rate_snr_dict, time, energy, old_time = analysis_result


            sky_image, w, bkg_data = get_sky_image(swift_ids[i], chime_ids[i], ra[i], dec[i], outdir, datadir,
                              np.min(old_time), my_trig_time, time_window, evt_file, ra_err=0, dec_err=0, clobber='True')

            #sky_image += bkg_data
            peak_snr, time_resolution, out_time_res, out_bans, out_snr = new_lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, my_trig_time,
                                                    energy_bans, time_window, outdir, chime_ids[i], swift_ids[i], old_time, sky_image, w, float(ra[i]), float(dec[i]))

            try:
                cwd = datadir + '/' + str(swift_ids[i]) + '/bat/event' + '/'
                fluence_lim, count_lim = wrapper_fluence_limit(evt_file, swift_ids[i], cwd, ra[i], dec[i], sky_image, w, bkg_data,
                                            stdout=None, stderr=None, clobber=True)
                count_lim = count_lim[0].astype('float')
            except Exception as e:
                error_msg += str(e) + ','
                fluence_lim = None
                count_lim =  None




            peak_snr, time_resolution, out_time_res, out_bans, out_snr = new_lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, my_trig_time,
                                                energy_bans, time_window, outdir, chime_ids[i], swift_ids[i], old_time, sky_image, w, float(ra[i]), float(dec[i]))



            data = {
                'swift_id': swift_ids[i],
                'chime_id': chime_ids[i],
                'timescale': time_resolution,
                'peak_snr': peak_snr,
                'out_time_res': out_time_res,
                'out_bans': out_bans,
                'out_snr': out_snr,
                'fluence_lim': fluence_lim,
                'count_lim': count_lim,
                'error_msgs': error_msg}

            result_dict[str(swift_ids[i])] = data




            with open(outdir + '/' + outcatalog_file, 'w') as f:
                json.dump(result_dict, f)
                print('updated output file!')


        except Exception as e:
            error_msg += str(e) + ','

            data = {
            'swift_id': swift_ids[i],
            'error_msgs': error_msg}

            result_dict[str(swift_ids[i])] = data

            with open(outdir + '/' + outcatalog_file, 'w') as f:
                json.dump(result_dict, f)
                print('updated output file!')

            continue



if __name__ == '__main__':
    main()
