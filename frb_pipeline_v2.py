#!/usr/bin/env/python
'''


# MAX TO DO lIST
 * update .rsp file generator
 * add function to fitt to models, and models as arguments
 * fix LC errror, add data downloader
 * update LC plots
 * add command line arguments
 * update / check relative paths 
 * write documentation

#=============================================================================#
# NAME:                                                                       #
# INPUTS:                                                                     #
# OUTPUTS:                                                                    #
# PURPOSE:                                                                    #
#                                                                             #
# LAST MODIFIED: 20-August-2022 by Maxwell A. Fine                            #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2022 Maxwell A. Fine and CHIME collaboration                  #
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


How to use:
-----------


Notes / Next Steps / WIP
------------------------
    * Need to change method of producing outdir, and moving files from
      subproccess to system
    * Need to update relative paths for producing outdir, and moving files
    * There should be a comparsion for computer time of making N lightcurves
      with HEASoft vs rebinning a single LC's data
    * The LC plots could be redone to be more usefull
    * Functions to produce sky images / DPI and search for sources left in
    script but not used in current version
        * search search_out_catalog needs to be finished
    * need new SNR cutoff for fluence limit?
    * add siwft event data downloder function
    * update docstring for light_curve_analysis to update search method
    * update doc strings for lc_plot and beyond
    * make executable
    * fix error with .rps
    * make easy results file, snr cutoff
    * try changing the names of dirs
    * make time_window a key word argument
    * update search widows to allow for different number of search times
    * make runnable from the command line with commands
    * write up documenation for inputs and general structure



Known Bugs
----------

#=============================================================================#
'''

# imports
import numpy as np
import matplotlib as pyplot
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib as mpl
import os, sys, time, shutil
import astropy.io.fits as pyfits
import numpy as np
import re
import subprocess # for running shell commands
import sys
from tqdm import tqdm
from swifttools.swift_too import Clock, ObsQuery, VisQuery
from datetime import datetime,timedelta
import json
from astropy.stats import poisson_conf_interval


#-----------------------------------------------------------------------------#
# heasoft functions
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

    hk_dir = datadir + '/' + swift_id + '/bat/hk'

    subprocess.run(['cp ' +str(hk_dir) + '/sw' + swift_id + 'bdecb.hk.gz '
                    + str(outdir) + '/sw' + swift_id + 'bdecb.hk.gz'  ],
                    shell=True, cwd = '../data', stdout = stdout,
                    stderr =stderr)

    subprocess.run(['cp data/' +str(swift_id) + '/auxil/sw' + swift_id
                    + 'sat.fits.gz '+ 'data/' + str(outdir) + '/sw'
                    + swift_id + 'sat.fits.gz' ],
                    shell=True, cwd = '..', stdout = stdout, stderr =stdout)

    subprocess.run(['cp data/' +str(swift_id) + '/bat/event/' + str(evt_file)
                + ' ' + 'data/' + str(outdir) + '/' + str(evt_file) ],
                    shell=True, cwd = '..', stdout = stdout, stderr =stderr)
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

def heasoft_lc(evt_file, energy_bans, timdel, cwd, stdout=None, stderr=None,
                         clobber = True):
    ''' Creates lightcurve of specified SWIFT/BAT evt_file using its assocaited
        heasoft_mask

        Parameters
        ----------
        evt_file: str
            name of SWIFT/BAT event file for target

        energy_bans: str
            binning for Energy bands for the light curve, In KeV, seperated by
            both - and , ex energy_bans='15-25,25-50,50-100,100-150'

        timedel: float
            time resolution to produce light curve

        outdir: str
            out directory where results & working files are placed for target
            additonaly, dir  where HEASoft commands are executed

        datadir: str
            directory which has the data directories for the SWIFT targets

        stdout:
            STDOUT for logging

        stderr:
            STDERR for logging

        clobber: True/False
            if True, overwrites files

        '''

    subprocess.run(['batbinevt detmask=frb.mask ' + str(evt_file)
                    + ' timedel = ' + str(timdel)  + ' weighted = no outtype ='
                    + ' lc energybins = ' + energy_bans
                    + ' outfile = frb.lc clobber ='
                    + str(clobber) ],
                    shell=True, cwd=cwd, stdout = stdout, stderr =stderr)
    return

def heasoft_evt_dpi(evt_file, tstart, tstop, cwd,
                      stdout=None, stderr=None, clobber = True):
    '''Runs 'batbinevt' to create event Dectector Plane Image dpi for the
    given evt_file from tstart to tstop. This is the same as heasoft_bgck_dpi
    besides output name

        Parameters
        ----------
        evt_file: str
            name of SWIFT/BAT event file for target

        tstart: float
            start time of DPI, in seconds and in MET

        tstop: float
            end time of DPI, in seconds and in MET

        cwd: str
            dir where HEASoft commands are executed

        datadir: str
            directory which has the data directories for the SWIFT targets

        stdout:
            STDOUT for logging

        stderr:
            STDERR for logging

        clobber: True/False
            if True, overwrites files

        '''


    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = evt.dpi energybins=- timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=dpi detmask=frb.mask weighted=no timedel=0.0'
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True,  )
    return

def heasoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=None, stderr=None, clobber = True):
    '''Runs 'batbinevt' to create event Dectector Plane Image dpi for the
    given evt_file from tstart to tstop. This is the same as heasoft_evt_dpi
    besides the LC output name

        Parameters
        ----------
        evt_file: str
            name of SWIFT/BAT event file for target

        tstart: float
            start time of DPI, in seconds and in MET

        tstop: float
            end time of DPI, in seconds and in MET

        cwd: str
            dir where HEASoft commands are executed

        datadir: str
            directory which has the data directories for the SWIFT targets

        stdout:
            STDOUT for logging

        stderr:
            STDERR for logging

        clobber: True/False
            if True, overwrites files


        '''
    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = bkg.dpi energybins=- timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=dpi detmask=frb.mask weighted=no timedel=0.0'
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True,  )

    return

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


    att = 'sw' + str(swift_id) +'sat.fits.gz'

    aux = 'sw' + str(swift_id) + 'bevtr.fits'

    subprocess.run(['batbinevt '  + str(evt_file[:-3])
                    + ' outfile = frb.pha energybins=CALDB:80 timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=pha detmask=frb.mask weighted=yes timedel=0.0'],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True, )

    subprocess.run(['batphasyserr frb.pha CALDB'],stdout=stdout , stderr=stderr, cwd=cwd , shell=True, )


    subprocess.run(['chmod u=rw,g=rw,o=rw ' + str(evt_file[:-3])],stdout=stdout ,
                   stderr=stderr, cwd=cwd , shell=True,)

    subprocess.run(['gunzip ' + str(evt_file) + ' clobber = True'],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)

    subprocess.run(['batmaskwtevt detmask=frb.mask clobber = True auxfile = auxfile infile=' +str(evt_file[:-3]) + ' attitude = ' +str(att)
                   + ' ra =' + str(ra) + ' dec =' +str(dec)],
                  stdout=stdout , stderr=stderr, cwd=cwd , shell=True,)


    subprocess.run(['batupdatephakw frb.pha auxfile clobber = True'],stdout=stdout ,
                   stderr=stderr, cwd=cwd , shell=True,)

    subprocess.run(['batdrmgen frb.pha frb.rsp clobber= True'],stdout=stdout ,
                   stderr=stderr, cwd=cwd , shell=True,)

    return

def heasoft_sky_image(swift_id, cwd, stdout=None, stderr=None, clobber = True):
    ''' Creates sky image using 'batfftimage' and the evt.dpi and bgck.dpi files

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    cwd: str
        dir where HEASoft commands are executed

    stdout:
        STDOUT for logging

    stderr:
        STDERR for logging

    clobber: True/False
        if True, overwrites files

    '''

    subprocess.run(['batfftimage evt.dpi attitude = sw' + swift_id + 'sat.fits.gz detmask=frb.mask'
                       + ' outfile = frb.sky bkgfile=bkg.dpi clobber =' + str(clobber)],
                       shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return

def heasoft_batcelldetect(sky_image, cwd, incatalog, outcatalog='cat.fits', snrthresh=3.5,
                           stdout=None, stderr=None, clobber = True):
    '''Runs 'batcelldetect' over sky image

    Parameters
    ----------
    sky_image : str
        name of sky image file to run batcelldetect over

    cwd: str
        dir where HEASoft commands are executed

    incatalog: str
        input catalog file for batcelldetect

    outcatalog: str
        name for output catalog file

    snrthresh: float
        snr cutoff for batcelldetect, lowest is 3.5

    stdout:
        STDOUT for logging

    stderr:
        STDERR for logging

    clobber: True/False
        if True, overwrites files

    '''

    subprocess.run(['batcelldetect ' + str(sky_image) + ' snrthresh = ' + str(snrthresh)
                   + ' incatalog = ' + str(incatalog) +' outfile=' + str(outcatalog)],
                  shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return


def wrapper_bgck_heasoft(swift_id, outdir, evt_file, tstart,
                         tstop, datadir, clobber = 'True', stdout=None, stderr = None):
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

    move_ess_heasoft_files(swift_id, evt_file, outdir, datadir, stdout=None, stderr = None)

    cwd=outdir

    heasoft_mask(swift_id, evt_file, outdir, stdout=stdout, stderr=None, clobber = clobber)

    heasoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    return


def wrapper_evt_detect_heasoft(swift_id, evt_file, tstart, tstop,
                               incatalog='test_cat', outcatalog = 'out.cat',
                               cwd='None', snrthresh='3.5', stdout=None,
                               stderr=None, clobber = True):

    '''Wrapper function: makes evt_dpi, sky image , and runs batcelldetect

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    evt_file: str
        name of SWIFT/BAT event file for target

    tstart: float
        start time of DPI, in seconds and in MET

    tstop: float
        end time of DPI, in seconds and in MET

    incatalog: str
        input catalog file for batcelldetect

    outcatalog: str
        name for output catalog file

    cwd: str
        dir where HEASoft commands are executed

    stdout:
        STDOUT for logging

    stderr:
        STDERR for logging

    clobber: True/False
        if True, overwrites files

    '''

    heasoft_evt_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    heasoft_sky_image(swift_id, cwd, stdout=stdout, stderr=stderr, clobber = clobber)

    heasoft_batcelldetect('frb.sky', cwd, incatalog, outcatalog, snrthresh=snrthresh,
                           stdout=stdout, stderr=stdout, clobber = clobber)
    return
#-----------------------------------------------------------------------------#
# non-heasoft functions

def run_search(incatalog, outdir, datadir, energy_bans, timedel,
               clobber='True', snrthresh='3.5',
               swift_evt_file_ending='bevshpo_uf.evt.gz'):
    '''Runs the function search_target on all targets in incatalog

    Parameters
    ----------
    incatalog: str
        input catalog file for batcelldetect

    outdir: str
        out directory where results & working files are placed for target
        additonaly, dir where HEASoft commands are executed

    swift_id : str
        SWIFT observation ID of target

    evt_file: str
        name of SWIFT/BAT event file for target

    tstart: float
        start time of DPI, in seconds and in MET

    tstop: float
        end time of DPI, in seconds and in MET

    datadir: str
        directory which has the data directories for the SWIFT targets

    energy_bans: str
        binning for Energy bands for the light curve, In KeV, seperated by
        both - and , ex energy_bans='15-25,25-50,50-100,100-150'

    timedel: float
        time resolution to produce light curve

    clobber: True/False
        if True, overwrites files

    snrthresh: float
        snr cutoff for batcelldetect, lowest is 3.5

    swift_evt_file_ending: str
        file ending for SWIFT event file

    outcatalog: str
        name for output catalog file

    cwd: str
        dir where HEASoft commands are executed

    stdout:
        STDOUT for logging

    stderr:
        STDERR for logging

    '''

    swift_id, ra, dec, ra_err, dec_err, swift_trig_time, chime_id = read_in_catalog_file(incatalog)

   # log('running search')
    for i in tqdm(range(len(swift_id)),desc='search over '+ str(len(swift_id)) + ' targets'):
        #i += 4
        # need to add
        # download swift/bat data for target if doesnt already exist

        resultsdir = outdir + '/' + swift_id[i]
        bck_tstart = str(float(swift_trig_time[i]) - 70)
        bck_tstop = str(float(swift_trig_time[i]) - 20)

        search_target(swift_id[i], ra[i], dec[i], energy_bans, timedel,
                             incatalog, resultsdir, datadir,
                             swift_trig_time[i], ra_err[i], dec_err[i],
                             bck_tstart, bck_tstop, snrthresh, outdir)

    return

def search_target(swift_id, ra, dec, energy_bans, timdel, incatalog, outdir, datadir,
                         trig_time, ra_error, dec_error, bck_tstart,
                         bck_tstop, snrthresh, clobber = 'True', swift_evt_file_ending='bevshpo_uf.evt.gz'):
    ''' Moved relevent SWIFT/BAT files into outdir for target, produces a
        lightcurve using HEASoft analysis software, searching the light_curve
        for a PEAK in SNR within the specified time window over logmarthicly
        spaced time scales. Uses resulting time scale, and peak to determine a
        fluence limit for the target. Outputs results into a .json file

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    ra: float
        ra of SWIFT target, in degrees

    dec: float
        dec of SWIFt targt, in degrees

    energy_bans: str
        binning for Energy bands for the light curve, In KeV, seperated by
        both - and , ex energy_bans='15-25,25-50,50-100,100-150'

    timedel: float
        time resolution to produce light curve

    incatalog: str
        input catalog file for batcelldetect

    outdir: str
        out directory where results & working files are placed for target
        additonaly, dir where HEASoft commands are executed

    datadir: str
        directory which has the data directories for the SWIFT targets

    trig_time: float
        trigger time of SWIFT target in MET, used for centering of the search
        window

    ra_error: float
        error in RA, degrees, currently not used

    dec_error: float
        error in dec, degrees, currently not used

    bck_tstart: float
        start time of background DPI, in seconds and in MET

    bck_tstop: float
        end time of background DPI, in seconds and in MET

    snrthresh: float
        snr cutoff for batcelldetect, lowest is 3.5

    clobber: True/False
        if True, overwrites files

    swift_evt_file_ending: str
        file ending for SWIFT event file

    '''

    results = outdir[0:-(len(swift_id) + 1)]
    create_heasoft_outdir(swift_id, outdir=results, stdout=None, stderr=None)

    with open(outdir + '/' + 'log.txt','wt',) as f:
        stdout = f
        stderr = f

        evt_file = get_evt_file(swift_id, swift_evt_file_ending='bevshpo_uf.evt.gz')

        wrapper_bgck_heasoft(swift_id, outdir, evt_file, bck_tstart,
                         bck_tstop, datadir, stdout=stdout, stderr = stderr)

        valid_evt = check_valid_evt_file(outdir, evt_file, 10, trig_time)
        #valid_evt = True
        if valid_evt == True:
            #print('data is valid')
            # check for existance of data around chime trigger
            #evt_file, bck_tstart, bck_tstop = inspect_evt_file()

            # move files into outdir, produce mask, bck_dpi
            wrapper_bgck_heasoft(swift_id, outdir, evt_file, bck_tstart,
                             bck_tstop, datadir, stdout=stdout, stderr = stderr)

            heasoft_rsp(evt_file, swift_id, cwd=outdir, ra=ra, dec=dec, stdout=stdout, stderr=stderr, clobber = clobber)

            # might want to move catalog file into outdir for convience?

            # produce lc
            heasoft_lc(evt_file, energy_bans, timdel,
                        cwd=outdir, stdout=stdout, stderr=stderr, clobber = clobber)

            # plots lc
            time_window='10'
            lc_file = 'frb.lc'
        #light_curve_analysis(swift_id, trig_time,
        #    time_window, lc_file, energy_bans, outdir)

            evt_start = trig_time
            evt_tstop = trig_time + 1
            # progress bar
            #for i in tqdm(range(len([evt_start])),leave=False,
            #              desc='searching ' + str(len([evt_start])) +' time windows'):
#
                # produces evt_dpi, and sky_image, looks for sources
            #    wrapper_evt_detect_heasoft(swift_id, evt_file, evt_start, evt_tstop,
            #                       incatalog=incatalog, outcatalog = 'out.cat',
            #                       cwd=outdir, snrthresh=snrthresh, stdout=stdout,
            ##                       stderr=stderr, clobber = clobber)
#

                # searches outcatalog from batcelldetect for sources found near target, if found writes to file / break
                # at the moment just prints if it found something
                #search_out_catalog(outdir=outdir, outcatalog='out.cat', RA_target=ra,
                        #           DEC_target=dec, radius_2=5**2)

            # if batcelldetect does not detect a source IE low SNR some function to establish flux limit?
            # make noise image (maybe make this as a default output?)
            return

        else:
         print('No data around trig_time')
         return


def  check_valid_evt_file(outdir, evt_file, time_window, trig_time):
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
        data = hdul[1].data
        #ra = header_dict['RA_PNT']
        #dec = header_dict['DEC_PNT']
        time_start = data[0][0]
        time = hdul[1].data["Time"]
    #print(np.min(time))
    #print(np.max(time))
    if trig_time <= time[0]:
        print('bad time window, trig time too low')
        return False
    if trig_time >= time[-1]:
        print('bad time window trig time too high')
        return False
    else:
        return True




def light_curve_analysis(swift_id, trig_time, time_window, lc_file, energy_bans, outdir):
    '''Searches for SNR peaks in the light curve in the search window by computing a
    running average with logscale time-window sizes
    returns SNR_peak_value, SNR_peak_time (compared to trig time), LC_analysis

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    trig_time: float
        trigger time of SWIFT target in MET, used for centering of the search
        window

    time_window: float
        window centered at trig_time that is searched, in seconds

    lc_file: str
        name of light curve file

    energy_bans: str
        binning for Energy bands for the light curve, In KeV, seperated by
        both - and , ex energy_bans='15-25,25-50,50-100,100-150'

    outdir: str
        out directory where results & working files are placed for target
        additonaly, dir where HEASoft commands are executed

    '''

    # this works but is kinda slow
    # is there a faster way?
    print('opening LC')
    with fits.open(outdir + '/' + lc_file) as hdul:
        #hdul.info()
        header_dict = hdul[0].header
        time = hdul[1].data["Time"]
        rate = hdul[1].data['rate']
        error = hdul[1].data['error']
        totcounts = hdul[1].data['totcounts']
        # fracexp = hdul[1].data['FRACEXP']
        # columns = hdul[1].data.columns
    #print (trig_time)
    #print(time[0])
    #print(time[-1])
    #if trig_time <= time[0]:
    #    print('bad time window, trig time too low')
    #    return
    #if trig_time >= time[-1]:
    #    print('bad time window trig time too high')
    #    return
    bans = energy_bans.split(',')
    #bans


    search_time_res = np.logspace(-4, 0, 5) # not user specified atm
    #s, trigger time +/- this is the search window

    # calcualte base SNR as running average fails for same time_res
    n = len(rate[0])
    peak_snr_arr = np.empty((len(search_time_res),n))
    peak_snr_arr_old = np.copy(peak_snr_arr)
    peak_snr_arr_new = np.copy(peak_snr_arr)
    time_res = time[1]-time[0]
    rate_dict = {}
    snr_dict = {}
    time_dict = {}
    results_dict = {}
    time_peak_arr = np.empty((len(search_time_res),n))

    for k in tqdm(range(len(search_time_res)),leave=False,
                  desc='searching ' + str(len(search_time_res)) +' time windows'):
        time_ratio = (time_res / search_time_res[k])
        new_bins = time_ratio * len(time)
        snr_dict[search_time_res[k]] = {}
        time_dict[search_time_res[k]] = {}
        rate_dict[search_time_res[k]] = {}
        # improve window_size code, this could be made to be more efficient
        window_size = time_ratio
        if np.floor(np.log10(search_time_res[k])) == -4:
            window_size = 1
        for i in range(n):
            new_rate_i = moving_average(rate[:,i], int(window_size**-1) )
            # new_time might be incorrect! - off by the egde case of moving_average
            new_time = np.linspace(time[0], time[-1], num=len(new_rate_i))
            mean = np.mean(new_rate_i)
            std = np.std(new_rate_i)
            snr = (new_rate_i - mean) / std

            rate_dict[search_time_res[k]][bans[i]] = new_rate_i
            snr_dict[search_time_res[k]][bans[i]] = snr
            time_dict[search_time_res[k]][bans[i]] = new_time


            #print('start window time + ' + str(new_time[(np.abs(new_time - (trig_time) - float(time_window) )).argmin()]))
            #print('end window time ' + str(new_time[(np.abs(new_time - (trig_time) + float(time_window) )).argmin()]))


            #peak_snr_arr_new[k, i] = np.max(snr[start_i:stop_i])

            #print(len(snr[start_i:stop_i]))

            trig_i = (np.abs(new_time - (trig_time)).argmin())
            time_ban = time[1]-time[0]
            index_ban = int(5 / time_ban)
            i_start = trig_i - index_ban
            i_stop = trig_i + index_ban


            peak_snr_arr[k,i] = np.max(snr[i_start:i_stop])



            idx_peak = (np.abs(peak_snr_arr[k,i] - snr[i_start:i_stop]).argmin())
            idx_peak +=i_start


            # idx peak must be wrong!
            time_peak_arr[k,i] = new_time[idx_peak]


    #print('time tests')
    #print(-1*new_time[trig_i] + new_time[i_stop])
    #print(new_time[trig_i] - new_time[i_start])
    #print(new_time[trig_i] - new_time[0])
    peak_snr = np.max(peak_snr_arr[:,0])
    # if that works then k, and I must be wrong?
    k = np.where(peak_snr == peak_snr_arr[:,0])[0][0]
    i = np.where(peak_snr == peak_snr_arr)[0][0]
    idx = (np.abs(snr[i_start:i_stop] - peak_snr).argmin())
    time_peak = time_peak_arr[k,i]


    #time_peak_snr = time_dict[search_time_res[k]][bans[i]][idx]

    results_dict['time_peak_snr'] = time_peak
    results_dict['peak_snr'] = peak_snr
    results_dict['peak_time_scale'] = search_time_res[k]
    results_dict['time_offset'] = time_peak - trig_time
    #results_dict['SNR_dict'] = snr_dict
    #results_dict['time_dict'] = time_dict
    results_dict['bans'] = bans


    # determine upper limit on counts, and background
    time_ratio = (time_res / search_time_res[k])
    new_bins = time_ratio * len(time)
    slice_size = int(len(totcounts) / len(time_dict[search_time_res[k]][bans[i]]))
    lens = (np.ones(int(new_bins), dtype = 'int64')) * int(time_ratio**-1)
    index = np.hstack(([0], lens[:-1])).cumsum()
    time_scale_counts = np.add.reduceat(totcounts, index, axis=0,)

    mean_counts = np.mean(time_scale_counts)
    std_counts = np.std(time_scale_counts)
#'time_peak_snr
    # add user imput for interval style
    count_limit = poisson_conf_interval(time_scale_counts, interval='root-n',
                                        sigma=1, background=0,
                                        confidence_level=None)

    new_time = np.linspace(time[0], time[-1], len(time_scale_counts))

    # select limit based on SNR peak

    count_idx = (np.abs(new_time - time_peak)).argmin()
    count_rate_limit = count_limit[1,count_idx]
    results_dict['count_rate_limit'] = count_rate_limit
    results_dict['background_counts'] = mean_counts
    results_dict['std_counts'] = std_counts

    # write results to file
    with open(outdir + '/result.json', 'w') as fp:
        json.dump(results_dict, fp) # this is



    plot_lc(swift_id, trig_time, rate_dict, snr_dict,
    search_time_res, time_dict, peak_snr,
    peak_snr_arr, bans, outdir, time_peak)


    return #results_dict, peak_snr_arr_new, peak_snr_arr

def plot_lc(swift_id, trig_time, rate_dict, snr_dict,
         search_time_res, time_dict, peak_snr,
         peak_snr_arr, bans, outdir, time_peak):
         '''Produces plots of the lc_file for a swift target at different time
         scales and energy bans

         Parameters
         ----------
         swift_id : str
             SWIFT observation ID of target

         trig_time: float
             trigger time of SWIFT target in MET, used for centering of the
             search window

         rate_dict: dict
             dictionary containing rate/s counts for each enery band and at the
             varring time scales

         snr_dict: dict
            dictionary containing snr arrs for each enery band and at the
            varring time scales

         search_time_res: arr
             numpy array containing the different time scales of the search

         time_dict: dict
            dictionary containing time arrs for each enery band and at the
            varring time scales

         peak_snr: float
            peak SNR value from the rate based SNR search

         peak_snr_arr: numpy arr
            numpy array containing the different SNR peaks for each
            time scales of the search

         bans:  # check

         time_peak: float
            time of SNR peak in MET

         outdir: str
             out directory where results & working files are placed for target
             additonaly, dir where HEASoft commands are executed'''

         mpl.rcParams['agg.path.chunksize'] = 10000
         plt.rcParams.update({'font.size': 20})
         n = len(bans)
         fig, ax = plt.subplots(nrows=n, ncols=len(search_time_res), sharex=True, figsize=(16*4, 6*5))
         k = np.where(peak_snr == peak_snr_arr[:,0])[0][0]
         # rate plot
         for k in range(len(search_time_res)):
             for i in range(len(bans)):
                 time = time_dict[search_time_res[k]][bans[i]]
                 rate_i = rate_dict[search_time_res[k]][bans[i]]
                 ax[i,k].plot(time-time[0], rate_i, label=bans[i] +' KeV')
                 trig_time_plt = find_nearest(time, trig_time)
                 ax[i,k].axvline(x=trig_time_plt-time[0], color ='black', linestyle=':', label='Trigger Time')
                 ax[i,k].set_ylabel('Rate ($s^{-1}$)')
                 ax[i,k].set_xlabel('Time ($s$)')
                 legend = ax[i,k].legend(bbox_to_anchor =(1, 1))#,loc='upper right')
                 legend.set_title(r'Time Resolution '+ str(search_time_res[k]) + '$s$')
         fig.suptitle('Light Curve for ' + str(swift_id) )
         fig.tight_layout()
         plt.savefig(outdir + '/lc_plot_' + str(swift_id)+ '.pdf')
         plt.close()

         # SNR plot
         plt.rcParams.update({'font.size': 20})
         fig, ax = plt.subplots(nrows=n, ncols=len(search_time_res), sharex=True, figsize=(16*4, 6*5))
         for k in range(len(search_time_res)):
             for i in range(len(bans)):
                 time = time_dict[search_time_res[k]][bans[i]]
                 snr_i = snr_dict[search_time_res[k]][bans[i]]
                 #plt.plot(time-time[0], snr_i, label=bans[i] +' KeV')
                 ax[i,k].plot(time-time[0], snr_i, label=bans[i] +' KeV')
                 trig_time_plt = find_nearest(time, trig_time)
                 ax[i,k].axvline(x=trig_time_plt-time[0], color ='black', linestyle=':', label='Trigger Time')

                 ax[i,k].set_ylabel('Rate ($s^{-1}$)')
                 ax[i,k].set_xlabel('Time ($s$)')


                 ax[i,k].set_ylabel('SNR')
                 ax[i,k].set_xlabel('Time ($s$)')
                #trig_time_plt = find_nearest(time, trig_time)
                #plt.axvline(x=trig_time_plt-time[0], color ='black', linestyle=':', label='Trigger Time')

                # this assumes peak is in the last index for i
                 snr_peak_time = find_element(snr_i, time, np.max(snr_i))
                 ax[i,k].axvline(x=time_peak-time[0], color ='purple',
                             linestyle='-.', label='Peak SNR Time (in window)')
               # legend = plt.legend(bbox_to_anchor =(1, 1))#,loc='upper right')
                 legend = ax[i,k].legend(bbox_to_anchor =(1, 1))#,loc='upper right')

                 legend.set_title(r'SNR '+ str(search_time_res[k]) + '$s$')
         fig.suptitle('SNR  ' + str(swift_id) )
         fig.tight_layout() # dies here
         plt.savefig(outdir + '/SNR_plot_' + str(swift_id)+ '.pdf')
         plt.close()
         return


def moving_average(x, w):
     '''Caculates a running average of width w, of array x

     Parameters
     ----------

     x: arr
        numpy arr that calculate the running average of arr x

     w: int
        size of window for running average

     '''
     return np.convolve(x, np.ones(w), 'valid') / w


def find_nearest(array, value):
    '''returns the nearest value in an array

    Parameters
    ----------

    array: arr

    value:

    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def find_element(arr1, arr2, value):
    'returns value of arr2 coresponding to closest value in arr1'
    arr1 = np.asarray(arr1)
    arr2 = np.asarray(arr2)
    idx = (np.abs(arr1 - value)).argmin()
    return arr2[idx]


def create_heasoft_outdir(swift_id, outdir, stdout=None, stderr=None):
    '''Creates outdir for heasoft outputs to be placed in

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    outdir: str
        out directory where results & working files are placed for target
        additonaly, dir  where HEASoft commands are executed

    stdout:
        STDOUT for logging

    stderr:
        STDERR for logging

    '''

    subprocess.run(['mkdir -p ' + str(swift_id)], cwd=outdir, shell=True, stdout=stdout, stderr=stderr)
    return


def read_in_catalog_file(cat_file):
    '''Reads input catalog file and returns arrarys of collums. Expected input
    is a .fits file with the following data collums:

    swift_id, ra, dec, ra_err, dec_err, swift_trig_time, chime_id

    Collums:
    --------

    swift_id : str
        SWIFT observation ID of target

    ra: float
        ra of target, degrees

    dec: float
        dec of target, degrees

    ra_err: flaot
        error of ra, degrees

    dec_err: flaot
        error of dec, degrees

    swift_trig_time: float
        trigger time of target in SWIFT MET

    chime_id: str,
        cross refrence ID

    '''
    with fits.open(cat_file) as hdul:
            #hdul.info()
            #header_dict = hdul[0].header

            data = hdul[1].data

    swift_id = []
    ra = np.zeros(len(data))
    dec = np.zeros(len(data))
    ra_err = np.zeros(len(data))
    dec_err = np.zeros(len(data))
    swift_trig_time = np.zeros(len(data))
    chime_id = []
    for i in range(len(data)):
        swift_id.append(str(data[i][0]))
        ra[i] = data[i][1]
        dec[i] = data[i][2]
        ra_err[i] = data[i][3]
        dec_err[i] = data[i][4]
        swift_trig_time[i] = data[i][5]
        chime_id.append((data[i][6]))


    # string as to not loose leading zeros
    swift_id = np.asarray(swift_id)
    # while this var is called chime_id its really the cross refrence name
    chime_id = np.asarray(chime_id)
    return swift_id, ra, dec, ra_err, dec_err, swift_trig_time, chime_id

# need to remake / not finished
def search_out_catalog(outdir, outcatalog, RA_target, DEC_target, radius_2=2**2, evt_window=1, id='test'):
    ''' This function is not complete, and for code not currently implemented'''
    with fits.open(str(outdir + '/' + outcatalog)) as hdul:
       # hdul.info()
        header_dict = hdul[0].header
        data = hdul[1].data
        #ra = header_dict['RA_PNT']
        #dec = header_dict['DEC_PNT']
    target_results_dict = {}

    for i in range(len(data)):
        RA_source = data[i][11]
        DEC_source = data[i][12]
        SNR_source = data[i][33]

        if distance_from_target_2(RA_target, DEC_target, RA_source, DEC_source) <= radius_2:
            if SNR_source>= snr_cutoff:
            # record detections close to target
                target_results_dict[str(i)] = {}

                target_results_dict[str(i)]['RA'] = data[i][11]
                target_results_dict[str(i)]['DEC'] = data[i][12]
                target_results_dict[str(i)]['SNR'] = data[i][33]
                #print('FOUND TARGET!')
                print('RA ' + str(data[i][11]) + ' DEC '  + str(data[i][12]) + ' SNR ' + str(data[i][33]) +' SWIFT ID ' + str(id) + ' Window_size ' +str(evt_window) )
                break

def distance_from_target_2(RA_target, DEC_target, RA_source, DEC_source):
        '''cauculates the distance of a target from a found source'''
        return np.abs((RA_target-RA_source)**2 + (DEC_target-DEC_source)**2)



def get_evt_file(swift_id, swift_evt_file_ending='bevshpo_uf.evt.gz'):
    '''Returns string of filename for SWIFT/bat event file based on its swift_id

    Parameters
    ----------
    swift_id : str
        SWIFT observation ID of target

    swift_evt_file_ending: str
        file ending for SWIFT event file

    '''
    evt_file = 'sw' + swift_id + swift_evt_file_ending
    return evt_file

#-----------------------------------------------------------------------------#
def main():
    import argparse
    """
    Runs the main program to be executed after if_name_ = main
    """
    # Help string to be shown using the -h option
    descStr = """
    to be written later"""

    epilog_text="""
    to be written later
    """
    # Parse the command line options
    parser = argparse.ArgumentParser(description=descStr,epilog=epilog_text,
                                 formatter_class=argparse.RawTextHelpFormatter)

    # add arguments like this
    #parser.add_argument("-U", dest="units", type=str, default="Jy/beam",
        #                help="Intensity units of the data. [Jy/beam]")
    #units          = args.units
    args = parser.parse_args()
    outdir='../results'
    datadir = '../data'
    #timedel = 1
    timedel = '1e-4'#'1e-4'
    incatalog = 'test_cat'
    energy_bans = '15-25,25-50,50-100,100-150'
    swift_evt_file_ending='bevshpo_uf.evt.gz'
    clobber = True
    snrthresh = '6'
    log=print
    time_window='10'
    lc_file = 'frb.lc'
    snr_cut_off = '7'

    run_search(incatalog, outdir, datadir, energy_bans, timedel,
               clobber=clobber, snrthresh=snrthresh,
               swift_evt_file_ending=swift_evt_file_ending)
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    main()
