#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2022 Maxwell A. Fine and CHIME collaboration                                        #
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
import matplotlib as pyplot
import matplotlib.pyplot as plt
from astropy.io import fits
import os, sys, time, shutil
import astropy.io.fits as pyfits
import numpy as np
import re
import subprocess # for running shell commands
import sys
from tqdm import tqdm
from swifttools.swift_too import Clock, ObsQuery, VisQuery
from datetime import datetime,timedelta

#-----------------------------------------------------------------------------#
# heasoft functions
# need to add a function to create a noise image
def create_heasoft_outdir(swift_id, cwd, stdout=None, stderr=None):
    '''Creates outdir for outputs to be placed in'''

    subprocess.run(['mkdir -p ' + str(swift_id)], cwd= cwd, shell=True,
                        stdout=stdout, stderr=stderr)
    return

def move_ess_heasoft_files(swift_id, evt_file, outdir,
                            stdout=None, stderr = None):
    '''Moves necessary SWIFT/BAT files into outdir for convience in generating
    heasoft outputs'IMPORTANT relative paths need to be updated / set to
    variables for final script!'''

# IMPORTANT relative paths need to be updated / set to variables for final script
    subprocess.run(['cp ' +str(hk_dir) + '/sw' + swift_id + 'bdecb.hk.gz '
                    + str(results_dir) + '/sw' + swift_id + 'bdecb.hk.gz'  ],
                    shell=True, cwd = '../data', stdout = stdout,
                    stderr =stderr)

    subprocess.run(['cp data/' +str(swift_id) + '/auxil/sw' + swift_id
                    + 'sat.fits.gz '+ 'data/' + str(results_dir) + '/sw'
                    + swift_id + 'sat.fits.gz' ],
                    shell=True, cwd = '..', stdout = stdout, stderr =stdout)

    subprocess.run(['cp data/' +str(swift_id) + '/bat/event/' + str(evt_file) + ' '
                    + 'data/' + str(outdir) + '/' + str(evt_file) ],
                    shell=True, cwd = '..', stdout = stdout, stderr =stderr)
    return


def heafsoft_mask(swift_id, evt_file, cwd, stdout=None,
                    stderr=None, clobber = clobber):
    '''Procduces SWIFT/bat detector mask by running 'bathotpix
        and 'batbinevt' '''

    subprocess.run(['batbinevt ' + str(evt_file) + ' weighted=no outunits=counts'
                    + ' outtype=dpi energy=- clobber=' + str(clobber)
                    + ' outfile= frb.dpi.mask'], cwd = cwd, shell=True,
                    stdout = stdout, stderr =stderr)

    subprocess.run(['bathotpix detmask=sw' +  swift_id + 'bdecb.hk.gz outfile ='
                    +' frb.mask infile= frb.dpi.mask clobber =' + str(clobber) ],
                    shell=True, cwd=cwd, stdout = stdout, stderr =stderr)
    return

def heafsoft_lc(evt_file, energybins, timdel, cwd, stdout=None, stderr=None, clobber = clobber):
    ''' Creates lightcurve of specified SWIFT/BAT evt_file using its assocaited mask

    Inputs:
    evt_file = string, file name of event file
    energybins = string, specified energy binning in KeV
    timdel = specified timeresolution '''

    subprocess.run(['batbinevt detmask=frb.mask ' + str(evt_file)
                    + ' timedel = ' + timdel  + ' weighted = no outtype ='
                    + ' lc energybins =' + energybins
                    + ' outfile = frb.lc clobber ='
                    + str(clobber) ],
                    shell=True, cwd=cwd, stdout = stdout, stderr =stderr)
    return

def heafsoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=None, stderr=None, clobber = clobber):
    '''Runs 'batbinevt' to create background dpi for the given evt_file
        from tstart to tstop'''

    # this is the same as heafsoft_evt_dpi besides output name
    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = bkg.dpi energybins=- timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    +  ' outtype=dpi detmask=frb.mask '
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True,  )
    return

def heafsoft_evt_dpi(evt_file, tstart, tstop, cwd,
                      stdout=None, stderr=None, clobber = clobber):
    '''Runs 'batbinevt' to create event dpi for the given evt_file
    from tstart to tstop'''

    # this is the same as heafsoft_bgck_dpi besides output name
    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = evt.dpi energybins=- timebinalg=u'
                    + ' clobber=' + str(clobber) +  ' outtype=dpi'
                    + ' detmask=frb.mask '
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True,  )
    return

def heafsoft_sky_image(swift_id, cwd, stdout=None, stderr=None, clobber = clobber):
    ''' Creates sky image using 'batfftimage' and the evt.dpi and bgck.dpi files '''

    subprocess.run(['batfftimage evt.dpi attitude = sw' + swift_id + 'sat.fits.gz detmask=frb.mask'
                       + ' outfile = frb.sky bkgfile=bkg.dpi clobber =' + str(clobber)],
                       shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return

def heafsoft_batcelldetect(sky_image, cwd, incatalog, outcatalog='cat.fits', snrthresh='3.5',
                           stdout=None, stderr=None, clobber = clobber):
    '''Runs 'batcelldetect' over sky image'''

    subprocess.run(['batcelldetect ' + str(sky_image) + ' snrthresh = ' + str(snrthresh)
                   + ' incatalog = ' + str(incatalog) +' outfile=' + str(outcatalog)],
                  shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return


def wrapper_bgck_heasoft(swift_id, outdir, evt_file, tstart,
                         tstop, stdout=None, stderr = None):
    '''Wrapper function: moves files, produces mask , and produces background dpi'''

    move_ess_heasoft_files(swift_id, evt_file, outdir, stdout=None, stderr = None)

    cwd=outdir

    heafsoft_mask(swift_id, evt_file, cwd=outdir, stdout=stdout, stderr=None, clobber = clobber)

    heafsoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    return


def wrapper_evt_detect_heasoft(swift_id, evt_file, tstart, tstop,
                               incatalog='test_cat', outcatalog = 'out.cat',
                               cwd=results_dir, snrthresh='3.5', stdout=None,
                               stderr=None, clobber = clobber):

    '''Wrapper function: makes evt_dpi, sky image , and runs batcelldetect'''

    heafsoft_evt_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    heafsoft_sky_image(swift_id, cwd, stdout=stdout, stderr=stderr, clobber = clobber)

    heafsoft_batcelldetect('frb.sky', cwd, incatalog, outcatalog, snrthresh='3.5',
                           stdout=stdout, stderr=stdout, clobber = clobber)
    return
#-----------------------------------------------------------------------------#

# cwd = results_dir
# outdir = results_dir
# timdel = '1'
# need to finish docstring, rebuilding
def run_search_on_target(swift_id, ra, dec, trigger_time, ra_error, dec_error):
    '''Main function runs search on single SWIFT/BAT target '''

    with open(results_dir + '/' + 'log.txt','wt',) as f:
        stdout = f
        stderr = f

        # check for existance of data around chime trigger
        evt_file, bck_tstart, bck_tstop = inspect_evt_file()

        # move files into outdir, produce mask, bck_dpi
        wrapper_bgck_heasoft(swift_id, outdir, evt_file, tstart,
                         tstop, stdout=stdout, stderr = stderr)

        # produce lc
        heafsoft_lc(evt_file, '15-150', '1', cwd, stdout=stdout, stderr=stderr, clobber = clobber)

        # plots lc
        plot_lc()

        # generates times to search for events
        evt_tstart, evt_tstop = get_time_windows()

        # print here? 'Running Search for swift_id'
        for i in range(len([evt_start])): # progress bar

            # produces evt_dpi, and sky_image, looks for sources
            wrapper_evt_detect_heasoft(swift_id, evt_file, tstart, tstop,
                               incatalog='test_cat', outcatalog = 'out.cat',
                               cwd=results_dir, snrthresh='3.5', stdout=stdout,
                               stderr=stderr, clobber = clobber)


            # searches outcatalog from batcelldetect for sources found near target, if found writes to file / break
            # at the moment just prints if it found something 
            search_out_catalog(outir=results_dir, outcatalog='out.cat', RA_target=ra,
                               DEC_target=dec, radius_2=5**2)

        # if batcelldetect does not detect a source IE low SNR some function to establish flux limit?
        # make noise image (maybe make this as a default output?)

    return
#-----------------------------------------------------------------------------#
# non-heasoft functions related to running search
def read_in_catalog_file(cat_file):

    # expected format
   #'swift_id' - dtype = 'S30', 'CHIME_id', 'S30', 'ra', '>f4', 'dec', '>f4', 'ra_err', '>f4','dec_err', '>f4',
   # 'name', 'S30'
    with fits.open(cat_file) as hdul:
            #hdul.info()
            header_dict = hdul[0].header

            data = hdul[1].data

    return data

# need to remake / not finished
def plot_lc(): # plots lc
    return None # placeholder

def search_out_catalog(outir, outcatalog, RA_target, DEC_target, radius_2=5**2, evt_window=1, id='test'):
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

    return np.abs((RA_target-RA_source)**2 + (DEC_target-DEC_source)**2)

def get_time_windows():
    return [tstart],[tend]

def inspect_evt_file():
    return evt_file, tstart, tstop

#-----------------------------------------------------------------------------#
def main():
    import argparse
    """
    Start the function to perform RM-synthesis if called from the command line.
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
    make_plots = True
    show_plots = True

    run_search()
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    main()
