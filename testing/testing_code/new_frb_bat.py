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


# will be command line entries eventually
results = '../results'
file_ext = '.pdf' # for the outputted plots
data_dir = '../data'
datadir = data_dir
sky_local = 'sky_frb.tsv'
guano_dump = 'GUANO dump inventory - GUANO triggers.tsv'
# args
clobber = True # if true remakes/overwrites files
time_res = 0.1 # seconds
energy_bins = 14-175
source_id = '00010374156'
swift_id = source_id
evt_file = 'sw00010374156bevshsl_uf.evt.gz'
snr_cutoff = 5
#tstart=' + str(evt_start)639072334.000000 tstop=639072896.800600'


events = '/bat/event'
event_dir = data_dir + '/' + source_id + events
hk_dir = data_dir + '/' +source_id + '/bat/hk'
results_dir = results + '/' + source_id
outdir = results_dir
#
cat_file = 'test_cat'

evt_start = 618680512.304
evt_end = 618680513.328
bck_start = evt_start -40
bck_end = evt_start -20

tstart = evt_start
tend =evt_end

tstop = evt_end

ra = 15.8908
dec = -73.8580
j = evt_file
#-----------------------------------------------------------------------------#
# heasoft functions
# need to add a function to create a noise image
def create_heasoft_outdir(swift_id, cwd, stdout=None, stderr=None):
    '''Creates outdir for outputs to be placed in'''

    subprocess.run(['mkdir -p ' + str(swift_id)], cwd= cwd, shell=True,
                        stdout=stdout, stderr=stderr)
    return

def move_ess_heasoft_files(swift_id, evt_file, outdir, datadir,
                            stdout=None, stderr = None):
    '''Moves necessary SWIFT/BAT files into outdir for convience in generating
    heasoft outputs'IMPORTANT relative paths need to be updated / set to
    variables for final script!'''

    hk_dir = datadir + '/' + swift_id + '/bat/hk'

# IMPORTANT relative paths need to be updated / set to variables for final script
    subprocess.run(['cp ' +str(hk_dir) + '/sw' + swift_id + 'bdecb.hk.gz '
                    + str(outdir) + '/sw' + swift_id + 'bdecb.hk.gz'  ],
                    shell=True, cwd = '../data', stdout = stdout,
                    stderr =stderr)

    subprocess.run(['cp data/' +str(swift_id) + '/auxil/sw' + swift_id
                    + 'sat.fits.gz '+ 'data/' + str(outdir) + '/sw'
                    + swift_id + 'sat.fits.gz' ],
                    shell=True, cwd = '..', stdout = stdout, stderr =stdout)

    subprocess.run(['cp data/' +str(swift_id) + '/bat/event/' + str(evt_file) + ' '
                    + 'data/' + str(outdir) + '/' + str(evt_file) ],
                    shell=True, cwd = '..', stdout = stdout, stderr =stderr)
    return


def heasoft_mask(swift_id, evt_file, outdir, stdout=None,
                    stderr=None, clobber = clobber):
    '''Procduces SWIFT/bat detector mask by running 'bathotpix
        and 'batbinevt' '''

    subprocess.run(['batbinevt ' + str(evt_file) + ' weighted=no outunits=counts'
                    + ' outtype=dpi energy=- clobber=' + str(clobber)
                    + ' outfile= frb.dpi.mask'], cwd = outdir, shell=True,
                    stdout = stdout, stderr =stderr)

    subprocess.run(['bathotpix detmask=sw' +  swift_id + 'bdecb.hk.gz outfile ='
                    +' frb.mask infile= frb.dpi.mask clobber =' + str(clobber) ],
                    shell=True, cwd=outdir, stdout = stdout, stderr =stderr)
    return

def heasoft_lc(evt_file, energy_bins, timdel, cwd=outdir, stdout=None, stderr=None, clobber = clobber):
    ''' Creates lightcurve of specified SWIFT/BAT evt_file using its assocaited mask

    Inputs:
    evt_file = string, file name of event file
    energybins = string, specified energy binning in KeV
    timdel = specified timeresolution '''

    subprocess.run(['batbinevt detmask=frb.mask ' + str(evt_file)
                    + ' timedel = ' + str(timdel)  + ' weighted = no outtype ='
                    + ' lc energybins = ' + energy_bins
                    + ' outfile = frb.lc clobber ='
                    + str(clobber) ],
                    shell=True, cwd=cwd, stdout = stdout, stderr =stderr)
    return

def heasoft_evt_dpi(evt_file, tstart, tstop, cwd,
                      stdout=None, stderr=None, clobber = clobber):
    '''Runs 'batbinevt' to create event dpi for the given evt_file
    from tstart to tstop'''

    # this is the same as heasoft_bgck_dpi besides output name
    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = evt.dpi energybins=- timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=dpi detmask=frb.mask weighted=no timedel=0.0'
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True,  )
    return

def heasoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=None, stderr=None, clobber = clobber):
    '''Runs 'batbinevt' to create background dpi for the given evt_file
        from tstart to tstop'''

    # this is the same as heasoft_evt_dpi besides output name
    subprocess.run(['batbinevt '  + str(evt_file)
                    + ' outfile = bkg.dpi energybins=- timebinalg=u'
                    + ' clobber= ' + str(clobber)
                    + ' outtype=dpi detmask=frb.mask weighted=no timedel=0.0'
                    + ' tstart=' + str(tstart) + ' tstop=' + str(tstop)],
                       stdout=stdout , stderr=stderr, cwd=cwd , shell=True,  )

    return

def heasoft_sky_image(swift_id, cwd, stdout=None, stderr=None, clobber = clobber):
    ''' Creates sky image using 'batfftimage' and the evt.dpi and bgck.dpi files '''

    subprocess.run(['batfftimage evt.dpi attitude = sw' + swift_id + 'sat.fits.gz detmask=frb.mask'
                       + ' outfile = frb.sky bkgfile=bkg.dpi clobber =' + str(clobber)],
                       shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return

def heasoft_batcelldetect(sky_image, cwd, incatalog, outcatalog='cat.fits', snrthresh='3.5',
                           stdout=None, stderr=None, clobber = clobber):
    '''Runs 'batcelldetect' over sky image'''

    subprocess.run(['batcelldetect ' + str(sky_image) + ' snrthresh = ' + str(snrthresh)
                   + ' incatalog = ' + str(incatalog) +' outfile=' + str(outcatalog)],
                  shell=True, cwd = cwd, stdout = stdout, stderr =stderr)
    return


def wrapper_bgck_heasoft(swift_id, outdir, evt_file, tstart,
                         tstop, datadir, stdout=None, stderr = None):
    '''Wrapper function: moves files, produces mask , and produces background dpi'''

    move_ess_heasoft_files(swift_id, evt_file, outdir, datadir, stdout=None, stderr = None)

    cwd=outdir

    heasoft_mask(swift_id, evt_file, outdir, stdout=stdout, stderr=None, clobber = clobber)

    heasoft_bgck_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    return


def wrapper_evt_detect_heasoft(swift_id, evt_file, tstart, tstop,
                               incatalog='test_cat', outcatalog = 'out.cat',
                               cwd=outdir, snrthresh='3.5', stdout=None,
                               stderr=None, clobber = clobber):

    '''Wrapper function: makes evt_dpi, sky image , and runs batcelldetect'''

    heasoft_evt_dpi(evt_file, tstart, tstop, cwd,
                      stdout=stdout, stderr=stderr, clobber = clobber)

    heasoft_sky_image(swift_id, cwd, stdout=stdout, stderr=stderr, clobber = clobber)

    heasoft_batcelldetect('frb.sky', cwd, incatalog, outcatalog, snrthresh=snrthresh,
                           stdout=stdout, stderr=stdout, clobber = clobber)
    return
#-----------------------------------------------------------------------------#
# non-heasoft functions related to running search

def run_search(incatalog, outdir, datadir, energy_bins, timedel,
               clobber=clobber, snrthresh='3.5', log=print,swift_evt_file_ending='bevshpo_uf.evt.gz'):
    '''Main program'''

    swift_id, ra, dec, ra_err, dec_err, swift_trig_time, chime_id = read_in_catalog_file(incatalog)

   # log('running search')
    for i in tqdm(range(len(swift_id)),desc='search over '+ str(len(swift_id)) + ' targets'):
        # need to add
        # download swift/bat data for target if doesnt already exist

        resultsdir = outdir + '/' + swift_id[i]
        bck_tstart = str(float(swift_trig_time[i]) - 70)
        bck_tstop = str(float(swift_trig_time[i]) - 20)
        print('On target')
        run_search_on_target(swift_id[i], ra[i], dec[i], energy_bins, timedel,
                             incatalog, resultsdir, datadir,
                             swift_trig_time[i], ra_err[i], dec_err[i],
                             bck_tstart, bck_tstop, snrthresh,log, outdir)

    return

def run_search_on_target(swift_id, ra, dec, energy_bins, timdel, incatalog, outdir, datadir,
                         trigger_time, ra_error, dec_error, bck_tstart,
                         bck_tstop, snrthresh, log=print, swift_evt_file_ending='bevshpo_uf.evt.gz'):
    '''Main function runs search on single SWIFT/BAT target '''

    results = outdir[0:-(len(swift_id) + 1)]
    create_heasoft_outdir(swift_id, outdir=results, stdout=None, stderr=None)

    with open(outdir + '/' + 'log.txt','wt',) as f:
        stdout = f
        stderr = f


        print('starting ')
        evt_file = get_evt_file(swift_id, swift_evt_file_ending='bevshpo_uf.evt.gz')


        # check for existance of data around chime trigger
        #evt_file, bck_tstart, bck_tstop = inspect_evt_file()

        # move files into outdir, produce mask, bck_dpi
        wrapper_bgck_heasoft(swift_id, outdir, evt_file, bck_tstart,
                         bck_tstop, datadir, stdout=stdout, stderr = stderr)

        # might want to move catalog file into outdir for convience?

        # produce lc
        heasoft_lc(evt_file, energy_bins, timdel,
                    cwd=outdir, stdout=stdout, stderr=stderr, clobber = clobber)

        # plots lc
        time_window='10'
        lc_file = 'frb.lc'
        light_curve_analysis(swift_id, trigger_time,
        time_window, lc_file, energy_bins, outdir)

        # generates times to search for events
          #evt_tstart, evt_tstop = get_time_windows(trigger_time)
        evt_start = trigger_time
        evt_tstop = trigger_time + 1


        # progress bar
        for i in tqdm(range(len([evt_start])),leave=False,
                      desc='searching ' + str(len([evt_start])) +' time windows'):

            # produces evt_dpi, and sky_image, looks for sources
            wrapper_evt_detect_heasoft(swift_id, evt_file, evt_start, evt_tstop,
                               incatalog=incatalog, outcatalog = 'out.cat',
                               cwd=outdir, snrthresh=snrthresh, stdout=stdout,
                               stderr=stderr, clobber = clobber)


            # searches outcatalog from batcelldetect for sources found near target, if found writes to file / break
            # at the moment just prints if it found something
            search_out_catalog(outdir=outdir, outcatalog='out.cat', RA_target=ra,
                               DEC_target=dec, radius_2=5**2)

        # if batcelldetect does not detect a source IE low SNR some function to establish flux limit?
        # make noise image (maybe make this as a default output?)

    return


def light_curve_analysis(swift_id, trig_time, time_window, lc_file, energy_bans, outdir):
    '''Searches for SNR peaks in the light curve in the search window by computing a
    running average with logscale time-window sizes
    returns SNR_peak_value, SNR_peak_time (compared to trig time), LC_analysis'''

    # this works but is kinda slow
    # is there a faster way?

    with fits.open(outdir + '/' + lc_file) as hdul:
        #hdul.info()
        header_dict = hdul[0].header
        time = hdul[1].data["Time"]
        rate = hdul[1].data['rate']
        error = hdul[1].data['error']
        totcounts = hdul[1].data['totcounts']
        # fracexp = hdul[1].data['FRACEXP']
        # columns = hdul[1].data.columns
    bans = energy_bans.split(',')
    bans

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
    for k in range(len(search_time_res)):
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

        #    start_i = (np.abs(new_time - (trig_time - time_window ))).argmin()
            #stop_i = (np.abs(new_time - (trig_time + time_window))).argmin()

            #peak_snr_arr_new[k, i] = np.max(snr[start_i:stop_i])

            peak_snr_arr[k,i] = np.max(snr)

    peak_snr = np.max(peak_snr_arr[:,0])
    k = np.where(peak_snr == peak_snr_arr[:,0])[0][0]
    i = np.where(peak_snr == peak_snr_arr)[0][0]
    idx = (np.abs(snr_dict[search_time_res[k]][bans[i]] - peak_snr)).argmin()
    time_peak_snr = time_dict[search_time_res[k]][bans[i]][idx]

    results_dict['time_peak_snr'] = time_peak_snr
    results_dict['peak_snr'] = peak_snr
    results_dict['peak_time_scale'] = search_time_res[k]
    results_dict['time_offset'] = trig_time - time_peak_snr

    plot_lc(swift_id, trig_time, rate_dict, snr_dict,
    search_time_res, time_dict, peak_snr,
    peak_snr_arr, bans, outdir)


    return results_dict, peak_snr_arr_new, peak_snr_arr

def plot_lc(swift_id, trig_time, rate_dict, snr_dict,
         search_time_res, time_dict, peak_snr,
         peak_snr_arr, bans, outdir):

    'Plots LC with best found time resolution (from searching on a log basis), and plots rate, and SNR'
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(12, 6), dpi=80)
    k = np.where(peak_snr == peak_snr_arr[:,0])[0][0]


    # rate plot
    for i in range(len(bans)):
        time = time_dict[search_time_res[k]][bans[i]]
        rate_i = rate_dict[search_time_res[k]][bans[i]]
        plt.plot(time-time[0], rate_i, label=bans[i] +' KeV')
    plt.ylabel('Rate ($s^{-1}$)')
    plt.xlabel('Time ($s$)')
    trig_time_plt = find_nearest(time, trig_time)
    plt.axvline(x=trig_time_plt-time[0], color ='black', linestyle=':', label='Trigger Time')
    legend = plt.legend(bbox_to_anchor =(1, 1))#,loc='upper right')
    legend.set_title(r'Time Resolution '+ str(search_time_res[k]) + '$s$')
    plt.title('Light Curve for ' + str(swift_id) )
    plt.savefig(outdir + '/lc_plot' + str(swift_id)+ '.pdf')

    # SNR plot
    plt.rcParams.update({'font.size': 20})
    plt.figure(figsize=(12, 6), dpi=80)
    for i in range(len(bans)):
        time = time_dict[search_time_res[k]][bans[i]]
        snr_i = snr_dict[search_time_res[k]][bans[i]]
        plt.plot(time-time[0], snr_i, label=bans[i] +' KeV')
    plt.ylabel('SNR')
    plt.xlabel('Time ($s$)')
    trig_time_plt = find_nearest(time, trig_time)
    plt.axvline(x=trig_time_plt-time[0], color ='black', linestyle=':', label='Trigger Time')

    # this assumes peak is in the last index for i
    snr_peak_time = find_element(snr_i, time, peak_snr)
    plt.axvline(x=snr_peak_time-time[0], color ='purple', linestyle='-.', label='Peak SNR Time (in window)')
    legend = plt.legend(bbox_to_anchor =(1, 1))#,loc='upper right')
    legend.set_title(r'Time Resolution '+ str(search_time_res[k]) + '$s$')
    plt.title('SNR for ' + str(swift_id) )
    plt.savefig(outdir +'/SNR_plot' + str(swift_id)+ '.pdf')

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def find_nearest(array, value):
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
    '''Creates outdir for outputs to be placed in'''

    subprocess.run(['mkdir -p ' + str(swift_id)], cwd=outdir, shell=True, stdout=stdout, stderr=stderr)
    return

def read_in_catalog_file(cat_file):
    'WIP for format, as at the moment batcelldetect is not working as intended'

    # I think part of the problem for batcelldetect is my dtype could be wrong?
    # swift_id ra, dec, ra_err, dec_err, swift_trig_time, chime_id
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

def get_evt_file(swift_id, swift_evt_file_ending='bevshpo_uf.evt.gz'):
    '''Returns string of filename for SWIFT/bat event file based on its swift id'''
    evt_file = 'sw' + swift_id + swift_evt_file_ending
    return evt_file

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
    outdir='../results'
    datadir = '../data'
    #timedel = 1
    timedel = '1e-4'#'1e-4'
    incatalog = 'test_cat'
    energy_bins = '5-25,25-50,50-100,100-350'
    swift_evt_file_ending='bevshpo_uf.evt.gz'
    clobber = True
    snrthresh = '6'
    log=print
    time_window='10'
    lc_file = 'frb.lc'

    run_search(incatalog, outdir, datadir, energy_bins, timedel,
               clobber=clobber, snrthresh=snrthresh, log=print,
               swift_evt_file_ending=swift_evt_file_ending)
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    main()
