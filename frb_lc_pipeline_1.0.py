 #!/usr/bin/env/python
#=============================================================================#
# NAME:                                                                       #
# INPUTS:                                                                     #
# OUTPUTS:                                                                    #
# PURPOSE:                                                                    #
#                                                                             #
# LAST MODIFIED: January 2023 by Maxwell A. Fine                              #
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

'''
expects a .csv file with the followiing

swift_id, trig_time, chime_id


'''

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


def make_outdir(outdir):
    '''makes outdir dir for results dirs to go into'''
    os.makedirs(outdir, exist_ok = True)


def read_in_catalog(incatalog):
    data = np.genfromtxt('frb_test_cat_jan_2023.csv', dtype=str, delimiter=',')
    swift_ids = data[:,0]
    trig_time = data[:,1]
    chime_ids = data[:,2]

    # convert trig time into swift time
    # there must be a faster way of doing this
    #swift_time_arr = np.zeros(len(trig_time))
    #for i in range(len(trig_time)):


    return swift_ids, trig_time, chime_ids


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

    #sw00096133003bevshpo_uf.evt.gz

    # datadir = ../data
    # outdir = ../results

    #/bigger dir/ # this file

    #bigger dir = ..




    subprocess.run(['cp ' + cwd + '/' +str(swift_id) + '/bat/event/' + str(evt_file)
                + ' ' + str(outdir) + '/' + str(evt_file) ],
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


# read csv in

# weight function

# lc analysis, open lc, lc plotting

# function to export results to a json
# snr thresh = 10
# terminal args
# snr thresth, catalog, energybans, time res, and number of samples
#

# def main
# gets csv data
# downloads swift data
# moves files
# opens evt_file and checks if its valid
# lc analysis
# plotting
# export results
# report summary


def get_lightcurve(outdir, evt_file, time_resolution, energy_bans, lc_file):
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


    #energy_bans = '15-25,25-50,50-100,100-150'

    # rate = 2d array of [ban1, ban2, ban3, ban4] counts

    #total counts

    bins=int(np.ptp(time)/time_resolution)
    total_counts, bins = np.histogram(time, bins=bins,)
    new_time = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]

    # now we need to do each energy ban
    # we need to split energy bans to get the ranges
    # split by , then by -
    # do a np.where along the energy arrary to sort them into the bans
    #
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
       # new_time = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]

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

    time_resolution = 1e-3
    print('swift id', swift_id)
    lc_file = 'sw' + swift_id + evt_file

    outdir = outdir +'/' + swift_id

    # open lc file, read out data

    time, totcounts, rate, old_time, energy = get_lightcurve(outdir,
                                                                              lc_file,
                                                                              time_resolution,
                                                                              energy_bans, lc_file)

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
    time_res = np.logspace(-2, np.log10(2*time_window), 10) #, 0.01, 0.1, 1
    for i in tqdm(range(len(time_res)),leave=False,
                  desc='searching ' + str(len(time_res)) +' time windows'):


        time_size_per = int(np.abs((time[-1] - time[0]))) / len(time) # time size per index
        window_len = int(time_res[i] / time_size_per)
        window = np.ones(window_len,dtype='float32')/ np.sqrt(window_len)


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

    return lc_analysis_dict, rate_snr_dict, time, energy


def expected_snr(n_samples):

    snr = special.erfinv((-1/n_samples) + 1) * np.sqrt(2)

    # this sometimes returns negative values
    # we only care about the + ones


    return np.abs(snr) # return postive value only


def lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, trig_time,
                energy_bans, time_window, outdir, chime_id, swift_id):
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

    fig, axs = plt.subplots(2, 2, figsize=(16,9),
        gridspec_kw={'width_ratios': [1, 1.75]})
    fs =15
    plt.yticks(fontsize=fs, color='k')
    plt.xticks(fontsize=fs, color='k')


    # top left
    counts_hist, bins_hist, etc = axs[0,0].hist(energy, bins=int(np.sqrt(len(energy))), color='green')
    axs[0, 0].set_title('Energy Distrubution', fontsize=fs, color='k')
    axs[0,0].set_xlabel('KeV', fontsize=fs, color='k')
    axs[0,0].set_ylabel('N', fontsize =fs, color='k')
    axs[0,0].set_xlim(0, 350) # look into xlim
#    axs[0,0].set_ylim(0, np.max(counts_hist)*1.05)

    # bottom left
    #axs[1, 0].plot(x, -y, 'tab:green')
    axs[1, 0].set_title('Axis [1, 0]')
    # lets loop through and grap all SNR vs time scale
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

    # now loop through energy bans
    my_keys = list(lc_analysis_dict['bans_dict'].keys())
    for i in range(len(my_keys)):
        energy_ban = my_keys[i]
        time_res_keys = list(lc_analysis_dict['bans_dict'][energy_ban].keys())

        time_res = []
        snr_val = []
        for j in range(len(time_res_keys)):
            time_res.append(float(time_res_keys[j]))
            snr_val.append(lc_analysis_dict['bans_dict'][energy_ban][time_res_keys[j]]) # changed to j from i
            # if statment for getting max SNR, and peak time scale

        axs[1, 0].scatter(time_res, snr_val, label=energy_ban + ' KeV')

        if np.max(snr_val) > peak_snr:
            peak_snr = np.max(snr_val)
            peak_time_scale = time_res(np.argmax(snr_val))


    time_resolution = peak_time_scale



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


    # right side


    legend_time_res = 'Time Resolution: '+ str(time_resolution) + ' s'

    snr_arr_tot = lc_analysis_dict['totcounts_dict']['snr_arr']

    #fig.suptitle()
    axs[0, 1].set_title('SNR of Lightcurve', fontsize=15, color='k')
    axs[1, 1].set_title('SNR of Search Window', fontsize=15, color='k')

    axs[0, 1].set_xlabel('Time (s)', fontsize=fs, color='k')
    axs[0, 1].set_ylabel('SNR', fontsize=fs, color='k')
    axs[0, 1].step(time - time[0],snr_arr_tot, color='purple')
    axs[0,1].axvline(x=trig_time - time[0], linestyle='--', color='k',
                    label="Trigger Time")
    axs[0,1].legend(fontsize=fs-2,title=legend_time_res)

    axs[1, 1].set_xlabel('Time (s)', fontsize=fs, color='k')
    axs[1, 1].set_ylabel('SNR', fontsize=fs, color='k')
    axs[1, 1].step((time - time[0])[start_i:stop_i],
    (snr_arr_tot)[start_i:stop_i], color='purple')
    axs[1,1].axvline(x=trig_time - time[0], linestyle='--', color='k',
                    label="Trigger Time")
    axs[1,1].legend(fontsize=fs-2,title=legend_time_res)


    fig.suptitle("SWIFT ID " + swift_id + '     Peak SNR: '+ str(peak_snr)[0:4]
    + '\nCHIME ID ' + chime_id +'    Peak Timescale ' + str(time_resolution)[0:4] +'s', size=fs)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    plt.savefig(outdir + 'diagnostic.pdf', bbox_inches='tight')
    plt.close()

    # returns peak snr, and peak time scale
    # want to also return an array of the peak snr, and time scales

    return peak_snr, time_resolutions


def outcatalog(outcatalog_file, swift_id, chime_id, timescale, peak_snr):

    with open(str(outcatalog_file), 'a') as file:
        # this doesn't delete and remake it every time I run it
        # hmm


    # Make a list of strings
        arr = [swift_id, chime_id, str(timescale), str(peak_snr)]

    # Using the writelines() function, write all the sentences of the array to the file.
        file.writelines(arr)
    return


def main():
    # read incatalog
    # need to add comand line interface
    energy_bans = '15-25,25-50,50-100,100-150'
    evt_file = 'bevshpo_uf.evt.gz'
    datadir = '../data_lc'
    outdir = '../lc_pipeline_res'
    incatalog = 'frb_test_cat_jan_2023.csv'
    time_window = 3
    swift_ids, trig_time, chime_ids = read_in_catalog(incatalog)
    outcatalog_file = 'test_outcatalog.txt'

    # delete outcatalog

    swift_ids = swift_ids

    # for each swift_id we make an outidr, and download the data
    for i in range(len(swift_ids)):

         # make outdir
         make_outdir(outdir + '/' + swift_ids[i])

    # download data
    all_downloaded  = True
    while all_downloaded is True:
        try:
            get_swift_bat_data(swift_ids, datadir, overwrite=False)
            all_downloaded = False
        except:
            all_downlaoded = True



    for i in range(len(swift_ids)):

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

            lc_analysis_dict, rate_snr_dict, time,energy = analysis_result

            my_outdir = outdir + '/' + str(swift_ids[i] + '/')

            peak_snr, time_resolution = lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, my_trig_time,
                        energy_bans, time_window, my_outdir, chime_ids[i], swift_ids[i])

            # write outcatalog

            outcatalog(outcatalog_file, swift_ids[i], chime_ids[i], time_resolution, peak_snr)

            if peak_snr > 12:
                print('SNR > 12!')



        except:
           print('failed lc analysis')


if __name__ == '__main__':
    main()
