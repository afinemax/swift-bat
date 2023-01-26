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
    swift_time_arr = np.zeros(len(trig_time))
    for i in range(len(trig_time[0:5])):
        cc = Clock()
        #datetime(2016,12,31,23,59,59)
        cc = Clock(swifttime=str(trig_time[i]))
        swift_time_arr[i] = cc.met
        #print(swift_time)

    return swift_ids, swift_time_arr, chime_ids


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
    print('trig_time', trig_time)
    lc_file = 'sw' + swift_id + evt_file

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
    time_res = np.logspace(-2, 2, 15) #, 0.01, 0.1, 1
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

def lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, trig_time,
                energy_bans, time_window, outdir):
    '''produces plots of the SNR for the best time scale, and plots SNR vs time scale '''

    print('start ')
    # determine best SNR
    # I need to filter out the low SNR energy channel prior to this
    # will be done in LC analysiss or based on the energy bands
    # plot

    # indexs to cut SNR curves to search window size
    #start_i = np.where(np.isclose(time - trig_time + time_window, 0, rtol=1e-05) == True)[0][0]
    #stop_i = np.where(np.isclose(time - trig_time - time_window, 0, rtol=1e-05) == True)[0][0]

    start_i = np.where(np.abs(time - trig_time + time_window) == np.min(np.abs(time-trig_time + time_window)))[0][0]
    stop_i = np.where(np.abs(time - trig_time - time_window) == np.min(np.abs(time-trig_time - time_window)))[0][0]



    something = plt.hist(energy, bins=int(np.sqrt(len(energy))))
    fs =15
    plt.title('Energy Distrubution of Counts', fontsize=fs, color='white')
    plt.xlabel('KeV', fontsize=fs, color='white')
    plt.ylabel('N', fontsize =fs, color='white')
    plt.xlim(0, 550)

    plt.savefig(outdir + 'Energy_distro.pdf')
    plt.close()



    # plot snr vs time scale
    # how to best plot this?
    # 1 plot with different color scatters? for different enenrgy bans?
    fs = 12
    #fig = plt.figure(figsize=(8, 6))

    # lets loop through and grap all SNR vs time scale

    time_res = []
    snr_val = []
    my_keys = list(lc_analysis_dict['totcounts_dict'].keys())
    for i in range(len(my_keys) -1):
        i += 1
        time_res.append(float(my_keys[i]))
        snr_val.append(lc_analysis_dict['totcounts_dict'][my_keys[i]])



    plt.scatter(time_res, snr_val, label='Total Counts')
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

        plt.scatter(time_res, snr_val, label=energy_ban + ' KeV')

    plt.legend(fontsize=fs, labelcolor='k',facecolor='white', framealpha=1 )
    #plt.legend(facecolor='white', framealpha=1)
    plt.xlabel('Time Scale (s)', fontsize=fs, color='k')
    plt.ylabel('Peak SNR', fontsize =fs, color='k')
    plt.yticks(fontsize=fs, color='k')
    plt.xticks(fontsize=fs, color='k')
    plt.title('Time Scale vs Peak SNR', fontsize= 15)
    plt.xscale('log')
    plt.savefig(outdir + 'peak_SNR_vs_time_scale.pdf',bbox_inches='tight')
#    plt.show()
    plt.close()


    # plot snr tot lc vs time cutout region for search window
    snr_arr_tot = lc_analysis_dict['totcounts_dict']['snr_arr']

    fs = 12
    fig, ax = plt.subplots(2,figsize=(12, 12))
    plt.yticks(fontsize=fs, color='k')
    plt.xticks(fontsize=fs, color='k')



    #fig.suptitle()
    ax[0].set_title('SNR of Total Lightcurve', fontsize=15, color='k')
    ax[1].set_title('SNR of Search Window', fontsize=15, color='k')

    ax[0].set_xlabel('Time (s)', fontsize=fs, color='k')
    ax[0].set_ylabel('SNR', fontsize=fs, color='k')
    ax[0].step(time - time[0],snr_arr_tot)

    ax[1].set_xlabel('Time (s)', fontsize=fs, color='k')
    ax[1].set_ylabel('SNR', fontsize=fs, color='k')
    ax[1].step((time - time[0])[start_i:stop_i],(snr_arr_tot)[start_i:stop_i])

    plt.savefig(outdir +'SNR_tot_lc.pdf',bbox_inches='tight')
#    plt.show()
    plt.close()

    # plot the different energy bans vs time cutout region for search window
    rate_snr_keys = list(rate_snr_dict.keys())
    for i in range(len(rate_snr_keys)):
        my_key = rate_snr_keys[i]

        snr_arr = rate_snr_dict[my_key]

        fs = 12
        fig, ax = plt.subplots(2,figsize=(12, 12))
        plt.yticks(fontsize=fs)
        plt.xticks(fontsize=fs)


        #fig.suptitle()
        ax[0].set_title('SNR of ' + my_key + 'KeV Band', fontsize=15)
        ax[1].set_title('SNR of Search Window', fontsize=15)

        ax[0].set_xlabel('Time (s)', fontsize=fs)
        ax[0].set_ylabel('SNR', fontsize=fs)
        ax[0].plot(time - time[0],snr_arr)

        ax[1].set_xlabel('Time (s)', fontsize=fs)
        ax[1].set_ylabel('SNR', fontsize=fs)
        ax[1].plot((time - time[0])[start_i:stop_i],(snr_arr)[start_i:stop_i])

        plt.savefig(outdir +'SNR_' + my_key + 'lc.pdf',bbox_inches='tight')
    #    plt.show()
        plt.close()

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

    swift_ids = swift_ids[0:30]

    # for each swift_id we make an outidr, and download the data
    for i in range(len(swift_ids)):

         # make outdir
         make_outdir(outdir + '/' + swift_ids[i])

    # download data
    get_swift_bat_data(swift_ids, datadir, overwrite=False)


    print(len(swift_ids))
    outcatalog_dict = {}
    for i in range(len(swift_ids)):
        print(i, swift_ids[i])

        try:
            move_ess_heasoft_files(swift_ids[i], evt_file, outdir, datadir,
                                    stdout=None, stderr = None)

            analysis_result = lc_analysis(swift_ids[i],
                    evt_file, trig_time[i], time_window,  energy_bans, outdir)

            lc_analysis_dict, rate_snr_dict, time,energy = analysis_result

            my_outdir = outdir + '/' + str(swift_ids[i] + '/')

            lc_plotting(lc_analysis_dict, rate_snr_dict, time, energy, trig_time[i],
                        energy_bans, time_window, my_outdir )
            #outcatalog_dict{str(swift_id[i])} = 



        except:
            print('failed lc analysis')
            #


if __name__ == '__main__':
    main()
