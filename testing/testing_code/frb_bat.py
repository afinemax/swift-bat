#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2022 Maxwell A. Fine                                          #
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
#-------------------------------------------------------------------------------
# arguments
results = '../results' # maybe I'll save stuff in the corresponding swift id dir? that seems reasonable
file_ext = '.pdf' # for the outputted plots
data_dir = '../data' # will loop through every dir for the bat data
sky_local = 'sky_frb.tsv'
guano_dump = 'GUANO dump inventory - GUANO triggers.tsv'
# args
clobber = True # if true remakes/overwrites files
time_res = 1 # seconds?
energy_bins = 14-175
source_id = '00010374156'
evt_file = 'sw00010374156bevshsl_uf.evt.gz'
#tstart=' + str(evt_start)639072334.000000 tstop=639072896.800600'


events = '/bat/event'
event_dir = data_dir + '/' + source_id + events
hk_dir = data_dir + '/' +source_id + '/bat/hk'
results_dir = results + '/' + source_id
#-----------------------------------------------------------------------------#
def get_swift_ids_data(data_dir):
    '''make list of swift IDS from the data dir'''
    output = subprocess.check_output(["ls",],universal_newlines=True, cwd= data_dir)
    source_ids = output.split() # list of strings the source IDs

    # pop off these dirs if they exist in the list
    source_ids.remove('GUANO')
    source_ids.remove('dump')
    source_ids.remove('inventory')
    source_ids.remove('-')
    source_ids.remove('GUANO')
    source_ids.remove('triggers.tsv')
    source_ids.remove('skycoords.dat')
    return source_ids
#-----------------------------------------------------------------------------#
def get_sky_local(data_dir, sky_local):
    sky_local_path =  data_dir + '/' + sky_local
    sky_local_data = np.genfromtxt(sky_local_path,delimiter='\t',skip_header=1)

    # format chime Event no, status, RA, RA Error,Dec ,Dec error
    # we want a dict of chime event with ra, dec and errors!

    chime_ids_raw = sky_local_data[:,0]
    chime_ids = chime_ids_raw.astype(int)

    # degrees
    ra = sky_local_data[:,2]
    ra_err = sky_local_data[:,3]
    dec = sky_local_data[:,4]
    dec_err = sky_local_data[:,5]

    # make dicts
    chime_ra = dict(zip(chime_ids,ra))
    chime_dec = dict(zip(chime_ids,dec))
    chime_ra_err = dict(zip(chime_ids,ra_err))
    chime_dec_err = dict(zip(chime_ids,dec_err))

    return chime_ra, chime_dec, chime_ra_err, chime_dec_err
#-----------------------------------------------------------------------------#
def get_gauno_info(data_dir,guano_dump):

    guano_path = data_dir + '/' + guano_dump


    gauno_dump_data = np.genfromtxt(guano_path,delimiter='\t', skip_header=1, dtype= None)
    swift_ids = np.zeros(len(gauno_dump_data))
    chime_trigger = []
    chime_ids = []
    tsar_bb = []
    index_del = []

    # we want the ones with tsar_good and bb == True right?
    for i in range(len(gauno_dump_data)):
        swift_ids[i] = gauno_dump_data[i][0]
        #chime_trigger.append(gauno_dump_data[i][1]) # this is a b'string' format
        #chime_ids.append(gauno_dump_data[i][2]) # this is a b'string' format
        tsar_bb.append(gauno_dump_data[i][-2]) # this is a b'string' format
        if tsar_bb[i] == b'TRUE':
            tsar_bb[i] = True
            chime_trigger.append(gauno_dump_data[i][1]) # this is a b'string' format
            chime_ids.append(gauno_dump_data[i][2])
        elif tsar_bb[i] == b'FALSE':
            tsar_bb[i] = False
        else:
            index_del.append(i)
        # I assume this stays in order but gauno_dump format should change
    tsar_arr = np.asarray(tsar_bb[0:min(index_del)])
    index_arr = np.where(tsar_arr == True)
    swift_ids = swift_ids[index_arr[0]]


    return swift_ids, tsar_arr, chime_ids, chime_trigger
#-----------------------------------------------------------------------------#
def str_swift_id(swift_ids):

    str_ids = []

    for i in range(len(swift_ids)):
        id_len = len(str(int(swift_ids[i])))
        str_ids.append('0' * (11-id_len) + str(int(swift_ids[i])))


    return str_ids

#-----------------------------------------------------------------------------#
def get_swift_time_from_chime(chime_trigger, swift_ids):
    chime_trigger_str =[]
    for i in range(len(chime_trigger)):
        chime_trigger_str.append(chime_trigger[i].decode('utf-8'))
    trigger = chime_trigger_str
    cc = Clock()
    cc.utctime = trigger
    cc.submit()
    cc.met
    res = dict(zip(swift_ids, cc.met))
    return res
#-----------------------------------------------------------------------------#
def swift_chime_dict(swift_ids, chime_ids):
    '''returns a dict mapping swift_ids to chime_ids'''

    # need to change chime_id into int from b_string
    chime_ids_int = np.zeros(len(chime_ids))
    for i in range(len(chime_ids_int)):
        chime_ids_int[i] = int(chime_ids[i].decode("utf-8"))

    return dict(zip(swift_ids,chime_ids_int))
#-----------------------------------------------------------------------------#
def swift_sky_local(data_dir, sky_local,swift_ids, chime_ids):
    '''returns dict mapping corresponding chime localizations to swift_ids
    - checks to be sure that data exists'''

    chime_ra, chime_dec, chime_ra_err, chime_dec_err = get_sky_local(data_dir, sky_local)
    swift_chime_map = swift_chime_dict(swift_ids, chime_ids)

    swift_ra = {}
    swift_ra_err = {}
    swift_dec = {}
    swift_dec_err = {}
    # need to check that chime_id exists in both before making new entry
    for i in range(len(swift_ids)):
        swift_id = swift_ids[i]
        try:
            chime_id = int(swift_chime_map[swift_id])
        except:
            continue
        else:
            try:
                swift_ra[swift_id] = chime_ra[chime_id]
                swift_ra_err[swift_id] = chime_dec[chime_id]
                swift_dec[swift_id] = chime_ra_err[chime_id]
                swift_dec_err[swift_id]= chime_dec_err[chime_id]
            except:
                continue


    return swift_ra, swift_ra_err, swift_dec, swift_dec_err
#-----------------------------------------------------------------------------#
def get_evt_files(event_dir):
    'makes a list of the evt files in a swift bat dir'
    bat_event_data = subprocess.check_output(["ls",],universal_newlines=True,
    cwd= event_dir)
    evt_files = bat_event_data.split()
    return evt_files
#-----------------------------------------------------------------------------#
def get_dir_trees(source_ids,data_dir,results):
    '''makes dicts of event, hk, and results dirs'''

    events = '/bat/event'


    event_dir = {}
    hk_dir = {}
    results_dir= {}
    for i in source_ids:
        source_id = i
        event_dir[i] = data_dir + '/' + source_id + events
        hk_dir[i] = data_dir + '/' +source_id + '/bat/hk'
        results_dir[i] = results + '/' + source_id

    return event_dir, hk_dir, results_dir
#-----------------------------------------------------------------------------#
def run_heasoft(source_id, evt_file, results, event_dir, hk_dir,
                results_dir, evt_start,evt_end,bck_start,bck_end,
                clobber=clobber ,energy_bins=14-175, time_res=1):
    '''Wrapper function that runs the main heasoft commands, and produces results

    Args:
    '''

    j = evt_file

    # make log file

    # make entry in results dir based on swift id
    with open('mkdir.txt','wt',) as f:
        subprocess.run(['mkdir ' + str(source_id)], cwd= results, shell=True, stdout = f, stderr =f )
    # think about what to do for log file for the mkdir


    with open(results_dir + '/' + j + '_log.txt','wt',) as f:
        #sys.stdout = open('output.txt','wt')

        #f = open(results_dir + "/log.txt", "wt",'a')
        # make entry in results dir based on swift id
        subprocess.run(['mkdir ' + str(source_id)], cwd= results, shell=True, stdout = f, stderr =f ) # change 0 to i when making wrapper!!
        #stderr=subprocess.STDOUT
        # make DPI
        subprocess.run(['batbinevt ' + str(j) + ' weighted=no '
        + 'outunits=counts outtype=dpi energy=-'
        +' clobber=' + str(clobber)
        +' outfile= frb.dpi'] , cwd = event_dir, shell=True, stdout = f, stderr =f)

        # move frb.dpi and (copy) det mask, event file. attitude file into results dir
        subprocess.run(['mv ' +str(event_dir) + '/frb.dpi '
                         + str(results_dir) + '/frb.dpi'  ],
                        shell=True, cwd = '../data', stdout = f, stderr =f)

        subprocess.run(['cp ' +str(hk_dir) + '/sw' + str(source_id) + 'bdecb.hk.gz '
                       + str(results_dir) + '/sw' + str(source_id) + 'bdecb.hk.gz'  ],
                       shell=True, cwd = '../data', stdout = f, stderr =f)

        subprocess.run(['cp data/' +str(source_id) + '/auxil/sw' + str(source_id) + 'sat.fits.gz '
                       + 'data/' + str(results_dir) + '/sw' + str(source_id) + 'sat.fits.gz' ],
                       shell=True, cwd = '..', stdout = f, stderr =f)

        subprocess.run(['cp data/' +str(source_id) + '/bat/event/' + str(j) + ' '
                       + 'data/' + str(results_dir) + '/' + str(j) ],
                       shell=True, cwd = '..', stdout = f, stderr =f)

        # make mask
        subprocess.run(['bathotpix detmask=sw' + source_id
                        + 'bdecb.hk.gz outfile = frb.mask infile = frb.dpi  clobber =' + str(clobber) ],
                       shell=True, cwd=results_dir, stdout = f, stderr =f)

        # make lc with heasoft
        # should I set energy bins here and time res?
        subprocess.run(['batbinevt detmask=frb.mask ' + str(j)
                    + ' timedel = ' + str(time_res)  + ' weighted=no outtype=lc energybins=- outfile =frb.heasoft.lc clobber ='
                    + str(clobber) ],
                       shell=True, cwd=results_dir, stdout = f, stderr =f)

       # make frb2_dpi
       # this one should contain our target source
       # need to add target_time and time size stuff
       # I am unsure if the timedel is correct here

       # the following steps require time windows!
        subprocess.run(['batbinevt '  + str(j)
                    + ' outfile =frb2.dpi energybins=' +str(energy_bins)
                    + ' clobber= ' + str(clobber) + ' outtype=dpi detmask=frb.mask timedel = 0'
                    + ' tstart=' + str(evt_start) + ' tstop=' + str(evt_end) ],
                   shell=True, cwd=results_dir, stdout = f, stderr =f)

        # make bkg.dpi
        # this one uses a off event time size to make a background dpi
        subprocess.run(['batbinevt '  + str(j)
                    + ' outfile =bkg.dpi energybins=' +str(energy_bins)
                    + ' clobber= ' + str(clobber) + ' outtype=dpi detmask=frb.mask timedel = 0'
                    + ' tstart=' + str(bck_start) + ' tstop=' + str(bck_end) ],
                   shell=True, cwd=results_dir, stdout = f, stderr =f)

        # make sky image
        subprocess.run(['batfftimage frb2.dpi attitude = sw' + str(source_id) + 'sat.fits.gz detmask=frb.mask'
                       + ' outfile = frb.sky bkgfile=bkg.dpi clobber =' + str(clobber)],
                       shell=True, cwd = results_dir, stdout = f, stderr =f)

        # run heasoft source finder
        subprocess.run(['batcelldetect frb.sky outfile=cat.fits'], shell=True, cwd = results_dir, stdout = f, stderr =f)

    return
#-----------------------------------------------------------------------------#
def run_search():
    """ runs search over all bat files in data dir that correspond to FRBS with
    localizations"""
    #swift_ids = get_swift_ids_data(data_dir)

    swift_ids, tsar_arr, chime_ids, chime_trigger = get_gauno_info(data_dir,guano_dump)
    swift_ids = str_swift_id(swift_ids) # makes swift ID into a str this is

    chime_ra, chime_dec, chime_ra_err, chime_dec_err = get_sky_local(data_dir, sky_local)

    swift_chime_trigger = get_swift_time_from_chime(chime_trigger, swift_ids)

    swift_ra, swift_ra_err, swift_dec, swift_dec_err = swift_sky_local(data_dir, sky_local,swift_ids, chime_ids)
    #swift_ids = get_swift_ids_data(data_dir)
    swift_ids = np.asarray(list(swift_ra.keys()))
    event_dir, hk_dir, results_dir = get_dir_trees(swift_ids,data_dir,results)

    print('Running HEASOFT analysis over BAT Files')
    for i in tqdm(range(len(swift_ids))): # tqdm is for progress bar
        #sleep(3)

        swift_id = swift_ids[i]

        evt_start = float(swift_chime_trigger[swift_id]) - 1.5
        evt_end = float(swift_chime_trigger[swift_id]) + 1.5

        bck_start = float(swift_chime_trigger[swift_id]) - 40
        bck_end = float(swift_chime_trigger[swift_id]) - 20

        evt_files = get_evt_files(event_dir[swift_id])
        for j in evt_files:
            run_heasoft(swift_id, j, results, event_dir[swift_id],
            hk_dir[swift_id], results_dir[swift_id],
            evt_start,evt_end,bck_start,bck_end , clobber=clobber,energy_bins=14-175, time_res=1,)

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

    run_search()
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    main()
