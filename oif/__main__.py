#!/opt/local/bin/python3.4
###########################################
# Asteroid Survey Simulation Software
# Developed at the Jet Propulsion Laboratory
# Authors: Shantanu Naidu
#          Davide Farnocchia
#          Steve Chesley
# Date:    Apr 06, 2018
##########################################


import os
import numpy as np
import spiceypy as sp
import time
import configparser
import glob
import sys
import argparse
import warnings
import pyoorb as oo
import shutil
import pandas as pd

from . import shared
from . import telescope as ts
from . import sso as ss

DATAPATH = os.path.join(os.path.dirname(__file__), 'data')
USER_DATAPATH = None
def resolve_path(fn):
    # if absolute, return as-is
    if os.path.isabs(fn):
        return fn

    # if relative, check if there's one in cwd
    if os.path.exists(fn):
        return os.path.abspath(fn)

    # next, check if there's one in the user-overriden path (if any)
    if USER_DATAPATH is not None:
        f = os.path.abspath(os.path.join(USER_DATAPATH, fn))
        if os.path.exists(f):
            return f

    # else, return path pointing to DATAPATH
    return os.path.abspath(os.path.join(DATAPATH, fn))

def resolve_oorb_ephem_path(planetary_ephem):
    if os.path.isabs(planetary_ephem):
        return planetary_ephem

    # 1. first resolve against our data path
    depath = resolve_path(planetary_ephem)
    if os.path.isfile(depath):
        return depath

    # 2. then respect in OORB_DATA envvar, if given
    if "OORB_DATA" in os.environ:
        return os.path.abspath(os.path.join(os.environ["OORB_DATA"], planetary_ephem))

    # 3. finally, look relative to the openorb binary (in .../share/openorb/foo.dat)
    import shutil
    oorb_exe = shutil.which("oorb")
    if oorb_exe is not None:
        return os.path.join(os.path.dirname(os.path.dirname(oorb_exe)), 'share', 'openorb', planetary_ephem)


def get_or_exit(config, section, key, message):
    try:
        return config[section][key]
    except KeyError:
        sys.exit(message)

def main():
    #Parsing command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input  configuration filename.", type=str)
    parser.add_argument("-f", help="Force the program to replace asteroid SPKs", action="store_true")
    args = parser.parse_args()
    
    inputfile=args.input

    configdict={}

    # Set up reasonable defaults
    config = configparser.SafeConfigParser()
    config.read_dict({
        'CONF': {
            'n_Proc':           '1',		# Number of parallel processes
            'Planetary ephem':  'de430.dat',	# Planetary ephemeris file for OpenOrb
            'SPICE metakernel': 'meta_kernel', # NAIF SPICE kernels (Earth orientation, leap sec, etc.) (input)
            'Cache dir':        '_cache',       # Directory where to place generated temporary files
        },
        'ASTEROID' : {
            'Make SPKs':         'T',
            'nDays':             '4000',
            'SPK step':          '20',
            'nbody':             'F',
            'Object1':           '1',
            'nObjects':          '-1',
            'Make SPKs':         'T',
            'Asteroid SPKs':     'ast',		# Base file name and path for storing and retrieving asteroid spks (input/output)
            'Asteroid SPK path': 'spks',
        },
        'SURVEY': {
            'Survey T0':         '0',
            'Field1':            '1',
            'nFields':           '1000',
            'Space':             'F',
            'MPCobscode file':   'obslist.dat',
        },
        'CAMERA': {
            'Threshold': '5',
            'SPICE IK':  'camera.ti',
        }
    })
    config.read(inputfile)

    # Update the built-in data location, if overridden by user
    try:
        global USER_DATAPATH
        USER_DATAPATH = config["CONF"]['Data path']
    except KeyError:
        pass

    #CONF section
    spice_mk        = config["CONF"]['SPICE metakernel']
    planetary_ephem = config["CONF"]['Planetary ephem']
    cachedir        = config["CONF"]['Cache dir']
    nProc           = int(config['CONF']['n_Proc'])

    #CAMERA section
    cameradef_file    = get_or_exit(config, 'CAMERA', 'Camera', 'Camera FOV definition file not provided')
    spiceik           = config['CAMERA']['SPICE IK']
    threshold         = config.getfloat('CAMERA', 'Threshold')
    if (threshold > 90.0):
        warnings.warn('Threshold was > 90 degrees. Setting it to 90',Warning)
        threshold=90

    #ASTEROID section
    population_model = get_or_exit(config, 'ASTEROID', 'Population model', 'Population model file not provided')
    object1          = int(config['ASTEROID']['Object1'])
    if (object1==0):
        sys.exit('Object1 should >=0')
    nObjects         = config.getint('ASTEROID','nObjects')
    asteroidspks     = get_or_exit(config, 'ASTEROID', 'Asteroid SPKs', 'Asteroid SPK file basename not provided')
    asteroidspkpath  = get_or_exit(config, 'ASTEROID', 'Asteroid SPK path', 'Directory containing asteroid SPKs not provided')
    makespks         = config.get('ASTEROID', 'Make SPKs')
    spkstart         = float(get_or_exit(config, 'ASTEROID', 'SPK T0', 'Start time for creating SPK(s) not known'))
    spkndays         = float(get_or_exit(config, 'ASTEROID', 'nDays', 'Number of days for creating SPK not known'))
    spkstep          = float(get_or_exit(config, 'ASTEROID', 'SPK step', 'SPK step unspecified'))
    if (spkndays/spkstep < 11 and makespks == 'T'):
        spkstep=spkndays/11
        warnings.warn('Not enough steps to create SPKs. Reducing SPK step to %f' %(spkstep),Warning)
    nbody            = config.get('ASTEROID','nbody')
    if config.has_option('ASTEROID', 'stopTimesFile'):
        stopTimesFile = config.get('ASTEROID', 'stopTimesFile')
    else:
        # Handle the missing option (log, set a default, etc.)
        stopTimesFile = False

    #SURVEY section
    surveydb         = get_or_exit(config, 'SURVEY','Survey database', 'Survey database not provided')
    starttime        = config.getint('SURVEY','Survey T0')
    Field1           = config.getint('SURVEY','Field1')
    nFields          = config.getint('SURVEY','nFields')
    spaceflag        = config.get('SURVEY','Space')
    asteroidspks     = os.path.join(asteroidspkpath, asteroidspks)

    # resolve file locations relative to built-in data paths,
    # taking account of user overrides
    spice_mk              = resolve_path(spice_mk)
    cameradef_file        = resolve_path(cameradef_file)
    population_model      = resolve_path(population_model)
    surveydb              = resolve_path(surveydb)
    stopTimesFile         = resolve_path(stopTimesFile)
    if spaceflag=='T':
        obscode       = config.get('SURVEY','SCID')
        scspk = obscode + '.bsp'
        scspk = resolve_path(scspk)
    # Done reading configuration file

    #If it made it this far, print header
    with open(inputfile,'r') as f:
        for row in f:
            if(not row.startswith("#") and not row.startswith(";") and row.strip()):
                print(row,end='')

    # Changing directory to data path
    os.makedirs(cachedir, exist_ok=True)
    os.chdir(cachedir)

    # Initialize OpenOrb
    depath = resolve_oorb_ephem_path(planetary_ephem)
    try:
        oo.pyoorb.oorb_init(ephemeris_fname=depath)
    except Exception as e:
        print(e)
        sys.exit("Unable to load planetary ephemerides for OpenOrb from %s" %(depath))

    # Fix up the meta-kernel paths, if it's been templated
    tmpdir = None
    txt = open(spice_mk).read()
    if "{{dirname}}" in txt:
        datalink = os.path.dirname(spice_mk)
        if len(datalink) > 80:
            # SPICE has an 80-character limitation on string variables (sigh), so
            # create a symlink to data dir from a (hopefully) shorter path
            import tempfile, atexit
            tmpdir = tempfile.TemporaryDirectory()
            atexit.register(tmpdir.cleanup)		# clean up the dir on exit
            datalink = os.path.join(tmpdir.name, 'd')
            os.symlink(os.path.dirname(spice_mk), datalink)

        # fill out the template
        txt = txt.replace("{{dirname}}", datalink)
        with open("tmp_meta_kernel", "w") as fp:
            fp.write(txt)
        spice_mk = "tmp_meta_kernel"

    # Loading Spice Meta Kernel
    try:
        sp.furnsh(spice_mk)
    except Exception as e:
        print(e)
        sys.exit("Unable to load SPICE metakernel from %s" %(spice_mk))

    #Read in the stopTimesFile
    if stopTimesFile:
        df_stopTimes = pd.read_csv(stopTimesFile, header=None, comment="#", sep='\s+')
        # Combine year, month, day into a single string
        df_stopTimes['date_str'] = df_stopTimes.apply(lambda row: f"{int(row[1]):04d}-{int(row[2]):02d}-{int(row[3]):02d}", axis=1)
        df_stopTimes['decimal_day'] = df_stopTimes[3].apply(lambda x: x - int(x))  # Extract decimal part of the day


        # Function to convert calendar date to MJD using only spiceypy
        def calendar_to_mjd(date_str, decimal_day):
            # Convert calendar date to ephemeris time (ET)
            et = sp.str2et(f"{date_str}")
            et += decimal_day * 86400

            # Convert ET to Julian Date in UTC
            jd = 2451545.0 + et/86400
            mjd = jd - 2400000.5
            #jd_utc_str = sp.et2utc(et, 'J', 10)  # JD format
            #jd = float(jd_utc_str.split()[1])  # Extract the JD value from the string
            # Convert JD to MJD
            #mjd = jd - 2400000.5
            return mjd
        # Apply the MJD conversion, including the decimal day
        df_stopTimes['MJD'] = df_stopTimes.apply(lambda row: calendar_to_mjd(row['date_str'], row['decimal_day']), axis=1)
        
        # Convert DataFrame to dictionary with column 0 as the key and MJD as the value
        stopTimeDict = df_stopTimes.set_index(0)['MJD'].to_dict()
        
    t0 = time.time()
    if spaceflag=='T':
        spaceflag = True
        if not glob.glob(scspk):
            sys.exit("Couldn't find spacecraft SPK file.")
        else:
            shutil.copy2(scspk,os.getcwd())
    else:
        spaceflag = False
        try:
            obscodefile   = config.get('SURVEY','MPCobscode file')
        except:
            sys.exit('MPC Obs codes file not provided')
        obscodefile       = resolve_path(obscodefile)
        try:
            obscode       = config.get('SURVEY','Telescope')
        except:
            sys.exit('Observatory code not provided')
        # Loading the MPC list of observatory codes and coordinates
        b=ts.telescopelist(obscodefile)
        # Creating SPICE SPK for an observatory
        b.createspk(obscode,spkstart-10,spkstart+spkndays+10)

    a=ss.asteroidlist(population_model, asteroidspks, spkstart-10, spkstart+spkndays+100, spkstep, stopTimeDict, object1, nObjects)

    if makespks=='T':
        if (glob.glob("%s" %(asteroidspks+'*.bsp')) and not args.f):
            sys.exit("Some of the SPKs might already exist. Run with -f flag to overwrite. Alternatively set Make SPKs = F in configuration file or change the value of Asteroid SPKs in the config file.")
        else:
            os.makedirs(asteroidspkpath,exist_ok=True)
            a.generatestates(nbody, args.f)
            del a
            a=ss.asteroidlist(population_model,asteroidspks,spkstart-10, spkstart+spkndays+100, spkstep, stopTimeDict, object1,nObjects)    

    # Loading camera FOV definition
    c=ts.camera(cameradef_file,spiceik)
    # Loading list of pointings from survey and creating SPICE kernels
    c.createckfk(obscode, surveydb, starttime, Field1, nFields, spice_mk)

    # starttime and ndays that covers the timespan of the survey
    tmptimes=c.fieldMJD
    starttime=tmptimes[0]
    endtime=tmptimes[-1]
    ndays=np.ceil(endtime-starttime)
    if (ndays<0):
        sys.exit('nFields exceeds the number of frames in the database')

    print("Survey length:")
    print("Field 1 : ", starttime)
    print("Field n : ", endtime)
    print("Days : ", ndays)
    print('END HEADER')

    threshold=np.radians(threshold)
    t1 = time.time()
    a.simulate(starttime, starttime+ndays, c, threshold, obscode)
    t2 = time.time()
    print("#Startup time (s): ", (t1-t0), " Simulation time (s): ", (t2-t1))
    #os.system('rm ckip fakesclk test.fk tmp.fk camera.ti cksetupfile tmp')
