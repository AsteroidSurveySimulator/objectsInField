#!/opt/local/bin/python3.4
###########################################
# Asteroid Survey Simulation Software
# Developed at the Jet Propulsion Laboratory
# Authors: Shantanu Naidu
#          Davide Farnocchia
#          Steve Chesley
# Date:    Apr 06, 2018
##########################################


import context
import os
import code.shared as shared
import code.telescope as ts
import numpy as np
import spiceypy as sp
import time
import configparser
import glob
import sys
import argparse
import warnings
import pyoorb as oo
import code.sso as ss


#Parsing command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input  configuration filename.", type=str)
parser.add_argument("-f", help="Force the program to replace asteroid SPKs", action="store_true")
args = parser.parse_args()

inputfile=args.input


configdict={}

config = configparser.SafeConfigParser({'Make SPKs':'F','nDays':'4000','SPK step':'20','nbody':'F','Survey T0':'0','Field1':'1','nFields':'1000','Space':'F','n_Proc':'1','Threshold':'5','SPICE IK':'camera.ti','Object1':'1','nObjects':'1','Make SPKs':'T'})
config.readfp(open(inputfile))

#CONF section
try:
    spice_mk          = config.get('CONF','SPICE metakernel')
except:
    sys.exit('SPICE meta kernel not provided')
try:
    planetary_ephem   = config.get('CONF','Planetary ephem')
except:
    sys.exit('Planetary ephemeris file not provided')
try:
    data_path         = config.get('CONF','Data path')
except:
    sys.exit('Data directory not provided')
nProc                 = config.getint('CONF','n_Proc')

#CAMERA section
try:
    cameradef_file    = config.get('CAMERA','Camera')
except:
    sys.exit('Camera FOV definition file not provided')

threshold             = config.getfloat('CAMERA','Threshold')
if (threshold > 90.0):
    warnings.warn('Threshold was > 90 degrees. Setting it to 90',Warning)
    threshold=90
    
spiceik               = config.get('CAMERA','SPICE IK')
                    
#ASTEROID section
try:
    population_model  = config.get('ASTEROID','Population model')
except:
    sys.exit('Population model file not provided')
object1               = config.getint('ASTEROID','Object1')
if (object1==0):
    sys.exit('Object1 should >=0')
nObjects              = config.getint('ASTEROID','nObjects')
try:
    asteroidspks      = config.get('ASTEROID','Asteroid SPKs')
except:
    sys.exit('Asteroid SPK file basename not provided')
try:
    asteroidspkpath   = config.get('ASTEROID','Asteroid SPK path')
except:
    sys.exit('Directory containing asteroid SPKs not provided')
makespks          = config.get('ASTEROID','Make SPKs')
try:
    spkstart          = config.getint('ASTEROID','SPK T0')
except:
    sys.exit('Start time for creating SPK(s) not known')
try:
    spkndays          = config.getfloat('ASTEROID','nDays')
except:
    sys.exit('Number of days for creating SPK not known')
spkstep               = config.getfloat('ASTEROID','SPK step')
if (spkndays/spkstep < 11 and makespks == 'T'):
    spkstep=spkndays/11
    warnings.warn('Not enough steps to create SPKs. Reducing SPK step to %f' %(spkstep),Warning)
nbody                = config.get('ASTEROID','nbody')

#SURVEY section
try:
    surveydb          = config.get('SURVEY','Survey database')
except:
    sys.exit('Survey database not provided')
starttime             = config.getint('SURVEY','Survey T0')
Field1                = config.getint('SURVEY','Field1')
nFields               = config.getint('SURVEY','nFields')
spaceflag             = config.get('SURVEY','Space')
try:
    asteroidspks      = os.path.join(asteroidspkpath,asteroidspks)
except:
    sys.exit('Asteroid SPK directory unknown')
    
spice_mk              = os.path.join(data_path,spice_mk)
cameradef_file        = os.path.join(data_path,cameradef_file)
spiceik               = os.path.join(data_path,spiceik)
population_model      = os.path.join(data_path,population_model)
surveydb              = os.path.join(data_path,surveydb)
# Done reading configuration file

#If it made it this far, print header
with open(inputfile,'r') as f:
    for row in f:
        if(not row.startswith("#") and not row.startswith(";") and row.strip()):
            print(row,end='')

# Loading Spice Meta Kernel
try:
    sp.furnsh(spice_mk)
except:
    sys.exit("Unable to load SPICE metakernel from %s" %(spice_mk))

# Loading files for OpenOrb and SPICE
depath=os.path.join(data_path,planetary_ephem)
if os.path.isfile(depath):
    oo.pyoorb.oorb_init(ephemeris_fname=depath)
else:
    sys.exit("Unable to load planetary ephemerides for OpenOrb from %s" %(depath))

# Changing directory to data path
os.chdir(data_path)

if spaceflag=='T':
    spaceflag = True
    obscode       = config.get('SURVEY','SCID')
    if not glob.glob(obscode+'.bsp'):
        sys.exit("Couldn't find spacecraft SPK file.")

else:
    spaceflag = False
    try:
        obscodefile   = config.get('SURVEY','MPCobscode file')
    except:
        sys.exit('MPC Obs codes file not provided')
    obscodefile       = os.path.join(data_path,obscodefile)
    try:
        obscode       = config.get('SURVEY','Telescope')
    except:
        sys.exit('Observatory code not provided')
    # Loading the MPC list of observatory codes and coordinates
    b=ts.telescopelist(obscodefile)
    # Creating SPICE SPK for an observatory
    b.createspk(obscode,spkstart-10,spkstart+spkndays+10)
    
a=ss.asteroidlist(population_model, asteroidspks, object1, nObjects)

if makespks=='T':
    if (glob.glob("%s" %(asteroidspks+'*.bsp')) and not args.f):
        sys.exit("Some of the SPKs might already exist. Run with -f flag to overwrite. Alternatively set Make SPKs = F in configuration file or change the value of Asteroid SPKs in the config file.")
    else:
        os.makedirs(asteroidspkpath,exist_ok=True)
        a.generatestates(nbody,spkstart-10, spkstart+spkndays+100,spkstep, args.f)
        del a
        a=ss.asteroidlist(population_model,asteroidspks,object1,nObjects)    

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
t0 = time.time()
a.simulate(starttime, starttime+ndays, c, threshold, obscode)
t1 = time.time()
print("Simulation time: ", (t1-t0))
os.system('rm ckip fakesclk test.fk tmp.fk camera.ti cksetupfile tmp')
