import numpy as np
import pyoorb as oo
import spiceypy as sp
import csv
import os
import sys
import pdb
import copy
import warnings
import time
from . import orbits
from . import ooephemerides
from . import shared

np.set_printoptions(precision=15)

############################################################################################################
#                                                Asteroid Class                                            #
############################################################################################################

class asteroids:

    def __init__(self, oorbElem, name):

        """Asteroids class stores all information relevant to an asteroid.
        Includes asteroid ID, elements, H & G values, state vectors, nightly states, etc.

        Parameters
        ----------

            oorbElem: array
                Array of elements and phys properties in the OpenOrb format
            name: string
                Users name/designation of the Asteroid 
        """
        self.id = int(oorbElem[0][0])
        self.name=name

        # Computing Spice ID from asteroid ID in list.
        self.spiceid=self.internal2spice(self.id)
        self.oorb_orbit = oorbElem

        self.ephcount=0

#-----------------------------------------------------------------------------------------------

    def internal2spice(self,id):

        """Computes a spice ID for an asteroid using its internal ID

        Parameters
        ----------
             id : int
                 Asteroid internal ID
        
        Returns
        -------
            spiceid : int
                 Asteroid ID for NAIF SPICE libraries
        
        """
        
        return(2444444+id)

#-----------------------------------------------------------------------------------------------

    def mpc2internal(self,Code):

        """Computes internal observatory code for SPICE
        
        """
        
        if (Code.lstrip('+-').isdigit()):
            internal_code=int(Code)
        else:
            internal_code=(ord(Code[0])-55)*100+int(Code[1:])

        internal_code = -internal_code
        return (internal_code)

#-----------------------------------------------------------------------------------------------
    
    def orbit2oorb(self, id, elements, elemtype, epoch, timescale, physprops):
        
        """ Store orbital elements in Open Orb standard.
        
        Open Orb elements are as follows:
        [id, q or a, e, i, upper case omega, lower case omega,
        t_p, element types, epoch, timescale for epoch, H, G]
        
        Element type: 2 - Cometary, 3 - Keplerian
        Timescale: 1 - UTC, 2 - UT1, 3 - TT, 4 - TAI

        (Davide: Use TT)

        Parameters
        ----------

            id : int
                Asteroid ID (not spice ID)
            elements : float array
                Orbital elements as described in docstring of __init__
            elemtype : int
                Type of orbital elements (cometary or Keplerian)
            epoch : float
                Epoch of elements (MJD)
            timescale : int
                Timescale of epoch (UTC, UT1, TT, or TAI) 
            physprops : float array
                [H, G]

        """

        oorb_elems = np.concatenate(([id], elements, [elemtype], [epoch],[timescale],physprops),0)
        oorb_elems = np.expand_dims(oorb_elems,0)
        
        return (oorb_elems)
        
#-----------------------------------------------------------------------------------------------

    def dis(self):
        
        """Placeholder to display asteroid properties

        """

        print (self.id, self.name)
        
#-----------------------------------------------------------------------------------------------

    def elem2vec(self):
        
        """Converting from orbital elements to state vectors

        """

        self.svec,err=oo.pyoorb.oorb_element_transformation(self.oorb_orbit,1)

#-----------------------------------------------------------------------------------------------
        
    def propagate1(self,epoch):

        """Propagate state vectors using N-body to epoch and store in memory

        """

        self.ephcount += 1
        tmpstart=np.expand_dims(self.svec[:][self.ephcount-1],0)
        tmp,err=oo.pyoorb.oorb_propagation(tmpstart,epoch, in_dynmodel='N')
        self.svec=np.append(self.svec,tmp,axis=0)

#-----------------------------------------------------------------------------------------------
        
    def propagate12b(self, epoch):

        """Propagate state vectors using 2-body to epoch and store in memory

        """

        self.ephcount += 1
        tmpstart=np.expand_dims(self.svec[:][self.ephcount-1],0)
        tmp,err=oo.pyoorb.oorb_propagation(tmpstart,epoch, in_dynmodel='2')
        self.svec=np.append(self.svec,tmp,axis=0)

#-----------------------------------------------------------------------------------------------
        
    def mjd2et(self, time):
        
        """Convert MJD to ET (seconds past J2000)

        """
        
        time=time+shared.mjd2jd
        return (sp.str2et('JD '+repr(time)))

#-----------------------------------------------------------------------------------------------
    
    def generateepochs(self,start,end,step):
        
        """Generate list of times using start and end times and number of steps

        Parameters
        ----------

            start : float
                Starting time (MJD)
            end : float
                Ending time (MJD)
                Same size as t
            step : float 
                Step size (days)  
         
        """

        start = float(start)
        end = float(end)
        step = float(step)
        self.step=step # Used in save module. Should get rid of this.
        n_times=1+(end-start)/step

        self.times = np.arange(start, end, step)
        self.timesjd=self.times+shared.mjd2jd
        self.timeset=np.zeros(len(self.times)+1)

        initialtime = self.svec[0][8]+shared.mjd2jd
        self.timeset[0]=sp.str2et('JD '+repr(initialtime))
        for i in range(1,len(self.timeset)):
            self.timeset[i]=sp.str2et('JD '+repr(self.timesjd[i-1]))

        # Formatting the times in the way openorb wants it
        twos=np.repeat(2,len(self.times))
        tmp=np.stack((self.times,twos))
        self.times=np.swapaxes(tmp,0,1)

        
#-----------------------------------------------------------------------------------------------
        
    def propagate(self, type, start,stop,step):

        #Type is the type of propagation
        #'F' is two-body
        #'T' is n-body
        self.generateepochs(start,stop,step)

        if type == 'T':
            
            for time in self.times:
                self.propagate1(time)

        if type == 'F':
            
            for time in self.times:
                self.propagate12b(time)

#-----------------------------------------------------------------------------------------------       

    def vismag(self, H, G, phase, delta, r):
        
        """ Compute the apparent visual magnitude
        Parameters
        ----------
            H : float
                Absolute magnitude
            G : float
                Slope parameter
            phase : float
                Phase angle (radians)
            delta : float
                Distance from the observer (km)
            r : float
                Distance from the Sun (km)

        Returns
        -------
            V : float
                Apparent visual magnitude of the asteroid
        
        """

        delta=delta/shared.au2km
        r=r/shared.au2km
        
        A=np.zeros(2)
        B=np.zeros(2)
        A[0] = 3.33
        A[1] = 1.87
        B[0] = 0.63
        B[1] = 1.22

        phi = np.exp(-A*np.power(np.tan(0.5*phase),B))
        phase_function = 2.5*np.log10((1-G)*phi[0] + G*phi[1])
        V = H + 5*np.log10(delta) + 5*np.log10(r) - phase_function
        return(V)
        
#-----------------------------------------------------------------------------------------------
            
    def save(self,base_fname, force):

        """ Save asteroid state vectors in a NAIF SPICE SPK file
         Base filename will be appended by asteroid internal ID + extension (.bsp)

        Parameters
        ----------
            base_fname: string
                Base file name. 
            force : boolean
                Force the program to rewrite asteroid SPK files
                
        """
        
        n_segments = 1 # This parameter should be an argument or should always be 1
        segmentsize = (len(self.svec)-1)//n_segments

        # File name to save SPK to
        self.spkname = base_fname+str(self.id)+".bsp"

        if force:
            if os.path.isfile(self.spkname):
                os.system("rm %s" %(self.spkname))
            
        handle = sp.spkopn(self.spkname,"spkfile",500)
        
        for i in np.arange(0,n_segments):

            # Computing first and last index for segment
            idx_first=i*segmentsize+1
            idx_last =(i+1)*segmentsize+1

            # first elements are initial values from initial orbit so ignoring those            
            segmenttimeet = self.timeset[idx_first:idx_last]
            segmenttime =   self.svec[idx_first:idx_last,8]
            segmentx =      self.svec[idx_first:idx_last,1]*shared.au2km
            segmenty =      self.svec[idx_first:idx_last,2]*shared.au2km
            segmentz =      self.svec[idx_first:idx_last,3]*shared.au2km
            segmentxd =     self.svec[idx_first:idx_last,4]*shared.au2km/shared.day2s
            segmentyd =     self.svec[idx_first:idx_last,5]*shared.au2km/shared.day2s
            segmentzd =     self.svec[idx_first:idx_last,6]*shared.au2km/shared.day2s
            
            tmpstates=[segmentx,segmenty,segmentz,segmentxd,segmentyd,segmentzd]
            spicestates=np.swapaxes(tmpstates,0,1)
            spicestates=spicestates.tolist()

            # Coordinate transformation from Ecliptic to Equatorial frame for saving as SPKs
            spicestateseq=copy.deepcopy(spicestates)
#            spicestateseq=spicestates[:]
            counter=0
            for state in spicestates:
                matrix=sp.sxform("ECLIPJ2000","J2000",segmenttimeet[counter])
                spicestateseq[counter]=sp.mxvg(matrix,spicestates[counter],6,6)
                counter=counter+1
                
#           SPK Type 9, Lagrange interpolation (Open Orb states are heliocentric, so central body is 10, which is the Sun)
#            sp.spkw09(handle, self.spiceid, 10, "J2000", segmenttimeet[0], segmenttimeet[-1], str(i), 10, segmentsize, spicestateseq, segmenttimeet)

#           SPK Type 5 (Two body interpolation)
            sp.spkw05(handle,2444444+self.id,10, "J2000", segmenttimeet[0], segmenttimeet[-1], str(i), shared.gms, segmentsize, spicestateseq, segmenttimeet)
            
        sp.spkcls(handle)
        
#-----------------------------------------------------------------------------------------------

    def nightlystates(self,timerange):
        
        """ Generate asteroid state at every midnight within the input time range
        
        Expects integer MJD dates (midnight). Generates geocentric coordinates
        (range, RA, DEC, rangerate, dRA, dDec) for the asteroid at each midnight
        between, and including, start and end.

        Parameters
        ----------
            timerange: 2 element float numpy array
                [Starting time, Ending time] (MJD)
                
        """

        # Converting from MJD to JD and computing number of days including endpoints
        start=int(timerange[0])+shared.mjd2jd
        stop =int(timerange[-1])+shared.mjd2jd
        ndays=int(stop-start)+1

        # Loading SPK
        try:
            sp.furnsh(self.spkname)
        except:
            sys.exit("Cannot load %s" %(self.spkname))

        # Initializing state variables
        self.nightlyt=np.zeros(ndays)
        self.nightlyr=np.zeros(ndays)
        self.nightlydr=np.zeros(ndays)
        self.nightlyra=np.zeros(ndays)
        self.nightlydra=np.zeros(ndays)
        self.nightlydec=np.zeros(ndays)
        self.nightlyddec=np.zeros(ndays)
 
        counter = -1

        for i in np.linspace(start,stop,num=ndays):

            counter += 1

            # Converting from UTC to ET
            timeet=sp.str2et('JD '+repr(i))

            # Finding direction of asteroids from geocenter
            # Finding states wrt station adds 4 days for 1e6 objects

            asteroidsvec=sp.spkezr(str(self.spiceid),timeet,"J2000","LT","Earth") # Units are km and km/s
            
            # Converting from Rectangular to Latitudinal coordinates.
            # tmp is a 6 element vector (R, Long, Lat, Dr, Dlong, Dlat) Units are radians and radians/s
            tmp=sp.xfmsta(asteroidsvec[0][0:6],"RECTANGULAR","LATITUDINAL"," ")

            self.nightlyt[counter]   = i
            self.nightlyr[counter]   = tmp[0]
            self.nightlydr[counter]  = tmp[3]
            self.nightlyra[counter]  = tmp[1]
            self.nightlydra[counter] = tmp[4]
            self.nightlydec[counter] = tmp[2]
            self.nightlyddec[counter]= tmp[5]

        # Unloading SPK
#        sp.unload(self.spkname)

#-----------------------------------------------------------------------------------------------

    def shortlist(self,camera,thresh_angle):

        """ Short list candidate frames for the asteroid

        Extrapolates asteroid position from nightly state to field time as:

        DEC = DEC_midnight+dDec*(field_time-midnight_time)
        RA = RA_midnight+dRA*cos(DEC_midnight)*(field_time-midnight_time)        

        Computes cosine of the anglular distance (A) between the field center (fra, fdec) and the asteroid (RA, DEC).
        Divides cosine of thresh_angle by cos(A) to assign a weight to the corresponding frame.
        Weight of 0 indicates that the asteroid is within threshold angular distance of the field center
        Weight of 1 or greater indicates that the asteroid is outside threshold angular distance of the field center.
        
        
        Parameters
        ----------
            camera : camera object

            thresh_angle : float
                Threshold angular distance (radians)

        Returns
        ---------
            field_times : float array 
                Array of candidate field times
            
            field_ids : float array
                Array of candidate field IDs
        
        """

        #Field times, id's, ra's, and dec's
        t=camera.fieldMJD
        ids=camera.obsHistID
        fra=camera.fieldRA
        fdec=camera.fieldDec

        # Cosine of threshold angle
        cos_thresh=np.cos(thresh_angle)

        # Computing indices of corresponding midnight states for each FOV
        indices=(t.astype(int)-t[0].astype(int)).astype(int)
        
        tt=t+shared.mjd2jd
        timediff = (tt - self.nightlyt[indices]) * shared.day2s

        # Declination extrapolated to field time
        deltadec = self.nightlyddec[indices] * timediff
        ftdec = self.nightlydec[indices] + deltadec

        # RA extrapolated to field time
        deltara = self.nightlydra[indices] * timediff
        ftra = self.nightlyra[indices] + deltara

        # Taking cosines and sines of necessary quantities
        cosfdec = np.cos(fdec)
        sinfdec = np.sin(fdec)
        cosftdec = np.cos(ftdec)
        sinftdec = np.sin(ftdec)
        radiff = fra-ftra
        cosradiff = np.cos(radiff)

        # Dot product between unit vectors pointing to asteroid and field center.
        cos_angle = cosfdec * cosftdec * cosradiff + sinfdec * sinftdec
        self.sep_weight = np.floor(cos_thresh/cos_angle)
        
        # Finding indices of candidate fields
#        self.candidatefields=np.where(self.sep_weight==0)[0]
        self.candidatefields=(self.sep_weight==0) | (self.nightlyr[indices] < 0.03*shared.au2km)

        # return a numpy array of candidate frame times and frame id's
        return(t[self.candidatefields], ids[self.candidatefields])

#-----------------------------------------------------------------------------------------------
    
    def checkvisspice(self,observer, camera, time, ids):

        """Use SPICE to check whether asteroid is in FOV at input times

        Parameters
        ----------
            observer : string
                Minor Planet Center Observatory Code
            camera : camera object

            time : float array
                Array of times (MJD UTC)

        """
        
        observerint=str(self.mpc2internal(observer))

        count = 0
        for i in range(0,len(time)):

            ttet= self.mjd2et(time[i])
            test=sp.fovtrg("-99999",str(self.spiceid),'POINT',' ','LT',observerint,ttet)
            if (test):# TESTING
                topostate,lttime = sp.spkezr(str(self.spiceid),ttet,"J2000","LT",observerint)
                radec=sp.xfmsta(topostate[0:6],"RECTANGULAR","LATITUDINAL"," ")
                astsunstate,tmp=sp.spkezr(str(self.spiceid),ttet-lttime,"J2000","None","10")
                
                #Computing phase angle
                phase_angle=np.dot(topostate[0:3],astsunstate[0:3])
                phase_angle=phase_angle/np.linalg.norm(topostate[0:3])/np.linalg.norm(astsunstate[0:3])
                phase_angle=np.arccos(phase_angle)

                #computing the visual magnitude
                astobsdelta=np.linalg.norm(topostate[0:3])
                astsundelta=np.linalg.norm(astsunstate[0:3])
                Vmag=self.vismag(self.oorb_orbit[0][10], self.oorb_orbit[0][11], phase_angle, astobsdelta, astsundelta)
                VH0=Vmag-self.oorb_orbit[0][10]

#                opstr="%-9d "  %(self.id)
                opstr="%-10s " %(self.name)
                opstr=opstr+"%9d " %(ids[i])
                opstr=opstr+"%12.6f " %(time[i])
#                opstr=opstr+"%10s " %(test)
                opstr=opstr+"%16.3f " %(radec[0])
                opstr=opstr+"%8.3f " %(radec[3])
                opstr=opstr+"%11.6f " %(np.degrees(radec[1])%360)
                opstr=opstr+"%9.6f " %(np.degrees(radec[4])*shared.day2s)
                opstr=opstr+"%10.6f " %(np.degrees(radec[2]))
                opstr=opstr+"%9.6f " %(np.degrees(radec[5])*shared.day2s)
                opstr=opstr+"%16.3f " %(astsunstate[0])
                opstr=opstr+"%16.3f " %(astsunstate[1])
                opstr=opstr+"%16.3f " %(astsunstate[2])
                opstr=opstr+"%11.6f " %(np.degrees(phase_angle))
                opstr=opstr+"%7.3f  " %(Vmag)
                opstr=opstr+"%7.3f  " %(VH0)
                print(opstr)
                count+=1
                # END TESTING

        sp.unload(self.spkname)
        
############################################################################################################
#                                            Asteroid List Class                                           #
############################################################################################################

class asteroidlist(asteroids):

    def __init__(self,inputfile,outputfile,object1,nObjects=-1):

        """ asteroidlist class contains a bunch of asteroid objects.
        
        Designed for convenience. 

        Parameters
        ----------
            inputfile : string
                File name of list of asteroid elements (Refer to SO.ssm from Grav et al. 2010)
            outputfile : string
                Base filename for NAIF SPICE SPK files of individual asteroids. 
                Base filename will be appended by asteroid ID, which is just the row number of the
                asteroid in the inputfile.
            object1 : int
                First object to be read from inputfile
            nObjects : int
                Number of objects to be read from inputfile

        """

        self.inputfile=inputfile
        self.outputfile=outputfile
        self.asteroids=[]

        #Initializing an orbits object and reading all the orbits to it.
        orbObj=orbits.Orbits()
        orbObj.readOrbits(inputfile)

        #Converting all the orbits from previous step into OpenOrb orbit format.
        ephemObj=ooephemerides.PyOrbEphemerides()
        ephemObj.setOrbits(orbObj)
        nOrbs=len(ephemObj.oorbElem)

        if nObjects == -1:
            nObjects = nOrbs-object1+1

        if (object1-1+nObjects > nOrbs):
            warnings.warn('nObjects (%d) goes beyond end of file (%s). Adjusting nObjects to %d.' %(nObjects, inputfile, nOrbs-object1+1), Warning)            
            nObjects=nOrbs-object1+1

        #Initializing asteroid objects for each orbit.
        for i in np.arange(object1-1,object1-1+nObjects):
            self.asteroids.append(asteroids([ephemObj.oorbElem[i]], orbObj.orbits.objId[i]))

        for i in self.asteroids:
            i.elem2vec()
            i.spkname=outputfile+str(i.id)+".bsp"
            
#-----------------------------------------------------------------------------------------------
            
    def generatestates(self,type,start,stop,step, force):

        """Generate state vectors for all asteroids contained in this object.

        Parameters
        ----------
            Type : float
                Type of propagation. 'F' - 2-body,
                                     'T' - n-body (8 planets and Sun)
            start : float
                Epoch for first set of state vectors (UTC MJD)
            stop : float
                Epoch for last set of state vectors (UTC MJD)
            step : float
                Step size for generating and storing state vectors (UTC MJD)
            force : boolean
                Force the program to rewrite asteroid SPK files
        """

        count=0

        # Propagate using openorb and save spks.
        while self.asteroids:
            i=self.asteroids[0]
            i.propagate(type,start,stop,step)
            i.save(self.outputfile, force)
            del i
            del self.asteroids[0]
            count=count+1
            
#-----------------------------------------------------------------------------------------------        

    def savestates(self, force):

        """Saving states of all asteroids in NAIF SPICE SPK format.

        Parameters
        ----------
            force : boolean
                Force the program to rewrite asteroid SPK files
        """
        
        for i in self.asteroids:
            i.save(self.outputfile, force)

#-----------------------------------------------------------------------------------------------

    def simulate(self, starttime, stoptime, camera, threshold, obscode):
        """
        """
#        loading all SPICE kernels required for simulation
        sp.furnsh(camera.ikfile)
        sp.furnsh(obscode+".bsp")
        sp.furnsh(camera.ckfile)
        sp.furnsh(camera.fkfile)
        sp.furnsh(camera.sclkfile)
        count=0

        #Print header
#        head="#AstID "
        head="ObjID "
        head=head+"FieldID "
        head=head+"FieldMJD "
        head=head+"AstRange(km) "
        head=head+"AstRangeRate(km/s) "
        head=head+"AstRA(deg) "
        head=head+"AstRARate(deg/day) "
        head=head+"AstDec(deg) "
        head=head+"AstDecRate(deg/day) "
        head=head+"Ast-Sun(J2000 x,y,z)(km) "
        head=head+"Sun-Ast-Obs(deg) "
        head=head+"V "
        head=head+"V(H=0) "
        print(head)
        
        while self.asteroids:
            i=self.asteroids[0]
            i.nightlystates([starttime, stoptime])
            [times,ids]= i.shortlist(camera,threshold)
            i.checkvisspice(obscode,camera,times,ids)
            del i
            del self.asteroids[0]
            count=count+1

        # Unloading all SPICE kernels required for simulation
        sp.unload(camera.ikfile)
        sp.unload(camera.ckfile)
        sp.unload(camera.fkfile)
        sp.unload(camera.sclkfile)
        sp.unload(obscode+".bsp")

#-----------------------------------------------------------------------------------------------

##################################################################################################
