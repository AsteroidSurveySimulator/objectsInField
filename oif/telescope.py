import filecmp
import shutil
import numpy as np
import spiceypy as sp
import csv
import pandas as pd
import sqlite3
import os
import sys
import matplotlib.pyplot as plt
from . import shared

############################################################################################################
#                                           Telescope List Class                                           #
############################################################################################################

class telescopelist:

    def __init__(self, mpcobslist):

        """Initialize internal observatory database.
        Database contains MPC observatory code. 
        
        Parameters
        ----------

            mpcobslist : string
                File name of file containing the Minor Planet Center observatory codes list

        """

        # Creating a list of dictionaries. Each dictionary holds information for an observatory
        try:
            tmp=pd.read_fwf(mpcobslist,widths=[4,10,9,9])
        except:
            sys.exit("Unable to load MPC Observatory list: %s" %(mpcobslist))
            
        self.obsdict=list(tmp.T.to_dict().values())

        count=0
        for i in self.obsdict:

            # Computing cartesian coordinates of observatories in Earth frame            
            rcos=i['cos']
            rsin=i['sin']
            long=np.radians(i['Long.'])
            x=rcos*np.cos(long)*shared.REarth
            y=rcos*np.sin(long)*shared.REarth
            z=rsin*shared.REarth

            # Inserting cartesian x, y, z, into list of dictionaries.
            self.obsdict[count]['x']=x
            self.obsdict[count]['y']=y
            self.obsdict[count]['z']=z            

            count +=1

#-----------------------------------------------------------------------------------------------

    def mpc2internal(self,Code):

        """Computes internal observatory code
        
        """
        
        if (Code.isdigit()):
            internal_code=int(Code)

        else:
            internal_code=(ord(Code[0])-55)*100+int(Code[1:])

        internal_code = -internal_code
        return (internal_code)

#-----------------------------------------------------------------------------------------------

    def mjd2et(self, time):

        """Convert MJD to ET

        """

        time=time+2400000.5
        return (sp.str2et('JD '+repr(time)))

#-----------------------------------------------------------------------------------------------
    
    def createspk(self,Code, start, stop):
        
        """Create NAIF SPICE SPK file for Earth-based observatory

        Parameters
        ----------
            Code : string
                Minor Planet Center observatory code
            start : float
                Beginning time in SPK (MJD UTC)
            stop : float
                Ending time in SPK (MJD UTC)
        
        """
        
        obsid=self.mpc2internal(Code)
        startet=self.mjd2et(start)
        stopet=self.mjd2et(stop)

        # Search for Observatory in List of Dictionaries
        for i in self.obsdict:
            if i['Code'] == Code:
                selected_observer=i

                # Location in km
                x=i['x']/1000.0
                y=i['y']/1000.0
                z=i['z']/1000.0

        # Writing SPK file
        with open ("tmp","w") as f:
            f.write("\\begindata\n\n")
            f.write("SITES = (\'%s\')\n" %(Code))
            f.write("%s_FRAME = \'EARTH_FIXED\'\n" %(Code))
            f.write("%s_IDCODE = %i\n" %(Code, obsid))
            f.write("%s_XYZ = ( %f, %f, %f )\n" %(Code, x, y, z))
            f.write("%s_EPOCH = @2000-Jan-1/12:00\n" %(Code))
            f.write("%s_CENTER = 399\n" %(Code))
            f.write("%s_BOUNDS = ( %f, %f)\n" %(Code, startet, stopet))
        f.close()

        # Removing existing spk for same observatory
        os.system(f'rm -f {Code}.bsp') # Use python command for deleting files. shutil.rmtree(path)
        os.system(f'pinpoint -def tmp -spk {Code}.bsp > /dev/null')


############################################################################################################
#                                                Camera Class                                              #
############################################################################################################


class camera:


    def __init__(self, configfile, ikfile):

        """Initialize a camera object

        Parameters
        ----------
            configfile : string
                Name of file containing camera description
            ikfile : string
                name of NAIF SPICE Instrument Kernel
        
        """
        
        self.configfile = configfile
        self.ikfile=ikfile
        self.ckfile="test.ck"
        self.fkfile="test.fk"
        self.sclkfile="fakesclk"

        with open(configfile, newline='') as f:
            reader = csv.reader(f,delimiter=' ', skipinitialspace=True)
            type=next(reader)
            if (type[0].lower() == 'polygon'):
                self.x=[]
                self.y=[]
                self.z=[]
                n_lines = 0
                for row in reader:
                    n_lines=n_lines+1
                    tmpx=np.tan(np.radians(float(row[0])))
                    tmpy=np.tan(np.radians(float(row[1])))
                    tmpz=1
                    self.x.append(tmpx)
                    self.y.append(tmpy)
                    self.z.append(tmpz)
                self.save_poly(self.ikfile,[0,0,1],self.x,self.y,self.z)
                
            elif (type[0].lower() == 'circle'):
                self.x=float(next(reader)[0])
                self.save_circ(self.ikfile,[0,0,1],self.x)

            else:
                sys.exit("In file %s, instrument FOV is invalid. Only \"Circle\" or \"Polygon\" allowed." %(configfile))



#-----------------------------------------------------------------------------------------------

    def rotationmatrix(self, axis, angle):

        """Generate rotation matrix

        Parameters
        ----------
            axis : int
                Axis of rotation.
                1 - x-axis
                2 - y-axis
                3 - z-axis
            angle : float
                Angle of rotation (radians)

        Returns
        -------
            R : float array (3 x 3)
                Rotation Matrix            

        """        
        c=np.cos(angle)
        s=np.sin(angle)

        if axis==1:
            R = [[1,0,0],[0,c, -s], [0,s, c]]
        if axis==2:
            R = [[c,0,s],[0,1, 0], [-s,0, c]]
        if axis==3:
            R = [[c,-s,0],[s,c, 0], [0,0, 1]]
        return(R)

#-----------------------------------------------------------------------------------------------

    def mjd2et(self, time):

        """Convert MJD to ET

        """
        time=time+2400000.5
        return (sp.str2et('JD '+repr(time)))

#-----------------------------------------------------------------------------------------------
    
    def internal2spice(self,asteroid):
        
        """Computes a spice ID for an asteroid using its internal ID

        Parameters
        ----------
             asteroid : asteroid object
        
        Returns
        -------
            spiceid : string
                 Asteroid ID for NAIF SPICE libraries
        
        """
        return(str(2444444+asteroid.id))

#-----------------------------------------------------------------------------------------------
    
    def mpc2internal(self,Code):

        """Computes internal observatory code
        
        """
        if (Code.isdigit()):
            internal_code=int(Code)
        else:
            internal_code=(ord(Code[0])-55)*100+int(Code[1:])
        internal_code = -internal_code
        return (internal_code)

#-----------------------------------------------------------------------------------------------
    
    def save_poly(self,ikfile,boresight,x,y,z):

        """Saves a NAIF SPICE Instrument Kernel for a polygonal FOV

        Parameters
        ----------
             ikfile : string
                 Filename of instrument kernel
             boresight : float array (3)
                 Instrument boresight vector
             x : float array 
                 array of x coordinates of instrument corner vectors
             y : float array
                  array of y coordinates of instrument corner vectors
             z : float array
                  array of z coordinates of instrument corner vectors

        """        
        with open(ikfile, "w") as f:
            f.write("KPL/IK \nComments. \nMore comments. \nMore comments.\n\\begindata\n")
            f.write("INS-99999_FOV_SHAPE = 'POLYGON'\n")
            f.write("INS-99999_FOV_FRAME = 'CAMERA_FRAME'\n")
            f.write("INS-99999_BORESIGHT = (%15.10f %15.10f %15.10f)\n" %(boresight[0],boresight[1],boresight[2]))
            f.write("INS-99999_FOV_BOUNDARY_CORNERS = (%15.10f %15.10f %15.10f \n" %(x[0],y[0],z[0]))
            for i in range(1,len(x)-1):
                f.write("%34s %15.10f %15.10f %15.10f \n" %(" ",x[i],y[i],z[i]))
            f.write("%34s %15.10f %15.10f %15.10f) \n" %(" ",x[len(x)-1],y[len(y)-1],z[len(z)-1]))
        f.close()


#-----------------------------------------------------------------------------------------------
    
    def save_circ(self,ikfile,boresight,radius):

        """Saves a NAIF SPICE Instrument Kernel for a circular FOV

        Parameters
        ----------
             ikfile : string
                 Filename of instrument kernel
             boresight : float array (3)
                 Instrument boresight vector
             radius : float 
                 radius of instrument field of view

        """        
        with open(ikfile, "w") as f:
            f.write("KPL/IK \nComments. \nMore comments. \nMore comments.\n\\begindata\n")
            f.write("INS-99999_FOV_SHAPE = 'CIRCLE'\n")
            f.write("INS-99999_FOV_FRAME = 'CAMERA_FRAME'\n")
            f.write("INS-99999_BORESIGHT = (%15.10f %15.10f %15.10f)\n" %(boresight[0],boresight[1],boresight[2]))
            f.write("INS-99999_FOV_CLASS_SPEC = 'ANGLES' \n")
            f.write("INS-99999_FOV_REF_VECTOR = (0.0  1.0  0.0) \n")
            f.write("INS-99999_FOV_REF_ANGLE  = %f \n" %(radius))
            f.write("INS-99999_FOV_ANGLE_UNITS= 'DEGREES' \n")
        f.close()

#-----------------------------------------------------------------------------------------------
        
    def readfields(self, dbname, line1, nlines, startdate, surveydbquery):

        """Reads fields from database

        Parameters
        ----------
             dbname : string
                 Name of Survey Sim database (LSST opsim format)
             line1  : int
                 First line to be read from the database (starting from 1)
             nlines : int
                 Number of lines to read from the database
             startdate : float
                 Shift all field times such that first one is startdate (0 if default)  
             surveydbquery : string
                 customize the sql query string to account for changes in the opsim database designations 

        """        

        conn=sqlite3.connect(dbname)
        c=conn.cursor()
        self.obsHistID=np.zeros(nlines)
        self.fieldMJD=np.zeros(nlines)
        self.fieldRA=np.zeros(nlines)
        self.fieldDec=np.zeros(nlines)
        self.rotSkyPos=np.zeros(nlines)

        count=0
#        exec_str='SELECT obsHistID,expMJD,fieldRA,fieldDec,rotSkyPos FROM Summary order by expMJD limit %d,%d' %(line1-1,nlines)
#        exec_str='SELECT observationId,observationStartMJD,ra,dec,angle FROM ObsHistory order by observationStartMJD limit %d,%d' %(line1-1,nlines)
        exec_str=surveydbquery+' limit %d,%d' %(line1-1,nlines)
        
        for row in c.execute(exec_str):
            self.obsHistID[count] = row[0]
            self.fieldMJD[count]  = row[1]
            self.fieldRA[count]   = np.deg2rad(row[2])
            self.fieldDec[count]  = np.deg2rad(row[3])
            self.rotSkyPos[count] = np.deg2rad(row[4])
            count +=1

        # startdate is 0 if not provided by user. In this case use the default MJDs.
        if (startdate > 1):
            self.fieldMJD=self.fieldMJD+(int(startdate)-int(self.fieldMJD[0]))
            
#-----------------------------------------------------------------------------------------------

    def createckfk(self, observer, dbname, t0, field1, nfields, mk, surveydbquery):

        """Creates NAIF SPICE FK frame and corresponding CK frame

        Parameters
        ----------
             observer : string
                 MPC observatory code
             dbname : string
                 Name of database (opsim format) 
             t0 : float
                 Time of first field (0 if we want to use default)
             field1 : int
                 First field to be read from the database
             nfields: int
                 Number of fields to be read from the database
             mk : string
                 Name of SPICE meta kernel
             surveydbquery : string
                 customize the SQL query string to account for changes in the opsim database designations 

        """      

        observerint=self.mpc2internal(observer)
        instrumentint=observerint*1000

        with open("cksetupfile", "w") as f:
            f.write("KPL/IK \nComments describing the keywords and values \nto follow, as well as any other pertinent \ninformation.\n\\begindata\n")
            f.write("LSK_FILE_NAME           = '%s'\n" %(mk))
            f.write("\n")
            f.write("INTERNAL_FILE_NAME      = 'Survey Sim Camera Orientation'\n")
            f.write("\n")
            f.write(f"MAKE_FAKE_SCLK          = 'tmpsclk'\n")
            f.write("CK_TYPE                 = 3\n")
            f.write("CK_SEGMENT_ID           = 'Instrument Orientation'\n")
            f.write("INSTRUMENT_ID           = %i \n" %(instrumentint))
            f.write("REFERENCE_FRAME_NAME    = 'J2000'\n")
            f.write("ANGULAR_RATE_PRESENT    = 'NO'\n")
            f.write("\n")
            f.write("INPUT_DATA_TYPE         = 'SPICE QUATERNIONS'\n")
            f.write("INPUT_TIME_TYPE         = 'UTC'\n")
            f.write("MAXIMUM_VALID_INTERVAL  = 60\n") 
            f.write("\n")
            f.write("PRODUCER_ID             = 'Survey Sim, JPL'\n")
            f.write("\\begintext")
        f.close()


        self.readfields(dbname,field1,nfields, t0, surveydbquery)
        with open("ckip","w") as f:

            for i in range(len(self.fieldRA)):
                quat=self.computerotmat(self.fieldRA[i], self.fieldDec[i], self.rotSkyPos[i])

                #This helps with duplicate entries. For example enigma_1189 can have same fieldID's under different propID's
                #Issue warning for duplicate time. Have a verbose mode for displaying that (true as default)
                if (self.fieldMJD[i] !=self.fieldMJD[i-1]):
                    JD=self.fieldMJD[i]+shared.mjd2jd
                    timestring= 'JD'+repr(JD)
                    f.write("%s %f %f %f %f\n" %(timestring,quat[0],quat[1],quat[2],quat[3]))
        f.close()
        try:
            os.system('rm -f tmp.ck tmpsclk test.ck fakesclk')
        except:
            pass

        os.system(f'msopck cksetupfile ckip tmp.ck >/dev/null ')

        os.system('rsync tmpsclk fakesclk > /dev/null')
        os.system('rsync tmp.ck test.ck > /dev/null')

        with open("tmp.fk","w") as f:
            f.write("\\begindata\n\n")
            f.write("FRAME_CAMERA_FRAME = %i\n" %(instrumentint))
            f.write("FRAME_%i_NAME = 'CAMERA_FRAME'\n" %(instrumentint))
            f.write("FRAME_%i_CLASS = 3\n" %(instrumentint))
            f.write("FRAME_%i_CLASS_ID = %i\n" %(instrumentint, instrumentint))
            f.write("FRAME_%i_CENTER = %i\n" %(instrumentint, observerint))
            f.write("CK_%i_SCLK = %i\n" %(instrumentint, observerint))
            f.write("CK_%i_SPK = %i\n\n" %(instrumentint, observerint))
            f.write("\\begintext\n")
        f.close()
   
        os.system('rsync tmp.fk test.fk')

#-----------------------------------------------------------------------------------------------
        
    def computerotmat(self, ra, dec, pa):

        """Computes a rotation matrix and corresponding SPICE quaternion
        to point intrument from default z-axis pointing to field pointing

        Parameters
        ----------
             ra : float
                 Right Ascension of field center (radians)
             dec : float
                 Declination of field center (radians)
             pa : float
                 Position angle of instrument (radians)
        """      

        #Rotating from equatorial J2000 to camera frame
        #PA is the angle from the north axis to the instrument y axis, measured towards east.
        #Note: We are rotating the coordinate frame axes and not points, so the angles of
        #      rotation passed to the function rotationmatrix are negatives.
        rotmat0=self.rotationmatrix(3,-ra-np.pi/2.0)
        rotmat1=self.rotationmatrix(1,dec-np.pi/2.0)
        rotmat2=self.rotationmatrix(3,pa)
        rotmat=np.dot(rotmat1,rotmat0)
        rotmat=np.dot(rotmat2,rotmat)
        quat=sp.m2q(rotmat)
        return (quat)

##################################################################################################

