from __future__ import print_function
import os
from itertools import repeat
import warnings
import numpy as np
import pandas as pd
import pyoorb as oo
from .orbits import Orbits

import time

__all__ = ['PyOrbEphemerides']


def dtime(time_prev):
    return (time.time() - time_prev, time.time())


class PyOrbEphemerides(object):
    """Generate ephemerides and propagate orbits using the python interface to Oorb.
    Inherits from Orbits and uses parent class to set orbital parameters.
    """
    def __init__(self, ephfile=None):
        # Set translation from timescale to OpenOrb numerical representation.
        # Note all orbits are assumed to be in TT timescale.
        # Also, all dates are expected to be in MJD.
        self.timeScales = {'UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4}
        self.elemType = {'CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5}

        # Set up oorb. Call this once.
        if ephfile is None and 'OORB_DATA' in os.environ:
            ephfile = os.path.join(os.getenv('OORB_DATA'), 'de405.dat')
        self.ephfile = ephfile
        self._init_oorb()
        self.oorbElem = None
        self.orb_format = None

    def _init_oorb(self):
        oo.pyoorb.oorb_init(ephemeris_fname=self.ephfile)

    def setOrbits(self, orbitObj):
        """Set the orbits, to be used to generate ephemerides.

        Immediately calls self._convertOorbElem to translate to the 'packed' oorb format.

        Parameters
        ----------
        orbitObj : Orbits
           The orbits to use to generate ephemerides.
        """
        if len(orbitObj) == 0:
            raise ValueError('There are no orbits in the Orbit object.')
        self._convertToOorbElem(orbitObj.orbits, orbitObj.orb_format)

    def _convertToOorbElem(self, orbitDataframe, orb_format):
        """Convert orbital elements into the numpy fortran-format array OpenOrb requires.

        The OpenOrb element format is a single array with elemenets:
        0 : orbitId (cannot be a string)
        1-6 : orbital elements, using radians for angles
        7 : element 'type' code (1 = CART, 2 = COM, 3 = KEP, 4 = DELauny, 5 = EQX (equinoctial))
        8 : epoch
        9 : timescale for epoch (1 = UTC, 2 = UT1, 3 = TT, 4 = TAI : always assumes TT)
        10 : magHv
        11 : g

        Sets self.oorbElem, the orbit parameters in an array formatted for OpenOrb.
        """
        oorbElem = np.zeros([len(orbitDataframe), 12], dtype=np.double, order='F')
        # Put in simple values for objid, or add method to test if any objId is a string.
        oorbElem[:,0] = np.arange(0, len(orbitDataframe), dtype=int) + 1
        # Add the appropriate element and epoch types:
        oorbElem[:,7] = np.zeros(len(orbitDataframe), float) + self.elemType[orb_format]
        oorbElem[:,9] = np.zeros(len(orbitDataframe), float) + self.timeScales['TT']
        # Convert other elements INCLUDING converting inclination, node, argperi to RADIANS
        if orb_format == 'KEP':
            oorbElem[:, 1] = orbitDataframe['a']
            oorbElem[:, 2] = orbitDataframe['e']
            oorbElem[:, 3] = np.radians(orbitDataframe['inc'])
            oorbElem[:, 4] = np.radians(orbitDataframe['Omega'])
            oorbElem[:, 5] = np.radians(orbitDataframe['argPeri'])
            oorbElem[:, 6] = np.radians(orbitDataframe['meanAnomaly'])
        elif orb_format == 'COM':
            oorbElem[:, 1] = orbitDataframe['q']
            oorbElem[:, 2] = orbitDataframe['e']
            oorbElem[:, 3] = np.radians(orbitDataframe['inc'])
            oorbElem[:, 4] = np.radians(orbitDataframe['Omega'])
            oorbElem[:, 5] = np.radians(orbitDataframe['argPeri'])
            oorbElem[:, 6] = orbitDataframe['tPeri']
        elif orb_format == 'CART':
            oorbElem[:, 1] = orbitDataframe['x']
            oorbElem[:, 2] = orbitDataframe['y']
            oorbElem[:, 3] = orbitDataframe['z']
            oorbElem[:, 4] = orbitDataframe['xdot']
            oorbElem[:, 5] = orbitDataframe['ydot']
            oorbElem[:, 6] = orbitDataframe['zdot']
        else:
            raise ValueError('Unknown orbit format %s: should be COM, KEP or CART.' % orb_format)
        oorbElem[:,8] = orbitDataframe['epoch']
        oorbElem[:,10] = orbitDataframe['H']
        oorbElem[:,11] = orbitDataframe['g']
        self.oorbElem = oorbElem
        self.orb_format = orb_format

    def _convertFromOorbElem(self, oorbElem):
        """Translate pyoorb-style orbital element array back into dataframe.

        Parameters
        ----------
        oorbElem : numpy.ndarray
            The orbital elements in OpenOrb format.

        Returns
        -------
        Orbits
            A new Orbits instance.
        """
        if self.orb_format == 'KEP':
            newOrbits = pd.DataFrame(self.oorbElem, columns=['objId', 'a', 'e', 'inc', 'Omega', 'argPeri',
                                                             'meanAnomaly', 'elem_type', 'epoch',
                                                             'epoch_type',
                                                             'H', 'g'])
            newOrbits['meanAnomaly'] = np.degrees(newOrbits['meanAnomaly'])
        elif self.orb_format == 'COM':
            newOrbits = pd.DataFrame(self.oorbElem, columns=['objId', 'q', 'e', 'inc', 'Omega', 'argPeri',
                                                             'tPeri', 'elem_type', 'epoch', 'epoch_type',
                                                             'H', 'g'])
        elif self.orb_format == 'CART':
            newOrbits = pd.DataFrame(self.oorbElem, columns = ['objId', 'x', 'y', 'z',
                                                               'xdot', 'ydot', 'zdot', 'elem_type', 'epoch',
                                                               'epoch_type', 'H', 'g'])
        else:
            raise ValueError('Unknown orbit format %s: should be COM, KEP or CART.' % self.orb_format)
        # Convert from radians to degrees.
        if self.orb_format == 'KEP' or self.orb_format =='COM':
            newOrbits['inc'] = np.degrees(newOrbits['inc'])
            newOrbits['Omega'] = np.degrees(newOrbits['Omega'])
            newOrbits['argPeri'] = np.degrees(newOrbits['argPeri'])
        # Drop columns we don't need and don't include in our standard columns.
        del newOrbits['elem_type']
        del newOrbits['epoch_type']
        # Have to swap orbit ids back to original values in Orbits object itself.
        newOrb = Orbits()
        newOrb.setOrbits(newOrbits)
        return newOrb

    def convertOrbitFormat(self, orb_format='CART'):
        """Convert orbital elements from the format in orbitObj into 'format'.

        Parameters
        ----------
        format : str, opt
            Format to convert orbital elements into.

        Returns
        -------
        """
        oorbElem, err = oo.pyoorb.oorb_element_transformation(in_orbits=self.oorbElem,
                                                              in_element_type=self.elemType[orb_format])
        if err != 0:
            raise RuntimeError('Oorb returned error %s' % (err))
        self.oorbElem = oorbElem
        self.orb_format = orb_format
        return

    def _convertTimes(self, times, timeScale='UTC'):
        """Generate an oorb-format array of the times desired for the ephemeris generation.

        Parameters
        ----------
        times : numpy.ndarray or float
            The ephemeris times (MJD) desired
        timeScale : str, optional
            The timescale (UTC, UT1, TT, TAI) of the ephemeris MJD values. Default = UTC, MJD.

        Returns
        -------
        numpy.ndarray
            The oorb-formatted 'ephTimes' array.
        """
        if isinstance(times, float):
            times = np.array([times])
        if len(times) == 0:
            raise ValueError('Got zero times to convert for OpenOrb')
        ephTimes = np.array(list(zip(times, repeat(self.timeScales[timeScale], len(times)))),
                            dtype='double', order='F')
        return ephTimes

    def _generateOorbEphs(self, ephTimes, obscode='I11'):
        """Generate ephemerides using OOrb (n-body).

        Parameters
        ----------
        ephtimes : numpy.ndarray
            Ephemeris times in oorb format (see self.convertTimes)
        obscode : int or str, optional
            The observatory code for ephemeris generation. Default=I11 (Cerro Pachon).

        Returns
        -------
        numpy.ndarray
            The oorb-formatted ephemeris array.
        """
        oorbEphems, err = oo.pyoorb.oorb_ephemeris(in_orbits=self.oorbElem, in_obscode=obscode,
                                                   in_date_ephems=ephTimes)
        if err != 0:
            raise RuntimeError('Oorb returned error %s' % (err))
        return oorbEphems

    def _generateOorbEphs2body(self, ephTimes, obscode='I11'):
        """Generate ephemerides using OOrb with two body mode.

        Parameters
        ----------
        ephtimes : numpy.ndarray
            Ephemeris times in oorb format (see self.convertTimes).
        obscode : int or str, optional
            The observatory code for ephemeris generation. Default=I11 (Cerro Pachon).

        Returns
        -------
        numpy.ndarray
            The oorb-formatted ephemeris array.
        """
        oorbEphems, err = oo.pyoorb.oorb_ephemeris_2b(in_orbits=self.oorbElem, in_obscode=obscode,
                                                      in_date_ephems=ephTimes)
        if err != 0:
            raise RuntimeError('Oorb returned error %s' % (err))
        return oorbEphems

    def _convertOorbEphs(self, oorbEphs, byObject=True):
        """Converts oorb ephemeris array to pandas dataframe, with labeled columns.

        The oorb ephemeris array is a 3-d array organized as: (object / times / eph@time)
        [objid][time][ephemeris information @ that time] with ephemeris elements
        0 : distance (geocentric distance)
        1 : ra (deg)
        2 : dec (deg)
        3 : mag
        4 : ephem mjd
        5 : ephem mjd timescale
        6 : dra/dt (deg/day) sky motion
        7 : ddec/dt (deg/day) sky motion
        8 : phase angle (deg)
        9 : solar elongation angle (deg)

        Here we convert to a numpy recarray, grouped either by object (default)
        or by time (if byObject=False).
        The resulting numpy recarray is composed of columns (of each ephemeris element),
        where each column is 2-d array with first axes either 'object' or 'time'.
        - if byObject = True : [ephemeris elements][object][time]
        (i.e. the 'ra' column = 2-d array, where the [0] axis (length) equals the number of ephTimes)
        - if byObject = False : [ephemeris elements][time][object]
        (i.e. the 'ra' column = 2-d arrays, where the [0] axis (length) equals the number of objects)

        Parameters
        ----------
        oorbEphs : numpy.ndarray
            The oorb-formatted ephemeris values
        byObject : boolean, optional
            If True (default), resulting converted ephemerides are grouped by object.
            If False, resulting converted ephemerides are grouped by time.

        Returns
        -------
        numpy.recarray
            The re-arranged ephemeris values, in a 3-d array.
        """
        ephs = np.swapaxes(oorbEphs, 2, 0)
        velocity = np.sqrt(ephs[6]**2 + ephs[7]**2)
        if byObject:
            ephs = np.swapaxes(ephs, 2, 1)
            velocity = np.swapaxes(velocity, 1, 0)
        # Create a numpy recarray.
        ephs = np.rec.fromarrays([ephs[0], ephs[1], ephs[2], ephs[3], ephs[4],
                                  ephs[6], ephs[7], ephs[8], ephs[9], velocity],
                                 names=['delta', 'ra', 'dec', 'magV', 'time', 'dradt',
                                        'ddecdt', 'phase', 'solarelon', 'velocity'])
        return ephs

    def generateEphemerides(self, times, timeScale='UTC', obscode='I11', byObject=True,
                            verbose=False):
        """Calculate ephemerides for all orbits at times `times`.

        This is a public method, wrapping self._convertTimes, self._generateOorbEphs
        and self._convertOorbEphs (which include dealing with oorb-formatting of arrays).

        The return ephemerides are in a numpy recarray, with axes
        - if byObject = True : [ephemeris values][object][@time]
        (i.e. the 'ra' column = 2-d array, where the [0] axis (length) equals the number of ephTimes)
        - if byObject = False : [ephemeris values][time][@object]
        (i.e. the 'ra' column = 2-d arrays, where the [0] axis (length) equals the number of objects)

        The ephemeris values returned to the user (== columns of the recarray) are:
        ['delta', 'ra', 'dec', 'magV', 'time', 'dradt', 'ddecdt', 'phase', 'solarelon', 'velocity']
        where positions/angles are all in degrees, velocities are deg/day, and delta is the
        distance between the Earth and the object in AU.

        Parameters
        ----------
        ephtimes : numpy.ndarray
            Ephemeris times in oorb format (see self.convertTimes)
        obscode : int, optional
            The observatory code for ephemeris generation. Default=807 (Cerro Tololo).
        byObject : boolean, optional
            If True (default), resulting converted ephemerides are grouped by object.
            If False, resulting converted ephemerides are grouped by time.
        verbose: boolean, optional
            If True, prints time required to calculate ephemerides. Default is False.

        Returns
        -------
        numpy.ndarray
            The ephemeris values, organized as chosen by the user.
        """
        t = time.time()
        ephTimes = self._convertTimes(times, timeScale=timeScale)
        oorbEphs = self._generateOorbEphs(ephTimes, obscode=obscode)
        ephs = self._convertOorbEphs(oorbEphs, byObject=byObject)
        dt, t = dtime(t)
        if verbose:
            print("# Calculating ephemerides for %d objects over %d times required %f seconds"
                  % (len(self.oorbElem), len(times), dt))
        return ephs

    def propagateOrbits(self, newEpoch):
        """Propagate orbits from self.orbits.epoch to new epoch (MJD TT).

        Parameters
        ----------
        new_epoch : float
            MJD TT time for new epoch.

        Returns
        -------
        PyOrbEphemerides
            New PyOrbEphemerides object, containing updated orbital elements for orbits specified by 'sso'.
        """
        newEpoch = self._convertTimes(newEpoch, timeScale='TT')
        old_orb_format = self.orb_format
        # COM format seems to crash propagation, so don't use that.
        if self.orb_format == 'COM':
            warnings.warn('Converting to CARTESIAN format elements')
            self.convertOrbitFormat(orb_format='CART')
        newOorbElem, err = oo.pyoorb.oorb_propagation_nb(in_orbits=self.oorbElem, in_epoch=newEpoch)
        if err != 0:
            raise RuntimeError('Orbit propagation returned error %d' % err)
        self.oorbElem = newOorbElem
        # Convert back to old format if necessary.
        if old_orb_format != self.orb_format:
            self.convertOrbitFormat(orb_format=old_orb_format)
        return
