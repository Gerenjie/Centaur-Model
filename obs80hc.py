#!/usr/bin/env python3
"""
annotates observations with heliocentric positions
"""

import collections
import sys

import numpy
import novas.compat as novas
import novas.compat.solsys as solsys
import novas.compat.eph_manager as eph_manager

__all__ = ["Heliocentric", "HCAnnotator"]

# Some constants
Rearth_km = 6378.1363
au_km = 149597870.700  # This is now a definition
Rearth_AU = Rearth_km / au_km

# fields copied from Optical/SpaceBased + heliocentric position
# of observing site or spacecraft
Heliocentric = collections.namedtuple('Heliocentric',
    'num desig disc note1 note2 jdutc ra dec mag band cod hx hy hz')
"""named tuple with fields in common with the obs80 named tuple types Optical
and SpaceBased, and annotated with heliocentric observing point coordinates.

Parameters
----------
num : string
    packed perm desig + comet/sat stuff
desig : string
    packed provisional or temporary desig
disc : string
    discovery asterisk
note1 : string
    overloaded note/program code
note2 : string
    overloaded record type/measurement technology
jdutc : float
    time of obs parsed to JD
ra : float
    RA parsed to hours
dec : float
    dec parsed to degrees
mag : float
    mag parsed to float
band : string
    mag band
cod : string
    3 character obs code
hx : float
    heliocentric x in AU
hy : float
    heliocentric y in AU
hz : float
    heliocentric z in AU
"""


def getEarthPosition(jd_tdb):
    """ The input time should be on the TDB time scale.
    TT - TDB is less than 2 milliseconds.
    Returns equatorial rectangular
    coordinates in AU referred to the ICRS.
    """
    pos, _ = solsys.solarsystem(jd_tdb, 3, 1)
    return pos


class HCAnnotator:
    """annotates observations with heliocentric positions.

    A JPL binary ephemeris file is required at the time the constructor is
    called.  There are apparently three ways to provide this file.

    1. Have the package novas_de405 installed and importable.
    2. Know the location of the binary ephemeris file and supply it as
       ephem_name to the constructor
    3. Know the location of the binary ephemeris file and set the environment
       variable EPHEMERIS_FILE.

    See readme, novas py readme and doc for novas eph_manager.ephem_open.

    Parameters
    ----------
    siteXYZ : dictionary
        dict as returned by module obscode, function siteXYZ.
    leapsec : object
        LeapSeconds object such as defined in module leapsec.
    ephem_name : string, optional
        File name of binary JPL ephemeris.
	
    """
    def __init__(self, siteXYZ, leapsec, ephem_name=None):
        
        eph_manager.ephem_open(ephem_name)
        self.siteXYZ = siteXYZ
        self.leapsec = leapsec
        self.cache = {}

    def annotate(self, parsed, other=None):
        """generator, annotates observations with heliocentric coordinates.

        Parameters
        ----------
        parsed : iterable
            iterable over objects including type Optical and SpaceBased,
            such as the generator obs80.read80
        other : list
            If a list is passed, objects other than Optical and SpaceBased
            are appended.  Pass None if the list is not of interest.

        Yields
        ------
        namedtuple
            Heliocentric named tuple for each Optical and SpaceBased
            observation, annotated with heliocentric positions.
        """
        for obs in parsed:
            t = type(obs).__name__
            if t == 'Optical':
                hx, hy, hz = self._xyzOpt(obs)
                n2 = obs.note2
                cod = obs.cod
            elif t == 'SpaceBased':
                hx, hy, hz = self._xyzSat(obs)
                n2 = 'S'
                cod = obs.cod
            elif t == 'Roving':
                hx, hy, hz = self._xyzRov(obs)
                n2 = 'V'
                cod = '247'
            else:
                try:
                    other.append(obs)
                except AttributeError:
                    pass
                continue
            yield Heliocentric(obs.num, obs.desig, obs.disc,
                obs.note1, n2, obs.jdutc, obs.ra, obs.dec,
                obs.mag, obs.band, cod, hx, hy, hz)

    def _xyzOpt(self, obs):
        """ calculates the heliocentric position of the observatory at UTC JD
        in equatorial cartesian coordinates.

        Note that for now it ignores polar motion, differences between TT and
        TDB, and differences between UTC and UT1.  The last is probably the
        biggest source of error.
        """
        key = (obs.cod, obs.jdutc)
        try:
            return self.cache[key]
        except KeyError:
            pass
        obsVec = self.siteXYZ[obs.cod]
        leaps = self.leapsec.getLeapSeconds(obs.jdutc)
        jd_tt = obs.jdutc + (32.184 + leaps) / (24.0 * 60 * 60)
        DUT1 = 0.0  # for now
        delta_t = 32.184 + leaps - DUT1
        jd_ut1 = obs.jdutc + DUT1
        xp, yp = 0.0, 0.0  # for now
        geocentric_vec = novas.ter2cel(jd_ut1, 0.0, delta_t, xp, yp, obsVec)
        geocentric_vec = Rearth_AU * numpy.array(geocentric_vec)
        pos = getEarthPosition(jd_tt)  # Might need to switch to tdb
        hc = (pos[0] + geocentric_vec[0],
            pos[1] + geocentric_vec[1], pos[2] + geocentric_vec[2])
        self.cache[key] = hc
        return hc

    def _xyzSat(self, obs):
        """ calculates the heliocentric position of an observing spacecraft at
        UTC JD in equatorial cartesian coordinates.
        """
        key = (obs.cod, obs.jdutc)
        try:
            return self.cache[key]
        except KeyError:
            pass
        leaps = self.leapsec.getLeapSeconds(obs.jdutc)
        jd_tt = obs.jdutc + (32.184 + leaps) / (24.0 * 60 * 60)
        pos = getEarthPosition(jd_tt)  # Might need to switch to tdb
        hc = pos[0] + obs.x, pos[1] + obs.y, pos[2] + obs.z
        self.cache[key] = hc
        return hc

    def _xyzRov(self, obs):
        """ calculates the heliocentric position of a roving observer at
        UTC JD in equatorial cartesian coordinates.
        """
        key = ('247', obs.jdutc)
        try:
            return self.cache[key]
        except KeyError:
            pass
        leaps = self.leapsec.getLeapSeconds(obs.jdutc)
        delta_t = 32.184 + leaps
        jd_tt = obs.jdutc + delta_t / (24.0 * 60 * 60)
        e = getEarthPosition(jd_tt)  # Might need to switch to tdb
        observer = novas.make_observer_on_surface(obs.lat, obs.lon, obs.alt,
            5, 990)
        # posvel returns GCRS
        o, _ = novas.geo_posvel(jd_tt, delta_t, observer)
        # add ICRS and GCRS?
        hc = (e[0] + o[0], e[1] + o[1], e[2] + o[2])
        self.cache[key] = hc
        return hc


if __name__ == '__main__':
    import pandas
    import obs80
    import obscode
    import leapsec

    try:
        hc = HCAnnotator(
            obscode.siteXYZ(obscode.read5c('obscode.dat')),
            leapsec.LeapSeconds('leap-seconds.list'))
        other = []
        flat = hc.annotate(obs80.parse80(open(sys.argv[1])), other)
        print(pandas.DataFrame(flat))
        print(len(other), "others")
        print("cache len:", len(hc.cache))
    except IndexError:
        print("usage: obs80hc <obs file>")
    except Exception as e:
        print(e)
