"""
    
--------------------------------------------------------------

Oct 2018

Payne


Intended classes.py to contain classes used by;

(i) functions in mpcutilities

(ii) functions used across a broad range of MPC-packages
      [when they have no obvious better home-package]

 Note that these are split into

(i) kepcart-related (conversion routines)

 I have currently commented-out the following classes,
 because they relate to undeveloped/unported functionality

(ii) detection-related

(iii) orbit-related

--------------------------------------------------------------

"""

# Import third-party packages
# --------------------------------------------------------------
import numpy as np
from ctypes import *
import healpy as hp
import sys
import collections
import random

# Import neighboring packages
# --------------------------------------------------------------
import MPC_library as MPCL
import phys_const as PHYS
from ele220 import Ele220
import kepcart as kc
from universal_kepler import universal_kepler as UK
import params
import angles as ANG


# (i) Classes related to kepcart (coordinate conversion)
# --------------------------------------------------------------

class State(Structure):
    """
        ctypes "CartState" Structure for passing cartesian
        coords back and forth between python & C
        """
    _fields_ = [
                ('x', c_double),
                ('y', c_double),
                ('z', c_double),
                ('xd', c_double),
                ('yd', c_double),
                ('zd', c_double)
                ]
    def get_xyz(self):
        return self.x,self.y,self.z
    def get_uvw(self):
        return self.xd,self.yd,self.zd
    def ec_to_eq(self):
        self.x,self.y,self.z = PHYS.rotate_ec_to_eq(self.get_xyz() )
        self.xd,self.yd,self.zd = PHYS.rotate_ec_to_eq( self.get_uvw() )

class CartState(Structure):
    """
        ctypes "CartState" Structure for passing cartesian
        coords back and forth between python & C
        """
    _fields_ = [
                ('x', c_double),
                ('y', c_double),
                ('z', c_double),
                ('xd', c_double),
                ('yd', c_double),
                ('zd', c_double)
    ]
    def get_xyz(self):
        return self.x,self.y,self.z
    def get_uvw(self):
        return self.xd,self.yd,self.zd
    def ec_to_eq(self):
        self.x,self.y,self.z = PHYS.rotate_ec_to_eq(self.get_xyz() )
        self.xd,self.yd,self.zd = PHYS.rotate_ec_to_eq( self.get_uvw() )

class KepState(Structure):
    """
        ctypes "KepState" Structure for passing Keplerian
        coords back and forth between python & C
        """
    _fields_ = [
                ('a', c_double),
                ('e', c_double),
                ('incl', c_double),
                ('longnode', c_double),
                ('argperi', c_double),
                ('meananom', c_double)
    ]

class CartStateEpoch(Structure):
    """
        ctypes "CartStateEpoch" Structure for passing cartesian
        coords (at a definined epoch) back and forth between python & C
        
        *** UNUSED ??? ***
        
        """
    _fields_ = [
        ('CartState', CartState),
        ('epoch', c_double)
    ]
    def get_xyz(self):
        return self.state.x,self.state.y,self.state.z
    def get_uvw(self):
        return self.state.xd,self.state.yd,self.state.zd
    def ec_to_eq(self):
        self.state.x,self.state.y,self.state.z = PHYS.rotate_ec_to_eq(self.get_xyz())
        self.state.xd,self.state.yd,self.state.zd = PHYS.rotate_ec_to_eq(self.get_uvw())

class KepStateEpoch(Structure):
    """
        ctypes "KepStateEpoch" Structure for passing cartesian
        coords (at a definined epoch) back and forth between python & C
        
        *** UNUSED ??? ***
        
        """
    _fields_ = [
        ('KepState', KepState),
        ('epoch', c_double)
    ]
        

# (ii) Classes related to MPChecker
# --------------------------------------------------------------
class OrbitObject(Ele220):
    '''
        Generalized container to hold the information regarding the orbit (at epoch) of an object
        '''
    
    # Initializing
    GMcenter    = None
    epoch       = None
    cart_state  = None
    desig       = None
    Hmag        = None

    # Pretty-print
    def __str__(self):
        return "OrbitObject_%s" % self.desig
    
    # Load crucial quantities
    def load_cart_state(self, cart_state):
        self.cart_state=UK.CartState(cart_state.x, cart_state.y, cart_state.z, cart_state.xd, cart_state.yd, cart_state.zd)
    
    def load_GMcenter(self, GMcenter):
        self.GMcenter=GMcenter

    def load_desig(self, desig):
        self.desig=desig

    def load_epoch(self, epoch):
        self.epoch=epoch

    def load_Hmag(self, Hmag):
        self.Hmag=Hmag

    # Generate cart-state from ele220-keplerian-elements
    # - Used when starting from GW's nbody output (in ele220 format)
    def generate_cart_state_via_ele220(self, GMcenter, e220 ):
        q,e,i,node,AP,Tp = e220.periDist(), e220.ecc(), e220.inc(), e220.node(), e220.argPeri(), e220.timePeri()
        a = q/(1.-e)
        n = (180./np.pi)*np.sqrt(GMcenter / abs(a)**3)
        meananom = (e220.maEpoch() - Tp)*n

        self.load_GMcenter(GMcenter)
        self.load_epoch(e220.maEpoch())
        self.load_cart_state(kc.cartesian(GMcenter, a, e, np.radians(i), np.radians(node), np.radians(AP), np.radians(meananom)))
        self.load_desig(e220.desig())
        try:
            self.load_Hmag(e220.h())
        except:
            print("e220", e220)
            print(e220.h() )
    # Generate the orbit object by directly specifying all necessary quantities
    # - Used when loading directly from file
    def direct_specification(self, hmag, center, epoch, state, desig):
        self.load_GMcenter(center)
        self.load_epoch(epoch)
        self.cart_state=UK.CartState(state[0], state[1], state[2], state[3], state[4], state[5] )
        self.load_desig(desig)
        self.load_Hmag(hmag)
    

    # Generate RA, DEC. In equatorial coordinates.
    def get_RADEC(self, array_of_synthetic_observation_times_TT, array_of_observatory_positions_wrt_center_ECLIPTIC ):
        RA_Deg_arr,DEC_Deg_arr = UK.RADEC_from_state_C(self.cart_state,
                                                       self.GMcenter,
                                                       self.epoch,
                                                       array_of_synthetic_observation_times_TT,
                                                       array_of_observatory_positions_wrt_center_ECLIPTIC)
        return np.asarray(RA_Deg_arr),np.asarray(DEC_Deg_arr)
    
    # Generate unit-vector(s) from RA, DEC. In equatorial coordinates.
    def get_unit_vector(self, array_of_synthetic_observation_times_TT, array_of_observatory_positions_wrt_center_ECLIPTIC ):
        RA_Deg_arr,DEC_Deg_arr = self.get_RADEC(array_of_synthetic_observation_times_TT,
                                                array_of_observatory_positions_wrt_center_ECLIPTIC )
                                                
        return np.transpose(np.array([  np.cos(np.radians(RA_Deg_arr))*np.cos(np.radians(DEC_Deg_arr)),
                                        np.sin(np.radians(RA_Deg_arr))*np.cos(np.radians(DEC_Deg_arr)),
                                        np.sin(np.radians(DEC_Deg_arr))]))

    # Generate healpix from RA, DEC. In equatorial coordinates.
    def get_healpix(self, array_of_synthetic_observation_times_TT, array_of_observatory_positions_wrt_center_ECLIPTIC ):
        UV = self.get_unit_vector( array_of_synthetic_observation_times_TT, array_of_observatory_positions_wrt_center_ECLIPTIC )
        return np.array([hp.vec2pix(params.hp_nside, *uv, nest=params.hp_nested) for uv in UV ])


    # Generate "verbose" output to supplement the RA,Dec data returned from get_RADEC()
    def get_RADEC_VERBOSE(self, array_of_synthetic_observation_times_TT, array_of_observatory_positions_wrt_center_ECLIPTIC):
        '''Want to populate a "Prediction" with: name Ra Dec Ra_Rate Dec_Rate H vec_helio vec_topo Vmag'''
        
        # Get the RA & Dec
        RA_arr, Dec_arr = self.get_RADEC(array_of_synthetic_observation_times_TT, array_of_observatory_positions_wrt_center_ECLIPTIC )

        # Get the UNIT-vector to the object (but remember I ultimately need to return the VECTOR)
        uv = ANG.radec2unit( RA_arr[0], Dec_arr[0])
        
        # Return a "Prediction" named-tuple
        return [\
                self.desig,  #name
                RA_arr[0],   # RA
                Dec_arr[0],  # Dec
                random.random(), # RA_Rate
                random.random(), # Dec_Rate
                self.Hmag, # H
                np.array([random.random(),random.random(),random.random()]), # vec_helio
                np.array([random.random(),random.random(),random.random()]), # vec_topo
                0.  # Vmag
           ]



Prediction = collections.namedtuple('Prediction', 'name RA Dec RA_Rate Dec_Rate H vec_helio vec_topo Vmag')
'''
    # Namedtuple structure to store the data that will be returned from calling get_RADEC / get_VERBOSE methods on OrbitObject
'''




class Pointing():
    '''
        Generalized container to hold the information for each pointing
        '''
    
    def __init__(self, *args, name=None):#obsCode, JDutc, RA_center_DEG, Dec_center_DEG, FoV_

        # Initializing
        self.UV         = None
        self.HPcenter   = None # HP of the center of the pointing
        self.HPfov      = None # HP covered by the pointing
        self.name       = None

        assert len(args) in [4,5], 'Initialized w/ incorrect # args:can be 4(obsCode, JDutc, RA, Dec) or 5(obsCode, JDutc, RA, Dec, FoV)'

        # Have to initialize with at least obsCode, JDutc, RA, Dec (CENTER of pointing : Assumed to be in degrees)
        self.obsCode, self.JDutc = args[0], args[1]
        self.load_RADEC(args[2], args[3] )

        # Can optionally also initialize with FOV diameters in degrees
        if len(args) >= 5:
            self.load_FoV(args[4])
        else:
            self.FOV             = None # Degrees (~DIAMETER)
            self.FOV_rad         = None # Radians (~DIAMETER)
        
        if name != None:
            self.load_name(name)

    # Pretty-print
    def __str__(self):
        if self.name == None:
            return  "%s_%r_%r_%r" % ( self.obsCode , self.JDutc , self.RADEC[0], self.RADEC[1] )
        else:
            return self.name

    # Convenience functions for loading
    def load_RADEC(self, RAdeg, DECdeg):
        self.RADEC = (RAdeg,DECdeg)
    def load_UV(self, UV):
        self.UV=UV
    def load_FoV(self, FOV):
        self.FOV        = FOV
        self.FOV_rad    = np.radians(FOV)
    def load_name(self, name):
        self.name = name

    # Useful functions related to angles ...
    def calc_UV(self):
        ''' Calculate the UV of the center '''
        self.load_UV( ANG.radec2unit(self.RADEC[0] , self.RADEC[1]).flatten((-1,3)) )
    
    def calc_angles_wrt_center_radec(self, RA_, DEC_):
        ''' Calculate the angles between the center & and an array of input RAs & DECs : **DEGREES** '''
        return ANG.angle_radec(self.RADEC[0] , self.RADEC[1],  RA_, DEC_ ).flatten()
    
    def select_RADEC_within_FoV(self, RA_, DEC_):
        ''' Calculate which of the listed RA_,DEC_ fall within the FoV of the center : INDICEES '''
        return np.where( self.calc_angles_wrt_center_radec( RA_, DEC_) < self.FOV/2. )[0]
    

    # Useful functions related to Healpix ...
    def HP_of_center(self):
        '''
            Get the healpix of the center
            '''
        if np.any(self.UV == None) :
            self.calc_UV()
        self.HPcenter = hp.vec2pix(params.hp_nside, self.UV[0],self.UV[1],self.UV[2] , nest=params.hp_nested)
        return self.HPcenter

    def HP_covered_by_FoV(self, MarginForErrorRadians=0. , PrecalculatedHP=None):
        '''
            Get the list of healpix touched by the FoV 
            N.B. FOV_rad assumed to be DIAMETER or SIDE
            
            MarginForErrorRadians
                Allowing for a margin of error around the edge of the circulare FoV
            
            PrecalculatedHP:
                Allowing the option of just passing in the HP from an external calculation
            '''
        # If pass-in nothing, then assume circular FoV, and claculate HP
        if PrecalculatedHP == None:
            if np.any(self.UV == None) :
                self.calc_UV()
            assert not np.any(self.FOV_rad == None), 'No FoV specified: cannot calculate coverage'
            self.HPfov = hp.query_disc(params.hp_nside, self.UV , MarginForErrorRadians + self.FOV_rad/2., nest=params.hp_nested, inclusive=True )
        # Else, just set the HP using whetever's passed-in
        else:
            self.HPfov = PrecalculatedHP
        
        return self.HPfov

    # Useful functions related to position of the observatory at the time the pointing was taken ...
    # ... this is EQUATORIAL
    def calc_heliocentric_position_of_observatory_in_equatorial_coords(self):
        self.helio_eq_posn = MPCL.Observatory().getObservatoryPosition(self.obsCode, self.JDutc)
        return self.helio_eq_posn
    
    def calc_heliocentric_position_of_observatory_in_ecliptic_coords(self):
        self.helio_ec_posn = np.dot(np.transpose(PHYS.rot_mat), self.calc_heliocentric_position_of_observatory_in_equatorial_coords() )
        return self.helio_ec_posn



# (iii) Classes related to kepler-step (2body orbit advance)
# --------------------------------------------------------------

class CartPartial(Structure):
    """
    ctypes "CartPartial" Structure for passing cartesian
    partial derivatives back and forth between python & C
    Upper-case letters are the components of the final cartesian state
    Lower-case letters are the components of the input cartesian state
    """
    _fields_ = [
                ('dXdx', c_double),
                ('dXdy', c_double),
                ('dXdz', c_double),
                ('dXdxd', c_double),
                ('dXdyd', c_double),
                ('dXdzd', c_double),
                ('dYdx', c_double),
                ('dYdy', c_double),
                ('dYdz', c_double),
                ('dYdxd', c_double),
                ('dYdyd', c_double),
                ('dYdzd', c_double),
                ('dZdx', c_double),
                ('dZdy', c_double),
                ('dZdz', c_double),
                ('dZdxd', c_double),
                ('dZdyd', c_double),
                ('dZdzd', c_double),
                ('dDXdx', c_double),
                ('dDXdy', c_double),
                ('dDXdz', c_double),
                ('dDXdxd', c_double),
                ('dDXdyd', c_double),
                ('dDXdzd', c_double),
                ('dDYdx', c_double),
                ('dDYdy', c_double),
                ('dDYdz', c_double),
                ('dDYdxd', c_double),
                ('dDYdyd', c_double),
                ('dDYdzd', c_double),
                ('dDZdx', c_double),
                ('dDZdy', c_double),
                ('dDZdz', c_double),
                ('dDZdxd', c_double),
                ('dDZdyd', c_double),
                ('dDZdzd', c_double)
                ]






"""

# (ii) Classes related to detections
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class UV(object):
    def __init__(self, RArad, DECrad):
        '''
        # Unit Vector from (RA,Dec) [input in radians]
        '''
        self.RA, self.DEC, self.theta  = RArad, DECrad , 0.5*np.pi-DECrad 
        self.UV                        = hp.ang2vec(self.theta , self.RA) 

class UVHP(UV):
    def __init__(self, RArad, DECrad):
        '''
        # Unit Vector & Healpix from (RA,Dec) [input in radians]
        '''
        # Make RA/Dec/Theta/UV available from UV
        #super().__init__(RArad, DECrad) 
        UV.__init__(self, RArad, DECrad) 
        # Make HP
        self.HP = hp.ang2pix(HPPar().nside, self.theta, self.RA, HPPar().NESTED)
        

    def get_UVHP(self):
        return self.UV, self.HP


class Pointing(UVHP):
    ''' 
    # This is intended to be a single pointing as supplied by the surveys
    # 
    # - At present it's not much more than elaborate UnitVector
    # - What's the point?
    # - Should it be broadened out to contain more info (connect to exposures, etc etc)?
    #
    # - This was taken from the Pointer stuff I wrote for Sonia
    '''

    def __init__(self, RArad, DECrad, JD_UTC, OBS_CODE, FOV_rad):      
        # Direct Inputs
        self.RA       = RArad
        self.DEC      = DECrad
        self.JD_UTC   = JD_UTC
        self.OBS_CODE = OBS_CODE
        self.FOV_rad  = FOV_rad

        # Exposure this pointing may be associated with  
        self.exp_str       = None                                   
    
        # It is likely to be useful to know where the observatory was at the time of the pointing 
        self.Helio_Posn = Space.pipeline_OC_Posn( self.OBS_CODE, self.JD_UTC )
        
        # Make RA/Dec/Theta/UV/HP available from UV
        UVHP.__init__(self, RArad, DECrad) 
                
    def __str__(self):
        return "%s_%s_%s_%s" % (self.RA, self.DEC, self.JD_UTC, self.OBS_CODE) 
        
    def set_FOV_rad(self, FOV_rad):
        self.FOV_rad = FOV_rad

    def add_pointing_association(self,exp):
        isinstance( exp , Exposure) 
        self.exp_str = str(exp)

class DetectionData(object):
    ''' For now I am going to make this able to hold an obs_80 namedtuple. 
        Would hope to expand/generalize to be able to hold data from anything in the future
    '''

    def __init__(self, single_valid_detection):
        if isinstance( single_valid_detection , obs80.Optical ):
            self.extract_obs80_Optical(single_valid_detection)

    def extract_obs80_Optical(self, obs_80):
        # Extract the rest of the quantities / make obs80 consistent
        self.RA       = obs_80.ra  * PHYS.hr2deg *np.pi/180. # because obs80.ra in [hr] not [deg]
        self.DEC      = obs_80.dec * np.pi/180. 
        self.JD_UTC   = obs_80.jdutc
        self.OBS_CODE = obs_80.cod



class Detection(DetectionData, UVHP):

    def __init__(self, single_valid_detection):

        # Make detection data available from parent
        DetectionData.__init__(self, single_valid_detection)

        # Make RA/Dec/Theta/UV/HP available from parent
        UVHP.__init__(self, self.RA,self.DEC) 

        # It is likely to be useful to know where the observatory was at the time of the detection 
        # Presumably equatorial at present
        self.Helio_Posn = Space.pipeline_OC_Posn( self.OBS_CODE, self.JD_UTC )

        # Other quantities to be dealt with later 
        self.Submission = None     # Submission (Batch) the detection is in
        self.Tracklet   = None     # Tracklet the detection is in
        self.Exposure   = None     # Exposure the detection is in
        self.OrbObj     = None 	   # Orbit-Object with which this detection is associated 

    def __str__(self):      # Unique "name" for detection 
        return "_".join((str(item) for item in (self.JD_UTC, self.RA, self.DEC, self.OBS_CODE)))
   
    def Set_Submission(self, sub):
        '''Set the Submission (batch) this detection is part of '''
        self.Submission = sub

    def Set_Tracklet(self,tracklet):
        '''Set the tracklet this detection is part of '''
        self.Tracklet = tracklet

    def Set_Exp(self, exp):
        '''Set the exposure this detection is part of '''
        self.Exp = exp

    def Set_OrbObj(self , Obj):
        '''Set the Orbit-Object this detection is part of '''
        self.OrbObj = Obj
        
 

class Tracklet(list):
    ''' I think I want to make a Tracklet intrinsically be a list object (same as VO_Bundle below)
        I think I want to broaden out the definition of a Tracklet to allow for mixed observatory-codes (thinking of TNOs)
    '''

    def __init__(self, Desig=None, Number=None):
      '''Initialize a new tracklet '''
      
      self.Desig      = Desig	                                           # observatory track designation
      self.Number     = Number                                             # enumerating after splitting track into tracklets  

    def __str__(self):
          return "%s_%s_%s" % (self.ObsCode , self.Desig , self.Number)    # Tracklet name

    def keys(self):
        return [str(det) for det in self]

    def OCs(self):
        return [det.OBS_CODE for det in self]

    def helios(self):
        return [det.Helio_Posn for det in self]

    def JDs(self):
        return [det.JD_UTC for det in self]

    def add_det(self, det):
        ''' 
        # Add Detection to Tracklet
        # Allow append if (i) currently empty, or (ii) detection is different to current entries 
        '''
        assert isinstance( det, Detection) , "Not a detection"
        if len(self) == 0 or (str(det) not in self.keys() ) :
            self.append(det)
     
    def add_dets(self, dets):
        assert isinstance( dets, (list,tuple) ) , "Not iterable"
        for det in dets: 
            self.add_det(det)
    
    

class Exposure(list):
    ''' 
    # Collect detections from the same time & obs-code
    # - Probably want to make it a little more similar to the tracklet object (fewer associated lists, just get data from the contained detections using a class-function call)
    '''

    def __init__(self, OBS_CODE=None, JD_UTC=None):
      '''Initialize a new exposure '''
      
      self.OBS_CODE      = None	                                  # observatory code 
      self.JD_UTC        = None                                   # jd-utc
      self.pnt_str       = None                                   # Pointing this exposure may be associated with  

      self.keys          = []                                     #  array of det key-strings
      self.HPs           = set()                                  #  *set* of det HPs
      self.HPs_extended  = set()                                  #  *set* of extended-HPs (derived from HP)

      self.det_match_ind = []                                     #  array to keep track of which dets have been matched to a known object 

    def __str__(self):
          return "%s_%s" % (self.ObsCode , self.DateTime)         # Exposure name
        
    def add_det(self, det):
      '''Add Detection to Exposure'''
      isinstance( det, Detection) 
      # Enforce that det.ObsCode & det.DateTime are the same as the other obs already in this Exposure
      if len(self) == 0 or ( det.OBS_CODE == self.OBS_CODE and det.JD_UTC == self.JD_UTC) :
          if len(self) == 0:
              self.OBS_CODE = det.OBS_CODE
              self.JD_UTC   = det.JD_UTC

          self.append(det)
          self.keys.append(str(det))
          self.HPs.add(det.HP)
          self.HPs_extended.add( det.HP )
          for item in hp.get_all_neighbours(HPPar().nside,  det.HP, nest=HPPar().NESTED):
              self.HPs_extended.add( item )

    def add_dets(self, dets):
        for det in dets: self.add_det(self, det)

    def add_det_match_indicees(self,ind):
        self.det_match_ind.append(ind)

    def add_pointing_association(self,pnt):
        isinstance( pnt, Pointing) 
        self.pnt_str = str(pnt)





# (iii) Classes related to object-orbits 
# --------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------

class VarOrb(object):
    ''' # This is intended to be a single variant orbit associated with a single NEO
        # - This was taken from the Pointer stuff I wrote for Sonia
        # Might want to generalize the handling of the orbit string (like for obs80) 
    '''

    def __init__(self, orb_str , ANG_TYPE="DEG"): # jd_utc , 
        
        # Split string
        self.Object,   self.H, self.G, self.Epoch, self.M, self.Peri, self.Node, self.Incl, self.e, self.n, self.a, self.NObs, self.NOpp, self.Arc, self.ArcUnit, self.rms, self.OrbitID = orb_str.split()

        # Make floats
        self.H     = float(self.H)
        self.G     = float(self.G)
        self.M     = float(self.M)
        self.Peri  = float(self.Peri)
        self.Node  = float(self.Node)
        self.Incl  = float(self.Incl)
        self.e     = float(self.e)
        self.n     = float(self.n)
        self.a     = float(self.a)
        self.NObs  = int(self.NObs)
        self.NOpp  = int(self.NOpp)
        self.Arc   = float(self.Arc)
        self.rms   = float(self.rms)

        # Ensure radians
        if ANG_TYPE == "DEG":
            self.M     *=  np.pi/180. 
            self.Peri  *=  np.pi/180. 
            self.Node  *=  np.pi/180. 
            self.Incl  *=  np.pi/180. 
            self.n     *=  np.pi/180.  # N.B. At this point the mean motion will be in radians per **DAY**
             
        # Epoch to jd
        self.jd_Epoch = Time.JDFromPacked(self.Epoch)

        '''
        # Julian Date array
        self.dates_ = CHEBY.cheby_times(jd_utc )

        # Ecliptic Cartesian coords at desired dates ::: cartesian(GM, a, e, incl, longnode, argperi, meananom)
        self.meananoms_ = np.array( [ self.M + self.n * ( date - self.jd_Epoch) for date in self.dates_ ] )
        self.carts_     = np.array( [ kc.cartesian(PAR.GMsun, self.a, self.e, self.Incl, self.Node, self.Peri, MA) for MA in self.meananoms_ ] )
        self.posns_     = np.array( [ [cart['x'],cart['y'],cart['z'] ] for cart in self.carts_ ] )
        
        # Transform to equatorial coordinates to facilitate comparison with equatorial RA/Dec observations 
        self.eq         = PAR.rot_mat.dot( self.posns_.T ).T
               
        # Generate Chebyshev coeffs for object (equatorial coords)
        self.cheby      = CHEBY.cheby_create(self.dates_ , self.eq )
                
        # Convenience Calc of LTT step (from the geocenter)
        self.geo_LTT    = CHEBY.cheby_LTT(self.dates_ , self.cheby)
        
         # Chebys for LTT  
        self.cheby_LTT  = CHEBY.cheby_create(self.dates_ , self.geo_LTT)
        '''
 
    def __str__(self):
          return "_".join((self.Object, self.OrbitID)) 


class VO_Bundle(list):
    ''' # This is intended to be a "bundle" of variant orbits associated with a single NEO
        # - This was taken from the Pointer stuff I wrote for Sonia
        # - Has been altered to inherit list properties
    '''

    def __init__(self, data=None ):
        '''Initialize a new exposure '''
      
        self.Object           = None
        self._orbit_keys      = set()                              # Convenient *set* of obs key-strings
        self._obs_match_ind   = []                                 # Convient array to keep track of which _obs have been matched to a known object 
        if data is not None:
            if   isinstance( data , VarOrb) : self.add_orbit(data)  
            elif isinstance( data , list )  : self.add_orbits(data)            
            else:                             sys.exit("Don't know how to handle")

    def __str__(self):
        return "%s" % (self.Object)                                # Bundle name == Constituent object name 
         
    def add_orbit(self, orbit):
        ''' # Add a variant orbit to this bundle 
        ''' 
        # If strings, make into var-orbit-objects
        if not isinstance( orbit , VarOrb ) : 
            orbit = VarOrb(orbit)
        # Allow append if (i) currently empty, or (ii) var-orbit is different to current entries while having the same object-id
        if len(self) == 0 or (str(orbit) not in self._orbit_keys and orbit.Object == self[0].Object) :
            self.append(orbit)
            self._orbit_keys.add(str(orbit))
        # Set the name
        if len(self) == 1:
            self.Object = self[0].Object

    def add_orbits(self, orbits):
        ''' # Add a list of variant orbits to the bundle 
            # - Would probably be more efficient to write as separate method ... 
        ''' 
        for orbit in orbits:
            self.add_orbit(orbit)



        
'''
class NBODY_ORB(object):

    def __init__(self, JD_ref , dt_TDB , dt_min , dt_max , n_coeff, orbit_spec):
        # JD_ref , dt_TDB , dt_min , dt_max , n_coeff , desig , RA_deg , Dec_deg , ltt_day , x_coeff , y_coeff , z_coeff , cart , elem , hp_ 
        # desig RA_deg Dec_deg ltt_day xcoeffs... ycoeffs... zcoeffs... x y z xd yd zd a e incl_deg longnode_deg argperi_deg meananom_deg hp_count hp0 hp1...
        desig , RA_deg , Dec_deg , ltt_day = arr_[0:4] 
        x_coeff = arr_[4:4+n_coeff]
        y_coeff = arr_[4+n_coeff   : 4+2*n_coeff      ]
        z_coeff = arr_[4+2*n_coeff : 4+3*n_coeff      ]
        cart    = arr_[4+3*n_coeff : 4+3*n_coeff + 6  ]
        elem    = arr_[3*n_coeff + 10 : 3*n_coeff + 16]
        n_hp    = arr_[3*n_coeff + 16 ]
        hp_     = [hp for hp in arr_[3*n_coeff + 17 : ] ]
'''



#For now I am going to make this able to hold 
#Will probably expand/generalize in the future 
'''
class OrbitData(object):

    def __init__(self, single_valid_orbit):
        if isinstance(   single_valid_orbit , NBODY_ORB ):
            self.extract_NBODY_ORB(single_valid_orbit)
        elif isinstance( single_valid_orbit , ele220.ele220 ): 
            self.extract_MPC_ORB(single_valid_orbit)
        else:
            sys.exit("OrbitData class doesn't recognize the type of orbit supplied")

    def extract_NBODY_ORB(self, orbit_spec):
        # For orbit in the form output by Holman/Payne integration scheme
        # JD_ref , dt_TDB , dt_min , dt_max , n_coeff , desig , RA_deg , Dec_deg , ltt_day , x_coeff , y_coeff , z_coeff , cart , elem , hp_ = ...
        pass

    def extract_MPC_ORB(self, orbit_spec):
        # For orbit in the MPC220 format, using SOnia's reader 
        pass


    # Need to define / link methods which advance orbits
    # 
    # But depends what's done / what's available
    #
    # The MPC-Orbits will need to be able to be advanced using (a) Keplerian, and perhaps (b) N-Body
    #
    # The NBODY-Orbits have been evolved *from* MPC-Orbits, and so will be in Chebyshev form
    # - But they will only be available / valid over a limited number of time-spans. 




'''  #An object-orbit object...
'''
class OrbObject(OrbitData):

    def __init__(self, data=None):
        self.desig = None 

        # Make orbit specification data available from parent
        if data is not None:
            OrbitData.__init__(self, data)
            
    def __str__(self):      # Unique "name" for observation 
        return self.desig


        

'''# Initialize a new object ...
'''
class Object(object):

    def __init__(self, JD_ref , dt_TDB , dt_min , dt_max , n_coeff , desig , RA_deg , Dec_deg , ltt_day , x_coeff , y_coeff , z_coeff , cart , elem , hp_ ):
        self.JD_ref      = float(JD_ref)
        self.dt_TDB      = float(dt_TDB)
        self.dt_min      = float(dt_min)
        self.dt_max      = float(dt_max)

        self.n_coeff     = int(n_coeff)
 
        self.desig       = str(desig)
        self.RA_deg      = float(RA_deg)
        self.Dec_deg     = float(Dec_deg)
        self.ltt_day     = float(ltt_day)
        self.x_coeff     = [float(item) for item in x_coeff]
        self.y_coeff     = [float(item) for item in y_coeff]
        self.z_coeff     = [float(item) for item in z_coeff]
        self.cart        = [float(item) for item in cart]
        self.elem        = [float(item) for item in elem]
        self.hp_         = [int(item)   for item in hp_]
        
        self.UnitVec     = {}


    def __str__(self):      # Unique "name" for observation 
        return self.desig

    def add_UnitVec(self, obs_code, exp_time , UnitVec):
        #Associate a unit vector with this object
        #I expect this to be used to hold an observer-centric UV. 
        #Am assuming the time is the exposure time in JD_UTC
      self.UnitVec[(obs_code, exp_time)] = UnitVec
      
    def pretty_print(self):
        print("Obj:       ", self.desig       )
        
        print("JD:        ",self.JD_ref      )  
        print("dt TDB:    ",self.dt_TDB      ) 
        print("dt-        ",self.dt_min      )
        print("dt+        ",self.dt_max      )

        print("Cheb-n:    ",self.n_coeff     )

        print("RA:        ",self.RA_deg      )
        print("Dec:       ",self.Dec_deg     )
        print("LTT(day):  ",self.ltt_day     )
        print("Cheb-x:    ",self.x_coeff     )
        print("Cheb-y:    ",self.y_coeff     )
        print("Cheb-z:    ",self.z_coeff     )
        print("Cartesians:",self.cart        )
        print("Elements:  ",self.elem        )
        print("HP:        ",self.hp_         )
        
        return True
'''        


"""


