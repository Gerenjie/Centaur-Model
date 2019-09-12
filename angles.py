"""

--------------------------------------------------------------

Mar 2019

Payne

Useful functions related to angles ...

--------------------------------------------------------------

"""
import healpy as hp

# Import third-party packages
# --------------------------------------------------------------
import numpy as np

# Import third-party packages
# --------------------------------------------------------------
def radec2unit( RAdeg, DECdeg):
    '''
        Function to convert RA, DEC to Unit-vector
        Assumed to be input in degrees
        '''
    RA_rad  = np.radians(RAdeg)
    DEC_rad = np.radians(DECdeg)
    return np.transpose(np.array([np.cos(RA_rad)*np.cos(DEC_rad), np.sin(RA_rad)*np.cos(DEC_rad),np.sin(DEC_rad)]))

def unit2radec(UV):
    '''
        Function to convert unit-vector to RA, DEC.
        Output in degrees
        '''
    
    # Make it at least 2d, even if a 1d (3-comp) array is supplied
    
    UV = np.atleast_2d(UV)
    assert UV.shape[1] == 3 , 'UV is not of the correct shape : %r' % UV.shape
    
    
    norm    = np.linalg.norm(UV, axis=1)
    DEC     = 0.5*np.pi - np.arccos(UV[:,2]/norm)
    
    # Initialize to zero
    RA      = np.zeros(len(UV))
    # Take care regarding corner case of x==y==0
    ind     = np.where( np.logical_or(UV[:,0]!=0, UV[:,1]!=0) )
    RA[ind] = np.arctan2(UV[:,1][ind],UV[:,0][ind])   # in ]-pi, pi]
    # Re-set range
    ind     = np.where(RA<0.)
    RA[ind] = RA[ind]+2.*np.pi                        # in  [0, 2pi[
    
    return np.degrees(RA) , np.degrees(DEC)

def angle_unitvectors(uv1,uv2):
    '''
        Calculate the angle between two unit vectors
        Returns angle in radians
        '''
    # Check the format is as required
    assert (uv1.shape == (3,) or uv1.shape == (1,3)), 'The behavior of this routine has only been tested when uv1 is a single unit vector'
    assert (uv2.shape == (3,) or uv2.shape[1] ==3 )
    
    # Do the dot-products to get the angle
    uv1.flatten()
    uv2 = np.atleast_2d(uv2)
    dot = np.dot( uv1,uv2.T )
    
    # Return results in radians
    return np.arccos( dot )

def angle_radec(RA1,Dec1,RA2,Dec2):
    '''
        Calculate the angle between two sets of RA & Dec
        Returns angle in degrees
        '''
    return np.degrees( angle_unitvectors( radec2unit(RA1,Dec1), radec2unit(RA2,Dec2)) )
