"""
# 
# --------------------------------------------------------------
# time.py
# July 2016 
# Payne & Holman

# Functions related to 
# (i) Physical Constants
# (ii) Unit-Conversions

# If efficient, might want to replace some of this with the JPL/SPICE stuff
# Or novas
# Or astropy
# Whatever

# --------------------------------------------------------------

"""

import numpy as np

# --------------------------------------------------------------
GMsun          = 2.9591220828559115e-04
# Obliquity of ecliptic at J2000
#ecl            = (84381.4118*(1./3600)*np.pi/180.)
ecl            = (84381.448*(1./3600)*np.pi/180.)
# This is now a definition
au_km          = 149597870.700                 
speed_of_light = 2.99792458e5 * 86400./au_km
# 1e-2 seconds -> 1e-7 days
LLT_Tol        = 1e-2 / (3600*24)

def rotate_matrix(ecl):
    ce = np.cos(ecl)
    se = np.sin(-ecl)
    rotmat = np.array([[1.0, 0.0, 0.0],
                       [0.0,  ce,  se],
                       [0.0, -se,  ce]])
    return rotmat
rot_mat        = rotate_matrix(ecl)
Rearth_km      = 6378.1363
Rearth_AU      = Rearth_km/au_km
hr2deg         = 15 
