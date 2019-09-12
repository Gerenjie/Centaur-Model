"""
# --------------------------------------------------------------
# Jul 2017 
# Payne & Holman

# Use C to do fast ...
# ... (i) Coordinate conversions (Keplerian <-> Cartesian)
# ...(ii) Orbit advance (kepler-stepper)
 
# --------------------------------------------------------------
"""

import os
import numpy as np
from ctypes import *

# Importing of local modules/packages
# --------------------------------------------------------------
import classes as Classes

# Import local files / dirs  
# --------------------------------------------------------------
lib = CDLL(os.path.join(os.path.dirname(__file__), 'libkepcart.so'))



# --------------------------------------------------------------
def keplerian(GM, state):
    """
    Computes the Keplerian orbital elements a, e, incl, longnode,
    argperi, and meananom, given a GM constant and an input state.
    
    *Returns*
        (a, e, incl, longnode, argperi, meananom) : tuple of floats
    
    """

    _keplerian = lib.keplerian
    _keplerian.argtypes = (c_double, Classes.State)
    _keplerian.restype = None

    a = c_double()
    e = c_double()
    incl = c_double()
    longnode = c_double()
    argperi = c_double()
    meananom = c_double()

    return_value = _keplerian(GM, state, byref(a), byref(e), byref(incl), byref(longnode), byref(argperi), byref(meananom))

    return (a.value, e.value, incl.value, longnode.value, argperi.value, meananom.value)

def keplerians(num, GM, state_arr):
    """
    Computes arrays of Keplerian orbital elements a, e, incl, longnode,
    argperi, and meananom, given a GM constant and an array of input states.
    
    *Returns*
    numpy arrays of a, e, incl, longnode, argperi, meananom
    
    """

    StateArray = Classes.State * num

    a_arr = np.zeros((num), dtype=np.double)
    e_arr = np.zeros((num), dtype=np.double)
    incl_arr = np.zeros((num), dtype=np.double)
    longnode_arr = np.zeros((num), dtype=np.double)
    argperi_arr = np.zeros((num), dtype=np.double)
    meananom_arr =np.zeros((num), dtype=np.double)

    _keplerians = lib.keplerians
    _keplerians.argtypes = (c_int, c_double, POINTER(StateArray), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double))
    _keplerians.restype = None


    return_value = _keplerians(num, GM, byref(state_arr),
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)))

    return a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr




def cartesian(GM, a, e, incl, longnode, argperi, meananom):
    """
    Computes the cartesian state given a GM constant and the orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
        (x, y, z, xd, yd, zd) : tuple of floats
    
    """
    
    _cartesian = lib.cartesian
    _cartesian.argtypes = (c_double, c_double, c_double, c_double, c_double, c_double, c_double, POINTER(Classes.State))
    _cartesian.restype = None

    state = Classes.State(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    return_value = _cartesian(GM, a, e, incl, longnode, argperi, meananom, byref(state))

    return state

def cartesians(num, GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Computes the cartesian states given a GM constant and arrays of the orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
    numpy arrays of x, y, z, xd, yd, zd states. : tuple of floats
    
    """

    StateArray = Classes.State * num
    state_arr = StateArray()

    _cartesians = lib.cartesians
    _cartesians.argtypes = (c_int, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(StateArray))
    _cartesians.restype = None


    return_value = _cartesians(num, GM,
                               a_arr.ctypes.data_as(POINTER(c_double)),
                               e_arr.ctypes.data_as(POINTER(c_double)),
                               incl_arr.ctypes.data_as(POINTER(c_double)),
                               longnode_arr.ctypes.data_as(POINTER(c_double)),
                               argperi_arr.ctypes.data_as(POINTER(c_double)),
                               meananom_arr.ctypes.data_as(POINTER(c_double)),
                               byref(state_arr))

    return state_arr


def cartesian_vectors(num, GM, a_arr, e_arr, incl_arr, longnode_arr, argperi_arr, meananom_arr):
    """
    Computes the cartesian position and velocity vectors given a GM constant and arrays of the orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
    arrays of x, y, z and xd, yd, zd values.
    
    """

    size = num*3
    array_of_size_doubles = c_double*size

    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()

    _cartesian_vectors = lib.cartesian_vectors
    _cartesian_vectors.argtypes = (c_int, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _cartesian_vectors.restype = None


    return_value = _cartesian_vectors(num, GM,
                                      a_arr.ctypes.data_as(POINTER(c_double)),
                                      e_arr.ctypes.data_as(POINTER(c_double)),
                                      incl_arr.ctypes.data_as(POINTER(c_double)),
                                      longnode_arr.ctypes.data_as(POINTER(c_double)),
                                      argperi_arr.ctypes.data_as(POINTER(c_double)),
                                      meananom_arr.ctypes.data_as(POINTER(c_double)),
                                      byref(pos_arr),
                                      byref(vel_arr))

    return pos_arr, vel_arr


def cartesian_elements(num, GM, elements_arr):
    """
    Computes the cartesian position and velocity vectors given a GM constant and an array of sets of orbital elemments 
    a, e, incl, longnode, argperi, and meananom.
    
    *Returns*
    arrays of x, y, z and xd, yd, zd values.
    
    """

    ElementsArray = Classes.Elements * num
    
    size = num*3
    array_of_size_doubles = c_double*size

    pos_arr = array_of_size_doubles()
    vel_arr = array_of_size_doubles()

    _cartesian_elements = lib.cartesian_elements
    _cartesian_elements.argtypes = (c_int, c_double, POINTER(ElementsArray), POINTER(array_of_size_doubles), POINTER(array_of_size_doubles))
    _cartesian_elements.restype = None

    return_value = _cartesian_elements(num, GM, byref(elements_arr),
                                      byref(pos_arr),
                                      byref(vel_arr))

    return pos_arr, vel_arr


