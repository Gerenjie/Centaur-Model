import healpy as hp
import numpy as np

# -----------------------------------
# ---- Healpy/Healpix Parameters ----
# -----------------------------------
hp_nested = True
hp_nside  = 256
hp_npix   = hp.nside2npix(hp_nside)
hp_area   = 4*np.pi*(180./np.pi)**2 / hp_npix
hp_side   = np.sqrt(hp_area)
