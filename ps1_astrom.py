# This is a code from Eddie Schlafly
# Matt Holman modified it to work with astropy

from astropy.io import fits
from scipy import optimize
import numpy as np
import os, fnmatch, re, pdb
import ps1_utils
from collections import defaultdict

def getcd(smffile, chip):
    h = fits.getheader(smffile, extname=chip+'.HDR')
    cd = np.array([[float(h['CD1_1A']), float(h['CD1_2A'])],
                      [float(h['CD2_1A']), float(h['CD2_2A'])]])
    return cd

def get_ps1_coords(h):
    coord = { }
    grab_coords = ['CDELT1', 'CDELT2', 'CRPIX1', 'CRPIX2',
                   'CRVAL1', 'CRVAL2',
                   'PC001001', 'PC001002', 'PC002001', 'PC002002',
                   'NPLYTERM']
    translate_dict = {'ctype1':'ctype',
                      'pc001001':'pc1_1', 'pc001002':'pc1_2',
                      'pc002001':'pc2_1', 'pc002002':'pc2_2',
                      'nplyterm':'npolyterms'}

    for s in grab_coords:
        translation = translate_dict.get(s.lower(), s.lower())
        try:
            coord[translation] = float(h[s])
        except KeyError as e:
            if translation == 'npolyterms':
                print ('Warning: nplyterm not set for chip %s, setting 0?'
                       % h['EXTDATA'][0:4])
                coord[translation] = 0
            else:
                coord[translation] = np.nan
    if np.isnan(coord['cdelt1']):
        return coord
    coord['ctype'] = h['CTYPE1']
    coord['polyterms'] = np.zeros(14)

    polytermnames = ['PCA1X2Y0', 'PCA2X2Y0', 'PCA1X1Y1', 'PCA2X1Y1',
                     'PCA1X0Y2', 'PCA2X0Y2', 'PCA1X3Y0', 'PCA2X3Y0',
                     'PCA1X2Y1', 'PCA2X2Y1', 'PCA1X1Y2', 'PCA2X1Y2',
                     'PCA1X0Y3', 'PCA2X0Y3']
    if coord['npolyterms'] > 1:
        for i in range(6):
            coord['polyterms'][i] = float(h[polytermnames[i]])
    if coord['npolyterms'] > 2:
        for i in range(8):
            coord['polyterms'][i+6] = float(h[polytermnames[i+6]])
    return coord

def rd_to_xy(r, d, image, mosaic, xguess=None, yguess=None):
    r = np.atleast_1d(r)
    d = np.atleast_1d(d)
    def chi2(x, fitrd):
        rp, dp = xy_to_rd(x[0], x[1], image, mosaic)
        return np.array([rp-fitrd[0], dp-fitrd[1]])
    
    outx = np.zeros_like(r)
    outy = np.zeros_like(r)
    for i in range(len(r)):
        if xguess==None and yguess==None:
            guess = np.array([mosaic['crval1'], mosaic['crval2']])
        else:
            guess = np.array([xguess[ig], yguess[ig]])
        fit = optimize.leastsq(chi2, guess, args=[r[i], d[i]])
        outx[i] = fit[0][0]
        outy[i] = fit[0][1]
    if len(r) == 1:
        outx = outx[0]
        outy = outy[0]
    return outx, outy

# MJH 30 March 2019
# 
# I thought I would add a little documentation, because I didn't
# understand how xy_to_lm, lm_to_rd were working together in xy_to_rd.
# In particular, I didn't understand why xy_to_rd was calling
# itself.  The resulting documentation is longer than I anticpated!
#
# Note: xy_to_rd is not generally called externally by a user. It is
# called by xy_to_rd_smfcoords, which supplies the correct WCS dictionaries.
#
# 1. xy_to_rd gets called with chip-level x,y values from
#    a specific detector.  The units are pixels.  The WCS for that
#    detector is in the 'image' dictionary.  Its 'ctype' ends with
#   'WRP', which stands for 'warp'.  The 'mosaic' dictionary contains
#    the WCS for the mosaic as a whole.
# 
# 2. The first call to xy_to_lm takes that those x,y values,
#    along with the WCS for the chip, in the 'image' dictionary.
#    xy_to_lm takes the chip-level x,y coordinates and applies a stretch,
#    a rotation, and various higher-order polynomial terms to
#    generate x,y coordinates such that those for each chip are
#    parallel and have the same scale.  However, the x,y values that
#    are returned have not yet been shifted to place them properly in
#    the overall mosaic.  The units are still pixels.
#
# 3. Next there is a call to lm_to_rd.  At this stage the 'ctype' is
#    that of the 'image' dictionary, i.e. 'WRP' rather than 'DIS'.  So
#    the 'not dis' clause of lm_to_rd is entered.  At this stage and in this
#    branch, the input l,m values are actually just the stretched,
#    rotated, and corrected x,y values from the call to xy_to_lm in the
#    previous step.   Those l,m values just get translated by the 'crval1'
#    and 'crval2' values.  These place the new chip coordinates in their
#    proper place within the mosaic.  The results are assigned to ra,dec
#    but really they are still x,y coordinates in the overall mosaic.
#    The pair called ra,dec (which again is still x,y) is returned at the
#    end of lm_to_rd.  The units are still pixels.  
#    
# 4. The results of the call to lm_to_rd are returned to xy_to_rd and stored
#    in r,d.
#
# 5. The 'ctype' field is still 'WRP', so warp is True and the 'if warp'
#    branch is entered.  At this stage xy_to_rd is called AGAIN, but this time
#    it is called with the r,d values from the previous step.  Those are now
#    global x,y values for the whole mosaic, i.e. 'warp' values.  And
#    xy_to_rd is now called with the 'mosaic' WCS, rather than the 
#    'image' (chip) WCS.
#
#    My guess is that the overall transformation that results from steps 1-5
#    is essentially the same from exposure to exposure.  The focal plane might
#    expand and contract with temperature, and it might flex as a function of
#    the Alt-Az and rotator positions.  But those terms are fairly small.  I
#    will check a bunch of smf files to see.  Of course, if a chip is replaced
#    or reoriented, that will be evident in those WCS terms.
#
# 6. We are now back at the top of xy_to_rd.  And xy_to_lm is called again, this
#    time with the global x,y values and the 'mosaic' WCS (which is now stored
#    in 'image').  As before, xy_to_lm does a stretch, rotation, and higher-order
#    corrections on what are now the 'warp' x,y values.  The stretch converts the
#    units from pixels to degrees.  The rotation aligns the axes with North up and
#    East left.  The higher-order terms, now applied in degrees, to  account for
#    other effects.  The results from this call to xy_to_lm are returned
#    to xy_to_rd and are stored in l,m.
# 
# 7. Next, lm_to_rd is also called again with the 'mosaic' WCS.  Within lm_to_rd
#    the 'ctype' value is now 'DIS (what is it in the 'PRIMARY' extension for the
#    overall mosaic.  So the 'else' branch is entered.  Within that branch a standard
#    transformation from l,m (tangent plane) to RA,Dec is performed.  The center of
#    the tangent plane projection is in 'CRVAL1' and 'CRVAL2'.  The results of the
#    transformation are stored in ra,dec.  The units are still degrees, and the ra
#    value is mod'ed to [0, 360).  At the end of lm_to_rd the ra,dec pair is returned
#    to xy_to_rd.  It is stored in r,d.
#
# 8. At this stage, the r,d pair is returned either to xy_to_rd_smfcoords or to
#    any other external call to xy_to_rd.  And that's the end.  Phew!
#
#
#    Note: I am still unsure about what the Pan-STARRS folks call a 'warp').
#    From the context of stacking images and sky cells, I believe that a warp is
#    the result of doing a local tangent plane projection with an adopted plate
#    scale.  
#
#   
def xy_to_rd(x, y, image, mosaic, debug=False):
    #x = x.astype('f8')
    #y = y.astype('f8')
    if np.isnan(image['crval1']) or np.isnan(mosaic['crval1']):
        raise ValueError('Image or mosaic has no astrometric solution; probably missing chip')
    l, m = xy_to_lm(x, y, image)
    r, d = lm_to_rd(l, m, image)
    warp = (image['ctype'][-3:] == 'WRP')
    if warp:
        dis = (mosaic['ctype'][-3:] == 'DIS')
        if not dis:
            raise ValueError('Must set mosaic with chip astrometry.')
        r, d = xy_to_rd(r, d, mosaic, mosaic, debug=debug)
    if debug:
        print(r, d)
    return r, d

def rd_to_xy_iter(r, d, image, mosaic, debug=False):
    r = r.astype('f8')
    d = d.astype('f8')
    if np.isnan(image['crval1']) or np.isnan(mosaic['crval1']):
        raise ValueError('Image or mosaic has no astrometric solution; probably missing chip')

    # First reverse the global r,d
    dis = (mosaic['ctype'][-3:] == 'DIS')
    if not dis:
        raise ValueError('Need to use global WCS.')
    l, m   = rd_to_lm(r, d, mosaic)
    xp, yp = lm_to_xy(l, m, mosaic)

    # Then reverse the local terms
    warp = (image['ctype'][-3:] == 'WRP')
    # do some check with warp

    l, m   = rd_to_lm(xp, yp, image)
    x2, y2 = lm_to_xy(l, m, image)

    if debug:
        print(x2, y2)
    return x2, y2

# This routine needs to be a little different
# from just reversing xy_to_rd.  The routine needs
# to determine which chip the detection would land in,
# if any.
# 1. Reverse the global r,d
# 2. Reverse the global mosaic
# 3. Determine which chip the detection lands in.
# 4. Reverse the local r, d for that chip.
# 5. Reverse the local terms for that chip.
#
def rd_to_xy_chips(r, d, smfcoords):
    r = r.astype('f8')
    d = d.astype('f8')
    mosaic = smfcoords['PRIMARY']
    if np.isnan(mosaic['crval1']):
        raise ValueError('Image or mosaic has no astrometric solution; probably missing chip')

    # First reverse the global r,d
    dis = (mosaic['ctype'][-3:] == 'DIS')
    if not dis:
        raise ValueError('Need to use global WCS.')
    l, m   = rd_to_lm(r, d, mosaic)
    xp, yp = lm_to_xy(l, m, mosaic)

    # At this step, either at lm or at the global xy,
    # the chip needs to be determined.

    names = ps1_utils.get_chip_names(xp, yp)
    chip_dict = defaultdict(list)
    for (idx, name), x, y in zip(names, xp, yp):
        chip_dict[name].append((idx, x, y))

    for k in chip_dict.keys():
        idx = np.array([t[0] for t in chip_dict[k]])
        x   = np.array([t[1] for t in chip_dict[k]])
        y   = np.array([t[2] for t in chip_dict[k]])
        chip_dict[k] = idx, x, y
    
    # Then reverse the local terms
    #warp = (image['ctype'][-3:] == 'WRP')
    # do some check with warp

    results_dict={}
    for k, v in chip_dict.items():
        if k == 'Corner':
            results_dict[k]=v
        else:
            idx, xp, yp = v
            image = smfcoords[k]
            l, m   = rd_to_lm(xp, yp, image)
            x2, y2 = lm_to_xy(l, m, image)
            results_dict[k] = (idx, x2, y2)

    return chip_dict, results_dict

def xy_to_lm(xpix, ypix, image):
    # 1. stretch
    xs = image['cdelt1']*(xpix-image['crpix1'])
    ys = image['cdelt2']*(ypix-image['crpix2'])

    # 2. rotation
    xr = (xs*image['pc1_1']+ys*image['pc1_2'])
    yr = (xs*image['pc2_1']+ys*image['pc2_2'])

    x, y = xs, ys
    l, m = xr, yr

    # 3. polynomial coorection
    npoly = image['npolyterms']
    pterms = image['polyterms'].reshape(7, 2)
    #npoly=0
    # the pterms are 
    if npoly > 1:
        l += x*x*pterms[0,0]+x*y*pterms[1,0]+y*y*pterms[2,0]
        m += x*x*pterms[0,1]+x*y*pterms[1,1]+y*y*pterms[2,1]
    if npoly > 2:
        l += x*x*x*pterms[3,0]+x*x*y*pterms[4,0]+x*y*y*pterms[5,0]+y*y*y*pterms[6,0]
        m += x*x*x*pterms[3,1]+x*x*y*pterms[4,1]+x*y*y*pterms[5,1]+y*y*y*pterms[6,1]
    return l, m

def polycorr(x, y, image):
    # These are the higher-order corrections to the pixel
    # coordinates.
    npoly = image['npolyterms']
    pterms = image['polyterms'].reshape(7, 2)
    # the pterms are
    #npoly=0
    dx, dy = 0.0, 0.0
    if npoly > 1:
        dx += x*x*pterms[0,0]+x*y*pterms[1,0]+y*y*pterms[2,0]
        dy += x*x*pterms[0,1]+x*y*pterms[1,1]+y*y*pterms[2,1]
    if npoly > 2:
        dx += x*x*x*pterms[3,0]+x*x*y*pterms[4,0]+x*y*y*pterms[5,0]+y*y*y*pterms[6,0]
        dy += x*x*x*pterms[3,1]+x*x*y*pterms[4,1]+x*y*y*pterms[5,1]+y*y*y*pterms[6,1]
    return dx, dy

# MJH: Added this April 2019 
def invert_polycorr(xp, yp, image, niters=3):
    x, y = xp, yp
    for i in range(niters):
        dx, dy = polycorr(x, y, image)
        x = xp - dx
        y = yp - dy
    return x, y

# MJH: Added this April 2019 
def invert_rotation(xp, yp, image):
    a = image['pc1_1']
    b = image['pc1_2']
    c = image['pc2_1']
    d = image['pc2_2']    
    inv_det = 1.0/np.abs(a*d-b*c)
    cp1_1 =  inv_det*d
    cp1_2 = -inv_det*b
    cp2_1 = -inv_det*c
    cp2_2 =  inv_det*a
    x = xp*cp1_1 + yp*cp1_2
    y = xp*cp2_1 + yp*cp2_2
    return x, y

# MJH: Added this 17 April 2019 
def invert_rot_polycorr(l, m, image, niters=3):
    dx, dy = 0.0, 0.0
    for i in range(niters):
        x, y = invert_rotation(l-dx, m-dy, image)
        dx, dy = polycorr(x, y, image)
    return x, y

# MJH: Added this 1 April 2019 
def lm_to_xy(l, m, image, niters=10):
    l = np.atleast_1d(l)
    m = np.atleast_1d(m)

    # Invert the rotation and poly corrections
    xs, ys = invert_rot_polycorr(l, m, image, niters=niters)

    # Reverse the stretch
    x = xs/image['cdelt1'] + image['crpix1']
    y = ys/image['cdelt2'] + image['crpix2']

    return x, y
        
def lm_to_rd(l, m, image):
    dis = (image['ctype'][-3:] == 'DIS')
    if not dis:
        ra  = l + image['crval1']
        dec = m + image['crval2']
    else:
        radeg = 180./np.pi
        r = np.sqrt(l**2.+m**2.)
        sphi =  l/(r+(r == 0)) # should go to 0 when r=0
        cphi = (-m+(r == 0))/(r+(r == 0))  # should go to 1 when r=0

        t = (radeg/(r+(r == 0)))*(r != 0)
        stht = (t+(t == 0))/np.sqrt(1.+t**2.)  # should go to 1 when t=0
        ctht = (1./np.sqrt(1.+t**2.))*(t != 0) # should go to 0 when t=0

        sdp  = np.sin(image['crval2']/radeg)
        cdp  = np.cos(image['crval2']/radeg)
      
        sdel = stht*sdp - ctht*cphi*cdp
        salp = ctht*sphi
        calp = stht*cdp + ctht*cphi*sdp
        alpha = np.arctan2(salp, calp)
        delta = np.arcsin(sdel)
      
        ra  = radeg*alpha + image['crval1']
        dec = radeg*delta
      
        ra = (ra % 360.)
    return ra, dec

# Added by Matt Holman 21 March 2019
# Following Gnomonic (Tangent Plane) projection
# https://lambda.gsfc.nasa.gov/product/iras/coordproj.cfm
def rd_to_lm(r, d, image):
    dis = (image['ctype'][-3:] == 'DIS')
    if not dis:
        l   = r - image['crval1']        
        m   = d - image['crval2']        
    else:
        radeg = 180./np.pi

        alpha = r/radeg
        delta = d/radeg
    
        alpha0 = image['crval1']/radeg
        delta0 = image['crval2']/radeg

        A = np.cos(delta)*np.cos(alpha-alpha0)
        scale = 1.0 # I believe this what we want here.
        F = scale/(np.sin(delta0)*np.sin(delta) + A*np.cos(delta0))
        l = -F*np.cos(delta)*np.sin(alpha-alpha0)
        m = -F * (np.cos(delta0) * np.sin(delta) - A*np.sin(delta0))
        l *= -radeg
        m *= -radeg
    return l, m

def xy_to_rd_smfcoords(x, y, chip, smfcoords):
    return xy_to_rd(x, y, smfcoords[chip], smfcoords['PRIMARY'])

def rd_to_xy_smfcoords(r, d, chip, smfcoords):
    return rd_to_xy(r, d, smfcoords[chip], smfcoords['PRIMARY'])

# MJH modified this so that it can optionally take an hdu list
# from a FITS file that has already been opens.  It's not a big
# change, but it saves time.
def readsmfcoords(smffile, hdulist=None):
    out = { }
    if hdulist==None:
        hdulist = fits.open(smffile)
    for i, hdu in enumerate(hdulist):
        #print
        hdu.name
        if ((hdu.name == 'PRIMARY') or
            (re.match(r'XY..\.HDR', hdu.name.strip()) is not None)
            or (re.match(r'XY..\.hdr', hdu.name.strip()) is not None)):
            name = hdu.name if hdu.name == 'PRIMARY' else hdu.name[0:-4]
            #print name
            out[name] = get_ps1_coords(hdu.header)
        if (hdu.name == 'XY__'): # screwed up name?
            name = (hdulist[i+1].name)[0:-4]
            out[name] = get_ps1_coords(hdu.header)
    return out

# stolen from internet, Simon Brunning
def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        files2 = [os.path.join(os.path.relpath(path, start=root), f)
                  for f in files]
        for filename in fnmatch.filter(files2, pattern):
            yield os.path.join(os.path.abspath(root), filename)

def get_smflist(smfdir=None):
    smfdir = (smfdir if smfdir is not None
              else os.getenv('PS_DATA_SMF_DIR', None))
    if not smfdir:
        raise ValueError('Must set PS_DATA_SMF_DIR')
    smf = list(locate('*.smf', smfdir))
    return list(smf)

def get_smf_filename(filename_base):
    smfdir = os.getenv('PS_DATA_SMF_DIR', None)
    if not smfdir:
        raise ValueError('Must set PS_DATA_SMF_DIR')
    matches = list(locate('*/'+filename_base+'*.smf', smfdir))
    if len(matches) == 0:
        raise ValueError("Could not find smf file for %s." % filename_base)
    if len(matches) > 1:
        print(matches)
        raise ValueError("Multiple possible smf files for %s." % filename_base)
    return matches[0]
