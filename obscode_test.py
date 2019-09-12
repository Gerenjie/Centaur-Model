import pytest

from obs80 import obscode


@pytest.fixture(scope='module')
def sites():
    # look a few places for a real obcode file
    for fn in ['obscode.dat', 'ObsCodes.html',
            '/var/www/html/mpc/public/iau/lists/ObsCodes.html']:
        try:
            return obscode.read5c(fn)
        except IOError:
            pass
    # otherwise use small sample file in source control
    return obscode.read5c('ocsample.dat')


def test_sites(sites):
    assert sites['G96'] == obscode.Site5c(blank=False, long=249.21128,
        rcos=0.845111, rsin=0.533614, site='Mt. Lemmon Survey')
    # Blank parallax constants
    assert sites['248'] == obscode.Site5c(blank=False, long=0.0, rcos=0.0,
        rsin=0.0, site='Hipparcos')
    # 0 parallax constants (not blank)
    assert sites['248'] == obscode.Site5c(blank=False, long=0.0, rcos=0.0,
        rsin=0.0, site='Hipparcos')
    # Southern hemisphere
    assert sites['413'] == obscode.Site5c(blank=False, long=149.06608,
        rcos=0.855595, rsin=-0.516262, site='Siding Spring Observatory')


def test_xyz(sites):
    xyz = obscode.siteXYZ(sites)
    assert xyz['W14'] == (0.04652656088240064, -0.8205960653868969, 0.567774)
    assert xyz['250'] == (None, None, None)
    assert xyz['248'] == (0.0, 0.0, 0.0)
    assert xyz['413'] == (-0.7338957925529292, 0.43981788242192676, -0.516262)
