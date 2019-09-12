import numpy
from obs80 import obscode
from obs80 import obs80
from obs80 import leapsec
from obs80 import obs80hc


def test_hc():
    hc = obs80hc.HCAnnotator(
        obscode.siteXYZ(obscode.read5c('obscode.dat')),
        leapsec.LeapSeconds('leap-seconds.list'))
    o = obs80.parseOpt('00433        2C2012 05 02.13978010 34 59.368-23 08 51.92         10.89Vg~15MD689')  # noqa E501
    xyz = hc._xyzOpt(o)
    # test values from JPL Horizons.  They pass with isclose defaults, but
    # there's still some discrepancy.
    assert numpy.isclose(xyz[0], -7.494457783619181E-01)
    assert numpy.isclose(xyz[1], -6.182950688586435E-01)
    assert numpy.isclose(xyz[2], -2.680185245982367E-01)
