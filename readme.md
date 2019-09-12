# obs80

Low level reader and parser for 80 column observation format.  Expects
relatively clean 80 column records.

## Prerequisites

### Python modules
* numpy
* novas
* pandas, used only in the test code.

### JPL Ephemeris

Module obs80hc uses novas functions that require access to a binary JPL
ephemeris.  See for example the
[novas readme](http://aa.usno.navy.mil/software/novas/novas_py/README)
in the section "Ephemerides".  Also, while currently undocumented, it appears
eph_manager.ephem_open() can find the ephemeris file in the package
`novas_de405` if that package is installed.

### File leap-seconds.list

A leap-seconds.list file is required as a function argument.  The format is
that prepared by Judah Levine at NIST.  An example source for this document is
https://www.ietf.org/timezones/data/leap-seconds.list.

### File obscode.dat

An obscode.dat file is required as a function argument.  Currently supported
is the current style 5 column obscode.dat format, the format currently publicly
accessible at http://www.minorplanetcenter.net/iau/lists/ObsCodes.html.

### Sphinx documentation

The Sphinx documentation preparation tool can be used to prepare html (or
other format) documentation from docstrings in the Python source.

* [Sphinx](http://www.sphinx-doc.org/) 1.2.3 is installed on mpcapp1. 
  It's a pretty old version and isn't working so well.   
* The file conf.py used by Sphinx loads the numpydoc extension. 
  See [NumPy HowTo](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt)
* Generate Sphinx documentation with the the command `make`, then try
  `doc/index.html` in a browser.

## Tests
Tests are written for the [pytest](https://docs.pytest.org/en/latest/) tool.

The obs80hc tests require a JPL ephemeris as mentioned above.  Thus you must
either have the package `novas_de405` installed into the environment that
pytest will see or else set the environment variable EPHEMERIS\_FILE.

Then to run tests, With pytest installed, just type `pytest` at the
command line:

```
$ pytest
============================= test session starts ==============================
platform linux -- Python 3.4.3 -- py-1.4.30 -- pytest-2.7.3
rootdir: /home/skeys/mpcdev/obs80, inifile: 
plugins: celery
collected 18 items 

leapsec_test.py .
obs80_test.py ..............
obs80hc_test.py .
obscode_test.py ..

========================== 18 passed in 1.43 seconds ===========================
```
