#!/usr/bin/python
'''
Generate the spectrum with time bin time 16s or 32 seconds
'''
'''
Usage:

python hebldgen.py blind_det16.FITS bldarr.fits DTime.FITS ehkname gtifile

'''

from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time

print "*********************************************************"
print "******************Running HXMT Bkg******************"
print "*********************************************************"

if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print "Example 1: hgti_abs lc.fits gti.txt out.txt"
        print "Example 2: Using interactive in prompt."
    sys.exit()
elif len(sys.argv)>=2:
    lcname=sys.argv[1]
    ingtiname=sys.argv[2]
    outgtifile=sys.argv[3]
else:
    lcname     = str(raw_input("Light curve: "))
    ingtiname  = str(raw_input("Relative time range file: "))
    outgtifile    = str(raw_input("Output GTI file:"))

hdulist = pf.open(lcname)
tb = hdulist[1].data
lctime = tb.field(0)
counts = tb.field(1)
hdulist.close()
lctime_t0 = np.min(lctime)

#plt.figure()
#plt.plot(lctime-lctime[0],counts)
#plt.show()
'''Read relative time range'''
gti0 = np.loadtxt(ingtiname)
gti1 = gti0 + lctime_t0

np.savetxt(outgtifile,gti1)

print gti0,gti1





