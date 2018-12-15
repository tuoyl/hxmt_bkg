from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import os
import commands
import sys
import time



if len(sys.argv)>=2:
    sp_lc_select  = sys.argv[1]
    sl_name       = sys.argv[2]
    outnam        = sys.argv[3]


def write_bkgspec(fname,channel,counts,expo,hdr_ext):
    spec_qua= np.zeros(256)
    spec_grp= np.ones(256)
    spec_col1 = pf.Column(name='CHANNEL', format='J', array=channel)
    spec_col2 = pf.Column(name='COUNTS', format='J', array=counts)
    spec_col3 = pf.Column(name='QUALITY', format='I', array=spec_qua)
    spec_col4 = pf.Column(name='GROUPING', format='I', array=spec_grp)
    cols = pf.ColDefs([spec_col1, spec_col2, spec_col3, spec_col4])
    hdr = pf.Header()
    hdr['EXTNAME']  = "SPECTRUM"
    hdr['PHAVERSN'] = "1992a"
    hdr['HDUCLASS'] = "OGIP"
    hdr['HDUCLAS1'] = "SPECTRUM"
    hdr['HDUCLAS2'] = "TOTAL"
    hdr['HDUCLAS3'] = "COUNT"
    hdr['HDUCLAS4'] = "TYPE:I"
    hdr['HDUVERS1'] = "1.2.1"
    hdr['CHANTYPE'] = "PI"
    hdr['DETCHANS'] = hdr_ext['DETCHANS']
    hdr['TELESCOP'] = 'HXMT'
    hdr['INSTRUME'] = 'HE'

    hdr["OBS_MODE"] = hdr_ext['OBS_MODE']
    hdr["DATE-OBS"] = hdr_ext['DATE-OBS']
    hdr["DATE-END"] = hdr_ext['DATE-END']
    hdr["OBJECT"]   = hdr_ext['OBJECT']
    hdr["TSTART"]   = hdr_ext['TSTART']
    hdr["TSTOP"]    = hdr_ext['TSTOP']
    hdr["RA_OBJ"]   = hdr_ext['RA_OBJ']
    hdr["DEC_OBJ"]  = hdr_ext['DEC_OBJ']

    hdr['CORRFILE'] = 'None'
    hdr['CORRSCAL'] = 1.0

    hdr['BACKFILE'] = 'NONE'
    hdr['BACKSCAL'] = 1.0

    hdr['RESPFILE'] = 'NONE'
    hdr['ANCRFILE'] = 'NONE'
    hdr['FILTER']   = 'NONE'

    hdr['AREASCAL'] = 1.0
    hdr['EXPOSURE'] = expo
    hdr['LIVETIME'] = expo
    hdr['DEADC']    = 1.0
    hdr['STATERR']  = True
    hdr['SYSERR']   = False
    hdr['POISSERR'] = True
    hdr['GROUPING'] = 0
    hdr['QUALITY']  = 0
    hdr['MJDREFI']  = hdr_ext['MJDREFI']
    hdr['MJDREFF']  = hdr_ext['MJDREFF']
    hdu = pf.BinTableHDU.from_columns(cols,header=hdr)
    hdu.writeto(fname,overwrite=True)

def write_lcurve(fname,time,counts,hdr_ext):
    lc_frac= np.ones(np.size(time))
    lc_col1 = pf.Column(name='Time', format='D', array=time)
    lc_col2 = pf.Column(name='Counts', format='J', array=counts)
    lc_col3 = pf.Column(name='FRACEXP', format='D', array=lc_frac)
    cols = pf.ColDefs([lc_col1, lc_col2, lc_col3])
    hdr = pf.Header()
    hdr['EXTNAME']  = "RATE"
    hdr['PHAVERSN'] = "1992a"
    hdr['HDUCLASS'] = "OGIP"
    hdr['HDUCLAS1'] = "LIGHTCURVE"
    hdr['HDUCLAS2'] = "ALL"
    hdr['HDUCLAS3'] = "COUNT"
    hdr['HDUCLAS4'] = "TYPE:I"
    hdr['HDUVERS1'] = "1.1.0"
    hdr['CHANTYPE'] = "PI"
    hdr['TELESCOP'] = hdr_ext['TELESCOP']
    hdr['INSTRUME'] = hdr_ext['INSTRUME']
    hdr['TIMEUNIT'] = hdr_ext['TIMEUNIT']
    hdr['MJDREFI']  = hdr_ext['MJDREFI']
    hdr['MJDREFF']  = hdr_ext['MJDREFF']
    hdr["OBS_MODE"] = hdr_ext['OBS_MODE']
    hdr["DATE-OBS"] = hdr_ext['DATE-OBS']
    hdr["DATE-END"] = hdr_ext['DATE-END']
    hdr["OBJECT"]   = hdr_ext['OBJECT']
    hdr["TSTART"]   = hdr_ext['TSTART']
    hdr["TSTOP"]    = hdr_ext['TSTOP']
    hdr["RA_OBJ"]   = hdr_ext['RA_OBJ']
    hdr["DEC_OBJ"]  = hdr_ext['DEC_OBJ']
    hdr['TIMEZERO']  = hdr_ext['TIMEZERO']
    hdr['TIMEDEL']  = hdr_ext['TIMEDEL']

    hdu = pf.BinTableHDU.from_columns(cols,header=hdr)
    hdu.writeto(fname,overwrite=True)



src_name=[]
sf = open(sl_name)
for sline in sf:
    src_name.append(sline)
sf.close()
#print(src_name)

dnum = np.size(src_name)
cnt_all = np.zeros(16384)
expo_all = 0
for ii in xrange(0,dnum):
    tmpstr = src_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For detector:", ii, " src file name: ", tmpstr[0:tmppos])
    spec_list    = pf.open(tmpstr[0:tmppos])
    spec_tab     = spec_list[1].data
    spec_ch = spec_tab.field(0)
    spec_counts  = spec_tab.field(1)
    spec_hdr     = spec_list[1].header
    spec_list.close()
    cnt_all[0:np.size(spec_counts)]  = cnt_all[0:np.size(spec_counts)] + spec_counts
    expo_all = expo_all + spec_hdr['exposure'] 

cnt_all = cnt_all[0:np.size(spec_counts)]
outname = outnam +'.pha'
write_bkgspec(outname,spec_ch,cnt_all,expo_all,spec_hdr)


plt.figure()
plt.plot(cnt_all)
plt.show()

