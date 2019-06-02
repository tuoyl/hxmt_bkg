#!coding=utf-8
'''
Model constructed by Background Group.
Lian Jinyuan, Zhangshu, Guo Chengcheng, Jin Jing, Zhangjuan, Zhang Shu, et al.
Mail liaojinyuan@ihep.ac.cn

'''
'''
This version was written by Ge Mingyu
Mail gemy@ihep.ac.cn

Usage:

hegti ehkfile.fits oldgtifile.fits newgti.fits

Using interactive method in prompt.	

'''
from astropy.io import fits as pf
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import os
import sys
import time
from scipy import interpolate

try:
    # Python 2
    xrange
except NameError:
    # Python 3, xrange is now named range
    xrange = range

print( "*********************************************************" )
print( "*********************************************************" )
print( "*********************************************************" )
print( "************ PRINT: hegti -h for usage   *************" )
print( "*********************************************************" )

uage_method1 = 'Method 1: hegti ehkfile.fits gtifile.fits newfile.fits'
uage_method2 = 'Method 2: Using interactive method in prompt.'

def print_usage(uage_method1,uage_method2):
    print(uage_method1)
    print(uage_method2)

if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print_usage(uage_method1,uage_method2)
    sys.exit()
elif len(sys.argv)>=2:
    ehkname       = sys.argv[1]
    gtifile       = sys.argv[2]
    gtioutname   = sys.argv[3]
else:
    ehkname     = str(raw_input("EHK file:"))
    gtifile     = str(raw_input("GTI file:"))
    gtioutname = str(raw_input("New gtifile:"))


#HEADAS=os.getenv('HEADAS')
HEADAS=os.getenv('REFPATH')
REFPATH=HEADAS+'/'


'''Check the time range in GTI'''
def is_ingti(START,STOP,tl,tu):
    num = len(START)
    is_gti = 0
    for ii in xrange(0,num):
        t0=START[ii]
        t1=STOP[ii]
        flag0 = (tl<=t0) & (tu>=t0)
        flag1 = (tl>=t0) & (tu<=t1)
        flag2 = (tl<=t1) & (tu>=t1)
        if flag0 | flag1 | flag2:
            is_gti = 1
            return is_gti
    return is_gti
'''Check the time in GTI '''
def is_ingti2(START,STOP,ti):
    num = np.size(START)
    is_gti = 0
    if num == 1:
        t0=START
        t1=STOP
        flag0 = (ti>=t0) & (ti<=t1)
        if flag0:
            is_gti = 1
            return is_gti
    if num >= 2:
        for ii in xrange(0,num):
            t0=START[ii]
            t1=STOP[ii]
            flag0 = (ti>=t0) & (ti<=t1)
            if flag0:
                is_gti = 1
                return is_gti
    return is_gti
'''Give the tag for specific time'''
def time_gtiflag(time,START,STOP,tflag):
    for ii in xrange(0,len(time)):
        tflag[ii] = is_ingti2(START,STOP,time[ii])

'''Give the index for specific time'''
def flag_selection(tflag,tindex):
    cnt = 0
    for ii in xrange(0,len(tflag)):
        if (tflag[ii] == 1):
            tindex[cnt] = ii
            cnt = cnt+1
    return cnt

'''Give the index for bkg map'''
def bkgmap_index(LON,LAT,ascend_flag ,step):
    lon_index = int(LON/step)
    lat_index = int((LAT+45)/step)
    map_index = lon_index+lat_index*72+ascend_flag*1296
    return map_index

'''Give the index for bkg map for an array'''
def bkgmap_time(lon_arr,lat_arr,ascend_flag,bkgmap_flag):
    fnum = np.size(lon_arr)
    for ii in xrange(0,fnum):
        bkgmap_flag[ii] = bkgmap_index(lon_arr[ii],lat_arr[ii],ascend_flag[ii],5.)

'''Give the data points for specific 5x5 size'''
def bkgmap_num(bkgmap_flag):
    dval = bkgmap_flag[1:(len(bkgmap_flag))] - bkgmap_flag[0:(len(bkgmap_flag)-1)]
    fnum = np.size((np.where(np.abs(dval) > 0))) + 1
    print("5x5 interval number:",fnum)
    return fnum

'''Give the start and stop indices for specific 5x5 size'''
def bkgmap_interval_index(bkgmap_flag,bkgmap_start_index,bkgmap_stop_index):
    fnum = np.size(bkgmap_flag)
    inum = np.size(bkgmap_start_index)
    dval = bkgmap_flag[1:fnum] - bkgmap_flag[0:(fnum-1)]
    nonzero_in = np.where(np.abs(dval) > 0)
    nonzero_num = np.size(nonzero_in)
    print( 'Non zeros index',nonzero_in[0:(fnum)] )
    bkgmap_start_index[1:(inum)] = nonzero_in + np.ones(nonzero_num,dtype=np.int)
    bkgmap_start_index[0] = 0
    #print( np.size(bkgmap_stop_index[0:(inum-2)]),np.size(nonzero_in) )
    bkgmap_stop_index[0:(inum-1)] =nonzero_in + np.zeros(nonzero_num,dtype=np.int)
    bkgmap_stop_index[inum-1] = fnum-1

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
    lc_col1 = pf.Column(name='Time', format='D', unit='s', array=time)
    lc_col2 = pf.Column(name='Rate', format='D', unit='counts/s', array=counts)
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

def creatnewgti(gtiold, gtiout, tstart, tstop):
    #copy imformation from old gtifile
    old = pf.open(gtiold)
    oldheader0 = old[0].header
    oldheader1 = old[1].header
    oldheader2 = old[-1].header
    olddata2 = old[-1].data
    old.close()
    ######## Create new FITS
    ######## extension 0
    primarydata = pf.PrimaryHDU(header=oldheader0)
    ######## extension 1 2
    c11 = pf.Column(name='START', array=tstart, format='1D')
    c12 = pf.Column(name='STOP', array=tstop, format='1D')
    data1 = pf.BinTableHDU.from_columns([c11,c12])
    c21=pf.Column(name='a', array=[0,0],format='1D')
    c22=pf.Column(name='b', array=[0,0],format='1D')
    data2 = pf.BinTableHDU.from_columns([c21,c22])
    ########
    newfits = pf.HDUList([primarydata,data1,data2])
#    newfits.writeto(gtiout,overwrite=True) #for server
    newfits.writeto(gtiout,overwrite=True)
    ######## fits Write new data
    naxis2 = int(len(tstart))
    oldheader1['NAXIS2'] = naxis2
    newfits = pf.open(gtiout,mode='update')
    newfits[0].header = oldheader0
    newfits[1].header = oldheader1
    newfits[2].header = oldheader2
    newfits[1].header['HISTORY'] = 'correct gti by legti.py'
    newfits[1].header['EXTNAME'] = 'GTI0'
    newfits[2].header['HISTORY'] = 'correct gti by legti.py'
    newfits[2].data = olddata2
    x1=olddata2.field(1)
    x1*=0
    newfits.close()
    print('done')
   #^^^^^^^END


'''Read gti extesion'''

hdulist = pf.open(gtifile)
tb = hdulist[1].data
START = tb.field(0)
STOP = tb.field(1)
hdulist.close()

'''Read EHK file and position devision by 5 x 5'''
ehklist = pf.open(ehkname)
ehk_tab = ehklist[1].data
ehk_time= ehk_tab.field(0)
ehk_LON = ehk_tab.field(8)
ehk_LAT = ehk_tab.field(9)
ehklist.close()

ehk_num = np.size(ehk_time)
#Ascend: 0
#Descend: 1
ehk_diff = ehk_LAT[1:(ehk_num)] - ehk_LAT[0:(ehk_num-1)] 
ehk_ascend_flag = np.zeros(ehk_num);
ehk_ascend_flag[np.where(ehk_diff >= 0)] = 0
ehk_ascend_flag[np.where(ehk_diff < 0)]  = 1
ehk_ascend_flag[ehk_num-1]  = ehk_ascend_flag[ehk_num-2]

tflag = np.zeros(ehk_num,dtype='int')
tindex = np.zeros(ehk_num,dtype='int')
time_gtiflag(ehk_time,START,STOP,tflag)
ehk_num = flag_selection(tflag,tindex)
ehk_time = ehk_time[tindex[0:(ehk_num)]]
ehk_LON = ehk_LON[tindex[0:(ehk_num)]]
ehk_LAT = ehk_LAT[tindex[0:(ehk_num)]]
ehk_ascend_flag= ehk_ascend_flag[tindex[0:(ehk_num)]]

bkgmap_flag = np.zeros(ehk_num,dtype='int')
bkgmap_time(ehk_LON,ehk_LAT,ehk_ascend_flag,bkgmap_flag)
bkgmap_flag = np.array(bkgmap_flag)

bkgspec_ID0 = np.zeros(ehk_num,dtype='float')

print(ehk_num)

'''Read background model'''
srcmapname=REFPATH+'HE_bkgmap.fits'
srclist = pf.open(srcmapname)
srcmap_tab = srclist[1].data
srcmap_IN  = srcmap_tab.field(0)
srcmap_LON = srcmap_tab.field(1)
srcmap_LAT = srcmap_tab.field(2)
srcmap_BKG = srcmap_tab.field(3)

'''Cal the spectra of all detectors from map'''
detid = 0
tin1 = detid*256+26
tin2 = (detid+1)*256
for ii in xrange(0,ehk_num):
    tmpindex = bkgmap_flag[ii]
    tmpspec  = srcmap_BKG[tmpindex,0:4608]
    bkgspec_ID0[ii] = np.sum(tmpspec[tin1:tin2])

'''Create new GTI'''

#plt.figure()
#plt.plot(ehk_time-ehk_time[0],bkgspec_ID0)
#plt.show()

din = np.where(bkgspec_ID0 <= 36)
din = np.array(din[0],dtype='int')
dtime = ehk_time[din]
din_num = np.size(dtime)

if din_num==0:
    print("No GTI!!!!")
    sys.exit()

ddtime = dtime[1:din_num] - dtime[0:(din_num-1)]
din2 = np.where(ddtime >= 30)
din2 = np.array(din2,dtype='int')
din2_num = np.size(din2)
int_num = din2_num + 1
print(int_num)
NSAA_START = np.zeros(int_num,dtype='int')
NSAA_STOP  = np.zeros(int_num,dtype='int')
NSAA_START[0] = 0
NSAA_STOP[int_num-1] = din_num-1
if (int_num>=2):
    NSAA_START[1:(int_num)]  = din2 + 1
    NSAA_STOP[0:(int_num-1)] = din2

gtistart1 = dtime[NSAA_START]
gtistop1  = dtime[NSAA_STOP]

print(STOP,START)
print(gtistop1,gtistart1)
creatnewgti(gtifile, gtioutname, gtistart1, gtistop1)


