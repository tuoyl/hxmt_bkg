#!/usr/bin/env python3
'''
Model constructed by Background Group.
Lian Jinyuan, Zhangshu, Guo Chengcheng, Jin Jing, Zhangjuan, Zhang Shu, et al.
Mail liaojinyuan@ihep.ac.cn

'''
'''
This version was written by Ge Mingyu
Mail gemy@ihep.ac.cn

Usage:

lebkgmap lc/spec screen.fits ehk.fits gtifile.fits lcname/specname chmin chmax outnam_prefix
    lc/spec: lc for background lightcurve and spec for background light curve 
    screen.fits: screened event file name, which should include the events for all blind detectors.
    gtifile.fits: the GTI file for LE
    lcname/specname: is ASCII file, which includes the name of the source file for small FOV
    chmin: minimum channel for light curve
    chmax: maximum channel for light curve
    outnam_prefix: the output prefix for the spectrum

lebkgmap sflag=lc/spec evtfile=screen.fits gtifile=gtifile.fits srcdat=lcname/specname chmin=chmin chmax=chmax outnam=outnam_prefix
    sflag: lc/spec, lc for background lightcurve and spec for background light curve 
    evtfile: screened event file name, which should include the events for all blind detectors.
    gtifile: the GTI file for LE
    srcdat: is ASCII file, which includes the name of the source file
    chmin:  minimum channel for light curve
    chmax:  maximum channel for light curve
    outnam: the output prefix for the spectrum


Using interactive method in prompt.	

History:

2019-04-10
Using the particle background (shape is not changing)
CXB is located in |b|>10
GXB is located in |b|<=10

The devision by LAT LON in EHK file is not considered in this version

The correction in the final step is canceled.
Upadte the error estimation.

2019-04-26
GXB map search Error:

From 'src_index1 = int(ll_src+5-ii)*20+int((bb_src+10+5.5-jj))'
to 'src_index1 = int(ll_src+5-ii)+int((bb_src+10+5.5-jj))*360'

GXB map region: 11x11

The sequence of GXBMAP save as: bb,ll

Update the exposure of GXBMAP for three LE boxes. 

"Spec_Model.txt"     --> "Spec_Model_20190322.txt"
"LE_gxmap_v2.5.fits" --> LE_gxmap_v2.6.fits

spec_gxb[1050:1536] = 0

2019-05-09:
Channel 1500-1535: set 0


'''

from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import time
from scipy.optimize import curve_fit
Ver = '2.0.7'

print( "*********************************************************" )
print( "******************  Running HXMT Bkg   ******************" )
print( "*********************************************************" )
print( "*********************************************************" )
print( "*********************************************************" )
print( "************ PRINT: lebkgmap -h for usage   *************" )
print( "*********************************************************" )
print( "HXMT background for Insight-HXMT/LE, ver-",Ver )
print( "The energy range for background lightcurve should be the same as source lightcurve")

uage_method1 = 'Method 1: lebkgmap lc/spec screen.fits gtifile.fits lcname/specname chmin chmax outnam_prefix'
uage_method2 = 'Method 2: Using interactive method in prompt.'
uage_method3 = 'Method 3: lebkgmap sflag=lc/spec evtfile=screen.FITS gtifile=gtifile.fits srcdat=lcname/specname chmin=chmin chmax=chmax outnam=outnam_prefix'

def print_usage(uage_method1,uage_method2,uage_method3):
    print(uage_method1)
    print(uage_method2)
    print(uage_method3)

def check_argument():
    sp_lc_select  = []
    evtfilename   = []
    gtifile       = []
    sl_name       = []
    chmin         = []
    chmax         = []
    outnam        = []
    slgti         = []
    len_arg = len(sys.argv)
    if len_arg <=1:
        raise IOError("Error input argument, RUN 'lebkgmap -h' for help")
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print_usage(uage_method1,uage_method2)
        return False, False, False, False
    if((len_arg>1)&(len_arg<8)):
        raise IOError("Error input argument, RUN 'lebkgmap -h' for help")
        sys.exit()
    sysflag = 0
    for i in xrange(len_arg):
        arg = sys.argv[i]
        if i == 0:continue
        arg_split = arg.split('=')
        if (len(arg_split)==1):continue
        argname = arg_split[0]
        argval  = arg_split[1]
        if argname == 'sflag': 
            sp_lc_select.append(argval)
            sysflag = sysflag + 1
        if argname == 'evtfile': 
            evtfilename.append(argval)
            sysflag = sysflag + 1
        if argname == 'gtifile':
            gtifile.append(argval)
            sysflag = sysflag + 1
        if argname == 'srcdat':
            sl_name.append(argval.strip())
            sysflag = sysflag + 1
        if argname == 'chmin':
            chmin.append(int(argval))
            sysflag = sysflag + 1
        if argname == 'chmax':
            chmax.append(int(argval))
            sysflag = sysflag + 1
        if argname == 'outnam':
            outnam.append(argval)
            sysflag = sysflag + 1
        if (len(sys.argv)==8):
            slgti.append('')
        if (len(sys.argv)>=9):
            if argname == 'newgti':
                slgti.append(argval)
                sysflag = sysflag + 1
    if((sysflag>0)&(sysflag<8)):
        raise IOError("Error input argument, RUN 'lebkgmap -h' for help")
        sys.exit()
    if (sysflag==0):
        sp_lc_select.append(sys.argv[1])
        evtfilename.append(sys.argv[2])
        gtifile.append(sys.argv[3])
        sl_name.append(sys.argv[4])
        chmin.append(int(sys.argv[5]))
        chmax.append(int(sys.argv[6]))
        outnam.append(sys.argv[7])
        if (len(sys.argv)==8):
            slgti.append('')
        if (len(sys.argv)>=9):
            slgti.append(sys.argv[8])
    return sp_lc_select[0],evtfilename[0],gtifile[0],sl_name[0],chmin[0],chmax[0],outnam[0],slgti[0]

if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print_usage(uage_method1,uage_method2,uage_method3)
    sys.exit()
elif len(sys.argv)>=2:
    sp_lc_select,evtfilename,gtifile,sl_name,chmin,chmax,outnam,slgti=check_argument()
else:
    sp_lc_select= str(raw_input("Selection(spec/lc):"))
    evtfilename = str(raw_input("Screened events file:"))
    gtifile     = str(raw_input("GTI file:"))
    sl_name     = str(raw_input("FileName inlude spectra or lightcurve name:"))
    chmin       = int(raw_input("Minimum channel:"))
    chmax       = int(raw_input("Maximum channel:"))
    outnam      = str(raw_input("The prefix of output file name:"))
    slgti       = str(raw_input("Specific time range file(NONE):"))

#HEADAS=os.getenv('HEADAS')
HEADAS=os.getenv('REFPATH')
if HEADAS==None:
    print("Environmental parameter REFPATH does not exist!")
    sys.exit()

REFPATH=HEADAS+'/'

if os.path.exists(REFPATH)==False:
    print("REFPATH does not exist!")
    sys.exit()

'''Check the time range in GTI'''
def is_ingti(START,STOP,tl,tu):
    num = np.size(START)
    is_gti = 0
    for ii in xrange(0,num):
        #print(ii,START)
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
def bkgmap_index(LON,LAT,step):
    lon_index = int(LON/step)
    lat_index = int((LAT+45)/step)
    map_index = lon_index+lat_index*72
    return map_index

'''Give the index for bkg map for an array'''
def bkgmap_time(lon_arr,lat_arr,bkgmap_flag):
    fnum = np.size(lon_arr)
    for ii in xrange(0,fnum):
        bkgmap_flag[ii] = bkgmap_index(lon_arr[ii],lat_arr[ii],5.)

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
    #print( 'Non zeros index',nonzero_in[0:(fnum)] )
    bkgmap_start_index[1:(inum)] = nonzero_in + np.ones(nonzero_num,dtype=np.int)
    bkgmap_start_index[0] = 0
    #print np.size(bkgmap_stop_index[0:(inum-2)]),np.size(nonzero_in)
    bkgmap_stop_index[0:(inum-1)] =nonzero_in + np.zeros(nonzero_num,dtype=np.int)
    bkgmap_stop_index[inum-1] = fnum-1



def write_bkgspec(fname,channel,counts,spec_err,expo,hdr_ext):
    spec_qua= np.zeros(1024)
    spec_grp= np.ones(1024)
    spec_col1 = pf.Column(name='CHANNEL', format='J', array=channel)
    spec_col2 = pf.Column(name='COUNTS', format='D', array=counts)
    spec_col3 = pf.Column(name='STAT_ERR', format='D', array=spec_err)
    spec_col4 = pf.Column(name='QUALITY', format='I', array=spec_qua)
    spec_col5 = pf.Column(name='GROUPING', format='I', array=spec_grp)
    cols = pf.ColDefs([spec_col1, spec_col2, spec_col3, spec_col4, spec_col5])
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
    hdr['INSTRUME'] = 'LE'

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
    hdr['POISSERR'] = False
    hdr['GROUPING'] = 0
    hdr['QUALITY']  = 0
    hdr['MJDREFI']  = hdr_ext['MJDREFI']
    hdr['MJDREFF']  = hdr_ext['MJDREFF']
    hdu = pf.BinTableHDU.from_columns(cols,header=hdr)
    hdu.writeto(fname,overwrite=True)

def write_lcurve(fname,time,counts,lc_err,hdr_ext):
    lc_frac= np.ones(np.size(time))
    lc_col1 = pf.Column(name='Time', format='D', unit='s',array=time)
    lc_col2 = pf.Column(name='RATE', format='D', unit='counts/s', array=counts)
    lc_col3 = pf.Column(name='Error', format='D', array=lc_err)
    lc_col4 = pf.Column(name='FRACEXP', format='D', array=lc_frac)
    cols = pf.ColDefs([lc_col1, lc_col2, lc_col3, lc_col4])
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
#From http://blog.sina.com.cn/s/blog_a046022d0101a6hs.html
# Ra: in degrees
# Dec: in degrees
# Output: in rad

def equatorial2galactic( RA, Dec):

    RA = RA * np.pi / 180
    Dec = Dec* np.pi / 180
       
    RAgp = 192.85948 * np.pi / 180.0
    Decgp = 27.12825 * np.pi / 180.0
    lcp = 122.932 * np.pi / 180.0
   
    sinb = np.sin( Dec ) * np.sin( Decgp ) + np.cos( Dec ) * np.cos( Decgp ) * np.cos( RA - RAgp )
    cosbsinlcp_l = np.cos( Dec ) * np.sin( RA - RAgp )
    cosbcoslcp_l = np.sin( Dec ) * np.cos( Decgp ) - np.cos( Dec ) * np.sin( Decgp ) * np.cos( RA - RAgp )
   
    b = np.arcsin( sinb )
    sinlcp_l = cosbsinlcp_l / np.cos( b )
    coslcp_l = cosbcoslcp_l / np.cos( b )
   
    if ( sinlcp_l > 0 ) : lcp_l = np.arccos( coslcp_l )
    else :
        if ( coslcp_l < 0 ) : lcp_l = 2 * np.pi - np.arccos( coslcp_l )
        else : lcp_l = np.arcsin( sinlcp_l )
    l = lcp - lcp_l
    if ( l < 0 ) : l = 2 * np.pi + l
    return [ l, b ]

def spec_bld(x,N,sig,L):
    y = N*np.exp(-((x-1550)/sig)**2) + L
    return y

def spec_bld2(x,a,b):
    y = a*x**b
    return y

outspec='./lebldspec'
timestep=64.0
chstep = 8
ledetchans=192

spec_ch = np.linspace(0,1535,1536)
spec_ch2= np.linspace(0,ledetchans-1,ledetchans)*chstep + chstep*0.5 
spec_ch2= np.linspace(0,ledetchans-1,ledetchans)*chstep + 3.5 


'''Read gti extesion'''

hdulist = pf.open(gtifile)
tb = hdulist[1].data
START = tb.field(0)
STOP = tb.field(1)
hdulist.close()
print("GTI START=",START)
print("GTI STOP=",STOP)

GTI_num = np.size(START)

if GTI_num<=0:
    print("Empty GTI!")
    sys.exit()

bldstart=np.min(START)
bldstop=np.max(STOP)

lcnum=int((bldstop-bldstart)/timestep)+1

timcen = np.linspace(0,lcnum-1,lcnum)*timestep + 0.5*timestep
timbot = np.linspace(0,lcnum-1,lcnum)*timestep
timrof = np.linspace(0,lcnum-1,lcnum)*timestep + timestep
trnum = np.size(timcen)

'''Read Blind map infomation and file'''
bldinfo_list  = pf.open(REFPATH+'LE_bldmap_info.fits')
bldinfo_tab   = bldinfo_list[1].data
bldmapstart  = bldinfo_tab.field('TSTART')
bldmapstop   = bldinfo_tab.field('TSTOP')
bldmapname  = bldinfo_tab.field('BLDMAP')
bldinfo_list.close()

bldmapnum = np.size(bldmapstart)
bldindex = -1
tmid = (bldstart+bldstop)/2.0
for ii in xrange(bldmapnum):
    if (tmid>=bldmapstart[ii])&(tmid<bldmapstop[ii]):
        bldindex = ii

bldmap_bkgname = REFPATH + bldmapname[bldindex]

if (bldindex==-1):
    print("The time interval for blind detector does not exist!")
    sys.exit()

'''Read Blind Detector map'''
#print(bldmap_bkgname)
srclist = pf.open(bldmap_bkgname[0])
srcmap_tab = srclist[1].data
srcmap_IN  = srcmap_tab.field(0)
srcmap_LON = srcmap_tab.field(1)
srcmap_LAT = srcmap_tab.field(2)
srcmap_BKG = srcmap_tab.field(3)

'''Read Coe information and file'''
coeinfo_list = pf.open(REFPATH+'LE_coemap_info.fits')
coeinfo_tab  = coeinfo_list[1].data
coemapstart  = coeinfo_tab.field('TSTART')
coemapstop   = coeinfo_tab.field('TSTOP')
coemapname   = coeinfo_tab.field('BLDMAP')
coeinfo_list.close()

coemapnum = np.size(coemapstart)
coeindex = -1
for ii in xrange(bldmapnum):
    if (tmid>=coemapstart[ii])&(tmid<coemapstop[ii]):
        coeindex = ii

coemap_bkgname = REFPATH + coemapname[coeindex]

if (coeindex==-1):
    print("The time interval for background does not exist!")
    sys.exit()

'''Read Blind Detector map'''
coelist = pf.open(coemap_bkgname[0])
coemap_tab = coelist[1].data
coemap_IN  = coemap_tab.field(0)
coemap_LON = coemap_tab.field(1)
coemap_LAT = coemap_tab.field(2)
coemap_BKG = coemap_tab.field(3)


'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++Read diffuse background for Galaxy plane+++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

src_name=[]
sf = open(sl_name)
for sline in sf:
    print( sline )
    src_name.append(sline)
sf.close()
tmpstr = src_name[0]
tmppos = tmpstr.find('\n')
spec_list    = pf.open(tmpstr[0:tmppos])
spec_tab     = spec_list[1].data
spec_hdr     = spec_list[1].header
spec_list.close()
ra_src =spec_hdr['RA_OBJ'] 
dec_src =spec_hdr['DEC_OBJ'] 

'''Read Simple model'''
print(REFPATH+'Spec_Model.txt')
#simbkgmod_list= np.loadtxt(REFPATH+'Spec_Model.txt')
simbkgmod_list= np.loadtxt(REFPATH+'Spec_Model_20190322.txt')
simbkgmod_map = simbkgmod_list[:,2]
simbkgmod_maperr = simbkgmod_list[:,3]

'''Read Simple binned model'''
print(REFPATH+'Spec_Model_192.txt')
simbkgmod_list= np.loadtxt(REFPATH+'Spec_Model_192.txt')
simbkgmod_map_192 = simbkgmod_list[:,2]
simbkgmod_maperr_192 = simbkgmod_list[:,3]


fac_all = 3.1548740 + 3.4587939 + 2.8518058;

'''Read CXB'''
#print(REFPATH+'CXB_LE_Whole_BA.txt')
#cxbkg_list= np.loadtxt(REFPATH+'CXB_LE_Whole_BA.txt')
#cxbkg_map = cxbkg_list[:,1]
print(REFPATH+'CXB_LE_BA_rebuild_1536_smooth.txt')
cxbkg_list= np.loadtxt(REFPATH+'CXB_LE_BA_rebuild_1536_smooth.txt')
cxbkg_map = cxbkg_list[:,1]

'''Read Galaxy X-ray Background'''
gxbkg_list= pf.open(REFPATH+'LE_gxmap_v2.6.fits')
gxbkg_tab = gxbkg_list[1].data
gxbkg_in  = gxbkg_tab.field(0)
gxbkg_ll  = gxbkg_tab.field(1)
gxbkg_bb  = gxbkg_tab.field(2)
gxbkg_map = gxbkg_tab.field(3)
gxbkg_expo0= gxbkg_tab.field(4)
gxbkg_expo1= gxbkg_tab.field(5)
gxbkg_expo2= gxbkg_tab.field(6)
gxbkg_list.close()

'''Cal Galaxy coordinates'''
[ll_src,bb_src] = equatorial2galactic( ra_src, dec_src)
ll_src = ll_src/3.1415926*180
bb_src = bb_src/3.1415926*180

print("For LE src file name: ", tmpstr[0:tmppos], " R.A. for source: ", ra_src," Decl. for source: ", dec_src, " ll: ",ll_src, " bb: ",bb_src)

spec_gxb = np.zeros(192)
spec_gxberr = np.zeros(192)
in1=0
in2=192
in3=192
in4=192*2
in5=192*2
in6=192*3

is_gxbflag = 0

if (np.abs(bb_src) <= 10):
    mspec0 = np.zeros(192)
    mspec0err = np.zeros(192)
    mspec1 = np.zeros(192)
    mspec1err = np.zeros(192)
    mspec2 = np.zeros(192)
    mspec2err = np.zeros(192)

    mspec = np.zeros(192)
    mspecerr = np.zeros(192)


    expo_10x10_0=0
    expo_10x10_1=0
    expo_10x10_2=0
    grid_num=0
    for ii in xrange(0,11):
        for jj in xrange(0,11):
            if bb_src+5.5-jj >= -10 and bb_src+5.5-jj <= +10:
                src_index1 = int(ll_src+5-ii)+int((bb_src+10+5.5-jj))*360
                mspec0 = mspec0 + gxbkg_map[src_index1,in1:in2]
                mspec1 = mspec1 + gxbkg_map[src_index1,in3:in4]
                mspec2 = mspec2 + gxbkg_map[src_index1,in5:in6]
                tmspecerr_0 = np.sqrt(gxbkg_map[src_index1,in1:in2])
                mspec0err   = np.sqrt(mspec0err**2 + tmspecerr_0**2)
                tmspecerr_1 = np.sqrt(gxbkg_map[src_index1,in3:in4])
                mspec1err= np.sqrt(mspec1err**2 + tmspecerr_1**2)
                tmspecerr_2 = np.sqrt(gxbkg_map[src_index1,in5:in6])
                mspec2err= np.sqrt(mspec2err**2 + tmspecerr_2**2)
                expo_10x10_0 = expo_10x10_0 + gxbkg_expo0[src_index1]
                expo_10x10_1 = expo_10x10_1 + gxbkg_expo1[src_index1]
                expo_10x10_2 = expo_10x10_2 + gxbkg_expo2[src_index1]
                grid_num = grid_num + 1
    mspec0 = mspec0/expo_10x10_0
    mspec0err=mspec0err/expo_10x10_0
    mspec1 = mspec1/expo_10x10_1
    mspec1err=mspec1err/expo_10x10_1
    mspec2 = mspec2/expo_10x10_2
    mspec2err=mspec2err/expo_10x10_2

    mspec = mspec0 + mspec1 + mspec2
    mspecerr=np.sqrt(mspec0err**2 + mspec1err**2 + mspec2err**2)

    #for uu in xrange(0,192):
    #    print(uu,mspec0[uu],mspec1[uu],mspec2[uu],expo_10x10_0,expo_10x10_1,expo_10x10_2)


    #plt.figure()
    #plt.plot(mspecerr)
    #plt.show()
    #print(expo_10x10)
    fac_gxball_part = np.sum(mspec[154:183])/np.sum(simbkgmod_map_192[154:183])
    mspec0=mspec
    mspec = mspec-simbkgmod_map_192*fac_gxball_part
    #for uu in xrange(0,192):
    #    print(uu,mspec0[uu],mspec[uu],simbkgmod_map_192[uu],fac_gxball_part,np.sum(mspec0[154:183]),np.sum(simbkgmod_map_192[154:183]))
    #src_index = int(ll_src)*20+int((bb_src+10))
    #tmpspec = gxbkg_map[src_index,in1:in2] + gxbkg_map[src_index,in3:in4]+ gxbkg_map[src_index,in3:in4]
    spec_gxb = mspec
    spec_gxberr = mspecerr
    #spec_gxb[154:192] = 0
    #spec_gxb[120:192] = 0
    #spec_gxb[0:4] = 0
    #spec_gxb[0:192] = mspec/np.sum(mspec[14:88])*np.sum(tmpspec[14:88])
    #spec_gxberr = spec_gxberr/np.sum(mspec[14:88])*np.sum(tmpspec[14:88])
    lin = np.where(np.fabs(spec_gxb) >=1000 )
    if (np.size(lin)>0):
        spec_gxb[lin] = 0
        spec_gxb = spec_gxb*0
    is_gxbflag = 1
    #ee = spec_ch2*13/1536.0+0.1
    #popt1, pcov1 = curve_fit(spec_bld2, ee[28:160], spec_gxb[28:160], maxfev=int(9e8))
    #yfit2 = spec_bld2(ee,popt1[0],popt1[1])
    #gxb_smooth = yfit2
    #gxb_smooth[0:28] = spec_gxb[0:28]
    
    #plt.figure()
    #plt.plot(spec_gxb)
    #plt.plot(yfit2)
    #plt.show()

    #spec_gxb = gxb_smooth

'''
plt.figure()
plt.plot(cxbkg_map)
plt.plot(spec_gxb/8)
plt.show()
'''

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++Read the blind detector file+++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''Read DETID==10,28,46,21,53,85 event data (Blind detecter)'''

evt_list  = pf.open(evtfilename)
evt_tab   = evt_list[1].data
evt_time  = evt_tab.field('Time')
evt_detid = evt_tab.field('Det_ID')
evt_cha   = evt_tab.field('PI')
evt_type  = evt_tab.field('Event_Type')
evt_grade  = evt_tab.field('GRADE')
evt_list.close()

bld_spec_arr = np.zeros((trnum,ledetchans))
bld_detid_index = np.where(((evt_detid == 13)|(evt_detid == 45)|(evt_detid == 77)|(evt_detid == 21)|(evt_detid == 53)|(evt_detid == 85)) & (evt_type == 0) & (evt_grade == 0))
evt_time   = evt_time[bld_detid_index]
evt_detid  = evt_detid[bld_detid_index]
evt_cha    = evt_cha[bld_detid_index]
evt_type   = evt_type[bld_detid_index]
evt_grade  = evt_grade[bld_detid_index]


cha_ran = np.linspace(0,ledetchans,ledetchans+1)*chstep#@
channel = np.linspace(0,ledetchans-1,ledetchans)*chstep#@

bldspec = np.zeros((trnum,ledetchans+2))
bldexpo = np.zeros((trnum,1))
bkgspec = np.zeros((trnum,1536+2))
bkgspec_err = np.zeros((trnum,1536+2))
speccnter=0

'''Obtain the blind spectra for each 5x5 degrees'''

bin0=192*0
bin1=192*1
bin2=192*2
bin3=192*3
bin4=192*4
bin5=192*5
bin6=192*6
bin7=192*7
bin8=192*8
bin9=192*9

bin_cor0 = 146
bin_cor1 = 184

stin1 = 1168
stin2 = 1464

totalexpo_bld = 0

cnt_obs_g10keV_whole = 0
tmpexpo_whole = 0
cnter = 0

for jj in xrange(0,trnum):
    t0 = timbot[jj] + bldstart
    t1 = timrof[jj] + bldstart
    is_gti = is_ingti(START,STOP,t0,t1)
    if(is_gti==0):
        continue;
    ttindex = np.where((evt_time >= t0) & (evt_time < t1))
    tmpcha = evt_cha[ttindex]
    tmptime= evt_time[ttindex]
    tmpnum = np.size(tmptime)
    if (tmpnum<10):
        continue
    dtmpt  = tmptime[1:(tmpnum)] - tmptime[0:(tmpnum-1)]
    dtmptin = np.where(dtmpt>=10)
    tmpexpo  = np.max(tmptime) - np.min(tmptime) - np.sum(dtmpt[dtmptin])
    tmpspec,bins = np.histogram(tmpcha,bins=cha_ran,range=[0,ledetchans])
    #tmpspec1536,bins = np.histogram(tmpcha,bins=cha_ran,range=[0,ledetchans])#@
    totalexpo_bld = totalexpo_bld + (t1-t0+1)
    bldspec[jj,0] = (t0 + t1)/2.0
    bldspec[jj,1] = tmpexpo
    bldspec[jj,2:(ledetchans+2)] = tmpspec
    cnt_obs_g10keV = np.sum(tmpspec[bin_cor0:bin_cor1])
    if cnt_obs_g10keV<=0:
        cnt_obs_g10keV=1
    cnt_obs_g10keV_whole = cnt_obs_g10keV_whole+cnt_obs_g10keV
    tmpexpo_whole = tmpexpo_whole + tmpexpo
    tmpsimbkgmod_map = simbkgmod_map
    spec_simbkg_10keV=np.sum(tmpsimbkgmod_map[stin1:stin2])
    spec_simbkg_10keV_err = np.sqrt(np.sum(tmpsimbkgmod_map[stin1:stin2]*tmpsimbkgmod_map[stin1:stin2]))
    fac_fac = cnt_obs_g10keV*fac_all/(spec_simbkg_10keV*tmpexpo)
    fac_facerr=np.sqrt(cnt_obs_g10keV)*fac_all/(spec_simbkg_10keV*tmpexpo)
    cnter = cnter + 1
    #print('fac_fac=',fac_fac,'fac_facerr=',fac_facerr,'moderr1=',simbkgmod_maperr[1100:1120]*fac_fac*tmpexpo
    #print('cnt_obs_g10keV=',cnt_obs_g10keV,'spec_simbkg_10keV=',spec_simbkg_10keV,'tmpexpo=',tmpexpo,'trnum=',trnum,'fac_all=',fac_all,'cnter=',cnter)
    spec_gxb2 = np.interp(spec_ch,spec_ch2,spec_gxb)
    spec_gxb2[1050:1536]=0
    spec_gxb2err = np.interp(spec_ch,spec_ch2,spec_gxberr)
    #cxbkg_map2 = np.interp(spec_ch,spec_ch2,cxbkg_map)
    bkgspec[speccnter,0] = (t0 + t1)/2.0
    bkgspec[speccnter,1] = tmpexpo
    bkgspec[speccnter,2:(1536+2)] = simbkgmod_map*fac_fac*tmpexpo + spec_gxb2*tmpexpo/8*is_gxbflag + cxbkg_map*tmpexpo*(1-is_gxbflag)
    bkgspec[speccnter,2:23] = 0
    bkgspec[speccnter,1502:1538] = 0
    #plt.figure()
    #plt.plot(simbkgmod_map*fac_fac)
    #plt.plot(spec_gxb2*tmpexpo/8*is_gxbflag)
    #plt.plot(cxbkg_map*tmpexpo*(1-is_gxbflag))
    #plt.show()
    bkgspec_err[speccnter,0] = (t0 + t1)/2.0
    bkgspec_err[speccnter,1] = tmpexpo
    bkgspec_err[speccnter,2:(1536+2)] = np.sqrt( (simbkgmod_maperr*fac_fac*tmpexpo)**2 + (simbkgmod_map*fac_facerr*tmpexpo)**2 + (spec_gxb2err*tmpexpo/np.sqrt(8))**2)

    #plt.figure()
    #plt.plot(spec_gxb2*tmpexpo/1/8)
    #plt.plot(simbkgmod_map*fac_fac*tmpexpo)
    #plt.show()
    #print( "Time", bldspec[jj,0] - bldspec[0,0],'Exposure: ', t1-t0, ' Photon number: ',np.size(tmpcha),np.sum(bkgspec[speccnter,2:(ledetchans+2)] ))
    speccnter = speccnter+1

spec_simbkg_10keV=np.sum(tmpsimbkgmod_map[stin1:stin2])
spec_simbkg_10keV_err = np.sqrt(np.sum(tmpsimbkgmod_map[stin1:stin2]*tmpsimbkgmod_map[stin1:stin2]))

fac_whole = cnt_obs_g10keV_whole*fac_all/(spec_simbkg_10keV*tmpexpo_whole)
fac_wholeerr=np.sqrt(cnt_obs_g10keV_whole)*fac_all/(spec_simbkg_10keV*tmpexpo_whole)

#print('fac_whole',fac_whole,'cnt_obs_g10keV_whole',cnt_obs_g10keV_whole,'spec_simbkg_10keV',spec_simbkg_10keV,'expo',tmpexpo_whole,'fac_all',fac_all,'fac_wholeerr',fac_wholeerr)

spec_whole = simbkgmod_map*fac_whole*tmpexpo_whole
spec_err_whole =  np.sqrt( (simbkgmod_maperr*fac_whole*tmpexpo_whole)**2 + 0*(simbkgmod_map*fac_wholeerr*tmpexpo_whole)**2 + (spec_gxb2err*tmpexpo_whole/np.sqrt(8))**2 )



'''
tmpspec,bins = np.histogram(evt_cha,bins=cha_ran,range=[0,ledetchans])
x = channel
par0 = [np.max(tmpspec[96:192]),400,0]
bounds0=([0,300,0],[2*np.max(tmpspec[96:192]),600,1000])
popt2, pcov2 = curve_fit(spec_bld, x[75:192], tmpspec[75:192], p0=par0, bounds=bounds0, maxfev=int(9e8))
print(popt2)
popt1, pcov1 = curve_fit(spec_bld2, x[25:75], tmpspec[25:75], maxfev=int(9e8))

yfit = spec_bld(x,popt2[0],popt2[1],popt2[2])
yfit2 = spec_bld2(x,popt1[0],popt1[1])

plt.figure()
plt.plot(x,tmpspec)
plt.plot(x,yfit)
plt.plot(x,yfit2)
plt.show()
spec_bld_all = np.zeros(192,dtype='float')
spec_bld_all[0:192] = tmpspec
spec_bld_all[75:192] = yfit[75:192]
spec_bld_all[25:75] = yfit2[25:75]
'''
bkgspec_time = bkgspec[0:speccnter,0]
bkgspec_exp  = bkgspec[0:speccnter,1]
bkgspec = bkgspec[0:speccnter,0:1538]
trnum = speccnter

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++ Calculate time range for map ++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

start_in_map=np.zeros(GTI_num)
stop_in_map =np.zeros(GTI_num)
cnt_map = 0
for ii in xrange(0,GTI_num):
    tmpst1 = START[ii]
    tmpst2 = STOP[ii]
    tflag  = np.zeros(trnum)
    #print(ii,tmpst2-tmpst1)
    time_gtiflag(bkgspec_time,tmpst1,tmpst2,tflag)
    tmpin = np.where(tflag == 1)
    if(np.size(tmpin)==0):
        continue
    start_in_map[cnt_map] = np.min(tmpin)
    stop_in_map[cnt_map] = np.max(tmpin)
    cnt_map = cnt_map+1

start_in_map=start_in_map[0:cnt_map]
stop_in_map=stop_in_map[0:cnt_map]
print("Start:",start_in_map)
print("Stop:",stop_in_map)

GTI_num = cnt_map

time_bkgmod = bkgspec_time
rate_bkgmod = np.zeros(trnum)

spec_ch = np.linspace(0,1535,1536)
spec_ch2= np.linspace(0,ledetchans-1,ledetchans)*chstep + chstep*0.5 
spec_ch2= np.linspace(0,ledetchans-1,ledetchans)*chstep + 3.5 
tmpchmin = chmin
tmpchmax = chmax+1
'''
for ii in xrange(0,trnum):
    tmpexpo_arr = bkgspec[ii,1]
    tmpspec_arr = bkgspec[ii,2:1538]
    spec_cnt = np.interp(spec_ch,spec_ch2,tmpspec_arr)
    spec_cnt2= np.sum(spec_cnt[tmpchmin:tmpchmax])
    tmpexpo = np.sum(tmpexpo_arr)
    rate_bkgmod[ii] = spec_cnt2/tmpexpo
'''

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++++  Read specific GTI  +++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
slgti_flag = 0
tflag2 = np.zeros(trnum,dtype='int')
tindex2 = np.zeros(trnum,dtype='int')

if (slgti!=''):
    newgti = np.loadtxt(slgti)
    newgti_shape = newgti.shape
    newgti_num = newgti_shape[0]
    if(np.size(newgti)==2):
        newgti_start = newgti[0]
        newgti_stop  = newgti[1]
    if(np.size(newgti)>=4):
        newgti_start = newgti[0:newgti_num,0]
        newgti_stop  = newgti[0:newgti_num,1]

    tmptime = time_bkgmod
    time_gtiflag(tmptime,newgti_start,newgti_stop,tflag2)
    ehk_num2 = flag_selection(tflag2,tindex2)
    slgti_flag = 1



'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++++++++++++++++++++Calculate spectrum fro BLD'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

if sp_lc_select == 'spec':
    print("Calculate background spectra now.")
    tmpexpo_arr = bkgspec[0:trnum,1]
    tmpspec_arr = bkgspec[0:trnum,2:1538]
    tmpspec_err = bkgspec_err[0:trnum,2:1538]
    tmpexpo = 0
    spec_cnt2 = np.zeros(1536,dtype='float')
    spec_err2 = np.zeros(1536,dtype='float')
    if(slgti_flag==0):
        tmpexpo=np.sum(tmpexpo_arr)
        spec_cnt2= np.sum(tmpspec_arr,axis=0)
        spec_err2= np.sqrt(np.sum(tmpspec_err*tmpspec_err,axis=0))
    if(slgti_flag==1):
        for jj in xrange(0,trnum):
            if(tflag2[jj]==1):
                tmpexpo = tmpexpo + tmpexpo_arr[jj]
                spec_cnt2= spec_cnt2 + tmpspec_arr[jj,0:1536]
                spec_err2= np.sqrt(spec_err*spec_err+tmpspec_err[jj,0:1536]*tmpspec_err[jj,0:1536])
    spec_cnt2[0:7] = 0

    spec_cnt = spec_cnt2
    spec_err = spec_err2
    #spec_cnt = np.interp(spec_ch,spec_ch2,spec_smooth)
    print(slgti_flag)
    #plt.figure()
    #plt.plot(ee,spec_smooth,'C1')
    #plt.plot(ee,lowess_y,'C3')
    #plt.plot(spec_cnt2/tmpexpo,'C2')
    #plt.plot(ee,spec_bld_all,'C3')
    #plt.plot(spec_bldmod,'C2')
    #plt.plot(spec_cnt-spec_bldmod)
    #plt.show()
    tmpstr = src_name[0]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For LE src file name: ", tmpstr[0:tmppos])
    spec_list    = pf.open(tmpstr[0:tmppos])
    spec_tab     = spec_list[1].data
    spec_channel = spec_tab.field(0)
    spec_counts  = spec_tab.field(1)
    spec_hdr     = spec_list[1].header
    spec_list.close()
    #rr = (spec_counts/spec_hdr['exposure'])/(spec_cnt/tmpexpo)
    #rr0=np.mean(rr[200:256])
    #print(rr0)
    outname = outnam +'.pha'
    #write_bkgspec(outname,spec_ch,spec_cnt,spec_err,tmpexpo,spec_hdr)
    write_bkgspec(outname,spec_ch,spec_cnt,spec_err_whole,tmpexpo,spec_hdr)
    #plt.figure()
    #print("Expo=",tmpexpo)
    #plt.plot(rr)
    #plt.show()


'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++  Calculate light curve    +++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''


if sp_lc_select == 'lc':
    print("Calculate background spectra now.")
    src_name=[]
    print(sl_name)
    sf = open(sl_name)
    for sline in sf:
        print( sline )
        src_name.append(sline)
    sf.close()
    print(src_name)
    tmpstr = src_name[0]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For LE src file name: ", tmpstr[0:tmppos])
    lc_list    = pf.open(tmpstr[0:tmppos])
    lc_tab     = lc_list[1].data
    lc_time    = lc_tab.field(0)
    lc_counts  = lc_tab.field(1)
    lc_hdr     = lc_list[1].header
    lc_list.close()
    lc_num = np.size(lc_time)
    lc_bkg     = np.zeros(lc_num)
    lc_bkg_map = np.zeros(trnum)
    lc_bkg_err = np.zeros(trnum)
    lc_tim_map = bkgspec[0:trnum,0]
    if(chmin <0):
        print("Illigal minmum channel, which will be set to 0")
        chmin=0
    if(chmax > 1536):
        print("Illigal maximum channel, which will be set to 1535")
        chmmax=1535

    for jj in xrange(0,trnum):
        tmpchmin = chmin + 2
        tmpchmax = chmax + 2
        tmpexpo_arr = bkgspec[jj,1]
        tmpspec_arr = bkgspec[jj,tmpchmin:tmpchmax]
        tmpspec_err = bkgspec_err[jj,tmpchmin:tmpchmax]
        spec_cnt= np.sum(tmpspec_arr)
        spec_err= np.sqrt(np.sum(tmpspec_err*tmpspec_err))/np.sum(tmpexpo_arr)
        #print(spec_cnt/np.sum(tmpexpo_arr),tmpexpo_arr)
        lc_bkg_map[jj] = lc_bkg_map[jj] + spec_cnt/np.sum(tmpexpo_arr)
        lc_bkg_err[jj] = np.sqrt(lc_bkg_err[jj]*lc_bkg_err[jj] + spec_err*spec_err)
    #fbkg=interpolate.interp1d(lc_tim_map,lc_bkg_map,kind='cubic')
    #lc_bkg=fbkg(lc_time)
    #print(lc_bkg_map)
    norm_in = np.where(lc_bkg_map >= -1e4)
    lc_bkg = np.interp(lc_time,lc_tim_map[norm_in],lc_bkg_map[norm_in])
    lc_err = np.interp(lc_time,lc_tim_map[norm_in],lc_bkg_err[norm_in])
    #plt.figure()
    #plt.plot(lc_tim_map[norm_in]-lc_tim_map[0],lc_bkg_map[norm_in])
    #plt.show()
    outname = outnam + '.lc'
    print("The background light curve: ",outname)
    write_lcurve(outname,lc_time,lc_bkg,lc_err,lc_hdr)

'''
if sp_lc_select == 'lc':
    print("Calculate background spectra now.")
    src_name=[]
    sf = open(sl_name)
    for sline in sf:
        print( sline )
        src_name.append(sline)
    sf.close()
    tmpstr = src_name[0]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For LE src file name: ", tmpstr[0:tmppos])
    lc_list    = pf.open(tmpstr[0:tmppos])
    lc_tab     = lc_list[1].data
    lc_time    = lc_tab.field(0)
    lc_counts  = lc_tab.field(1)
    lc_hdr     = lc_list[1].header
    lc_list.close()
    lc_num = np.size(lc_time)
    lc_bkg     = np.zeros(lc_num)
    lc_bkg_all = np.zeros(lc_num)
    spec_ch = np.linspace(0,1535,1536)
    spec_ch2= np.linspace(0,ledetchans-1,ledetchans)*chstep + chstep*0.5 
    if(chmin <0):
        print("Illigal minmum channel, which will be set to 0")
        chmin=0
    if(chmax > 1536):
        print("Illigal maximum channel, which will be set to 255")
        chmax=1535
    for ii in xrange(0,GTI_num):
        tmpin1 = int(start_in_map[ii])
        tmpin2 = int(stop_in_map[ii])
        tmpin3 = int(stop_in_map[ii])
        if(tmpin1-tmpin2==0):
            tmpin2 = int( tmpin1+1 )
            tmpin3 = int( tmpin1+1 )
            if tmpin2 >= np.max(stop_in_map):
                tmpin2 = int( np.max(stop_in_map) )
        tmpst1 = time_bkgmod[tmpin1]-timestep*0.5
        tmpst2 = time_bkgmod[tmpin2]+timestep*0.5
        tmptime_bkgmod = time_bkgmod[tmpin1:tmpin3]
        tmprate_bkgmod = rate_bkgmod[tmpin1:tmpin3]
        tmpflag  =np.zeros(lc_num)
        time_gtiflag(lc_time,tmpst1,tmpst2,tmpflag)
        tmpin = np.where(tmpflag == 1)
        #print(tmpin)
        lc_bkg[tmpin] = np.interp(lc_time[tmpin],tmptime_bkgmod,tmprate_bkgmod)
    outname = outnam + '.lc'
    write_lcurve(outname,lc_time,lc_bkg,lc_hdr)
    #plt.figure()
    #plt.plot(lc_time-lc_time[0],lc_bkg)
    #plt.show()
'''

print("Finish.")




