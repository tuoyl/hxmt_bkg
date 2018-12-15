#!/usr/bin/python
'''
Model constructed by Background Group.
Lian Jinyuan, Zhangshu, Guo Chengcheng, Jin Jing, Zhangjuan, Zhang Shu, et al.
Mail liaojinyuan@ihep.ac.cn

'''
'''
This version was written by Ge Mingyu
Mail gemy@ihep.ac.cn

Usage:

lebkgmap lc/spec blind_det.FITS gtifile.fits lcname/specname chmin chmax outnam_prefix (spec_time_arnge)
    lc/spec: lc for background lightcurve and spec for background light curve 
    blind_det_events.FITS: should only include the events for blind detecters.
    gtifile.fits: the GTI file for ME
    lcname/specname: is ASCII file, which includes the name of the source file for small FOV
    chmin: minimum channel for light curve
    chmax: maximum channel for light curve
    outnam_prefix: the output prefix for the spectrum

Using interactive method in prompt.	

'''

from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import os
import commands
import sys
import time

Ver = '2.0.5'

print "*********************************************************"
print "******************  Running HXMT Bkg   ******************"
print "*********************************************************"
print "*********************************************************"
print "*********************************************************"
print "************ PRINT: lebkgmap -h for usage   *************"
print "*********************************************************"
print "HXMT background for Insight-HXMT/LE, ver-",Ver

uage_method1 = 'Method 1: lebkgmap lc/spec evt_bld.fits gtifile.fits lcname/specname chmin chmax outnam_prefix (spec_time_arnge)'
uage_method2 = 'Method 2: Using interactive method in prompt.'

def print_usage(uage_method1,uage_method2):
    print(uage_method1)
    print(uage_method2)

if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print_usage(uage_method1,uage_method2)
    sys.exit()
elif len(sys.argv)>=2:
    sp_lc_select  = sys.argv[1]
    evtfilename   = sys.argv[2]
    gtifile       = sys.argv[3]
    sl_name       = sys.argv[4]
    chmin         = int(sys.argv[5])
    chmax         = int(sys.argv[6])
    outnam        = sys.argv[7]
    if (len(sys.argv)==9):
        slgti     = sys.argv[8]
else:
    sp_lc_select= str(raw_input("Selection(spec/lc):"))
    evtfilename = str(raw_input("Blind detector file:"))
    gtifile     = str(raw_input("GTI file:"))
    sl_name     = str(raw_input("FileName inlude spectra or lightcurve name:"))
    chmin       = int(raw_input("Minimum channel:"))
    chmax       = int(raw_input("Maximum channel:"))
    outnam      = str(raw_input("The prefix of output file name:"))
    slgti       = str(raw_input("Specific time range file(NONE):"))

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
    print 'Non zeros index',nonzero_in[0:(fnum)]
    bkgmap_start_index[1:(inum)] = nonzero_in + np.ones(nonzero_num,dtype=np.int)
    bkgmap_start_index[0] = 0
    #print np.size(bkgmap_stop_index[0:(inum-2)]),np.size(nonzero_in)
    bkgmap_stop_index[0:(inum-1)] =nonzero_in + np.zeros(nonzero_num,dtype=np.int)
    bkgmap_stop_index[inum-1] = fnum-1



def write_bkgspec(fname,channel,counts,expo,hdr_ext):
    spec_qua= np.zeros(1024)
    spec_grp= np.ones(1024)
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
    lc_col2 = pf.Column(name='Counts', format='D', array=counts)
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

outspec='./lebldspec'
timestep=64.0
chstep = 8
ledetchans=192

'''For hespecgen '''
lespecgen='lespecgen '
evtfile_par=' evtfile='
outfile_par=' outfile='
userdetid_par=' userdetid="13 45 77"'
eventtype_par=' eventtype=1'
starttime_par=' starttime='
stoptime_par=' stoptime='
clobber_par=' clobber=yes'

'''Read gti extesion'''

hdulist = pf.open(gtifile)
tb = hdulist[1].data
START = tb.field(0)
STOP = tb.field(1)
hdulist.close()
print("GTI START=",START)
print("GTI STOP=",STOP)
GTI_num = np.size(START)

bldstart=np.min(START)
bldstop=np.max(STOP)

lcnum=int((bldstop-bldstart)/timestep)+1

timcen = np.linspace(0,lcnum-1,lcnum)*timestep + 0.5*timestep
timbot = np.linspace(0,lcnum-1,lcnum)*timestep
timrof = np.linspace(0,lcnum-1,lcnum)*timestep + timestep
trnum = np.size(timcen)

'''Read Coeffiecients'''

coe_list  = pf.open(REFPATH+'LE_coe.fits')
coe_tab   = coe_list[1].data
coe_A0  = coe_tab.field('A_00')
coe_A1  = coe_tab.field('A_01')
coe_A2  = coe_tab.field('A_02')
coe_B0  = coe_tab.field('B_00')
coe_B1  = coe_tab.field('B_01')
coe_B2  = coe_tab.field('B_02')
coe_list.close()

coe_A = coe_A0 + coe_A1 + coe_A2
coe_B = coe_B0 + coe_B1 + coe_B2

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++Read diffuse background for Galaxy plane+++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

src_name=[]
print(sl_name)
sf = open(sl_name)
for sline in sf:
    print sline
    src_name.append(sline)
sf.close()
print(src_name)
tmpstr = src_name[0]
tmppos = tmpstr.find('\n')
spec_list    = pf.open(tmpstr[0:tmppos])
spec_tab     = spec_list[1].data
spec_hdr     = spec_list[1].header
spec_list.close()
ra_src =spec_hdr['RA_OBJ'] 
dec_src =spec_hdr['DEC_OBJ'] 

gxbkg_list= pf.open(REFPATH+'LE_gxmap.fits')
gxbkg_tab = gxbkg_list[1].data
gxbkg_in  = gxbkg_tab.field(0)
gxbkg_ll  = gxbkg_tab.field(1)
gxbkg_bb  = gxbkg_tab.field(2)
gxbkg_map = gxbkg_tab.field(3)
gxbkg_list.close()

[ll_src,bb_src] = equatorial2galactic( ra_src, dec_src)
ll_src = ll_src/3.1415926*180
bb_src = bb_src/3.1415926*180

print("For LE src file name: ", tmpstr[0:tmppos], " R.A. for source: ", ra_src," Decl. for source: ", dec_src, " ll: ",ll_src, " bb: ",bb_src)

spec_gxb = np.zeros(192)

if (np.abs(bb_src) <= 10):
    src_index = int(ll_src/2)*10+int((bb_src+10)/2)
    tmpspec = gxbkg_map[src_index,0:192]
    spec_gxb[0:192] = tmpspec

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++Read the blind detector file+++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''Read DETID==10,28,46 event data (Blind detecter)'''

evt_list  = pf.open(evtfilename)
evt_tab   = evt_list[1].data
evt_time  = evt_tab.field('Time')
evt_detid = evt_tab.field('Det_ID')
evt_cha   = evt_tab.field('PI')
evt_type  = evt_tab.field('Event_Type')
evt_grade  = evt_tab.field('GRADE')
evt_list.close()

#bld_spec_arr = np.zeros((trnum,ledetchans))
#bld_detid_index = np.where((evt_detid == 13)|(evt_detid == 45)|(evt_detid == 77) & (evt_type == 0))
#evt_time  = evt_time[detid17_index]
#evt_cha   = evt_cha[detid17_index]

cha_ran = np.linspace(0,ledetchans,ledetchans+1)*chstep#@
channel = np.linspace(0,ledetchans-1,ledetchans)*chstep#@

bldspec = np.zeros((trnum,ledetchans+2))
bldexpo = np.zeros((trnum,1))
bkgspec = np.zeros((trnum,ledetchans+2))
speccnter=0

'''Obtain the blind spectra for each 5x5 degrees'''
totalexpo_bld = 0
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
    tmpspec,bins = np.histogram(tmpcha,bins=cha_ran,range=[0,ledetchans])#@
    totalexpo_bld = totalexpo_bld + (t1-t0+1)
    bldspec[jj,0] = timcen[jj] + bldstart
    bldspec[jj,1] = tmpexpo
    bldspec[jj,2:(ledetchans+2)] = tmpspec
    bkgspec[speccnter,0] = timcen[jj] + bldstart
    bkgspec[speccnter,1] = tmpexpo
    bkgspec[speccnter,2:(ledetchans+2)] = (coe_A*tmpexpo + tmpspec*coe_B)/chstep + spec_gxb*tmpexpo/chstep
    print( "Time", timcen[jj] + bldstart ,'Exposure: ', t1-t0, ' Photon number: ',np.size(tmpcha),np.sum(bkgspec[speccnter,2:(ledetchans+2)] ))
    speccnter = speccnter+1
    #plt.figure()
    #plt.plot(channel,tmpspec)
    #plt.show()


bkgspec_time = bkgspec[0:speccnter,0]
bkgspec_exp  = bkgspec[0:speccnter,1]
bkgspec = bkgspec[0:speccnter,0:194]
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
    print(ii,tmpst2-tmpst1)
    time_gtiflag(bkgspec_time,tmpst1,tmpst2,tflag)
    tmpin = np.where(tflag == 1)
    if(np.size(tmpin)==0):
        continue
    start_in_map[cnt_map] = np.min(tmpin)
    stop_in_map[cnt_map] = np.max(tmpin)
    cnt_map = cnt_map+1

print("Start:",start_in_map)
print("Stop:",stop_in_map)

GTI_num = cnt_map

time_bkgmod = bkgspec_time
rate_bkgmod = np.zeros(trnum)

spec_ch = np.linspace(0,1535,1536)
spec_ch2= np.linspace(0,ledetchans-1,ledetchans)*chstep + chstep*0.5 

tmpchmin = chmin
tmpchmax = chmax+1

for ii in xrange(0,trnum):
    tmpexpo_arr = bkgspec[ii,1]
    tmpspec_arr = bkgspec[ii,2:194]
    spec_cnt = np.interp(spec_ch,spec_ch2,tmpspec_arr)
    spec_cnt2= np.sum(spec_cnt[tmpchmin:tmpchmax])
    tmpexpo = np.sum(tmpexpo_arr)
    rate_bkgmod[ii] = spec_cnt2/tmpexpo

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++++++++++++++++++++Calculate spectrum fro BLD'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

if sp_lc_select == 'spec':

    print("Calculate background spectra now.")
    tmpexpo_arr = bkgspec[0:trnum,1]
    tmpspec_arr = bkgspec[0:trnum,2:194]
    spec_cnt2= np.sum(tmpspec_arr,axis=0)
    print(np.size(spec_ch2),np.size(spec_cnt2))
    spec_cnt = np.interp(spec_ch,spec_ch2,spec_cnt2)
    #plt.figure()
    #plt.plot(spec_cnt,'C1')
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
    tmpexpo=np.sum(tmpexpo_arr)
    #rr = (spec_counts/spec_hdr['exposure'])/(spec_cnt/tmpexpo)
    #rr0=np.mean(rr[200:256])
    #print(rr0)
    rr0 = 1.
    outname = outnam +'.pha'
    write_bkgspec(outname,spec_ch,spec_cnt*rr0,tmpexpo,spec_hdr)
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
    sf = open(sl_name)
    for sline in sf:
        print sline
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
        if(tmpin1-tmpin2==0):
            tmpin2 = tmpin1+1
        tmpst1 = time_bkgmod[tmpin1]-timestep*0.5
        tmpst2 = time_bkgmod[tmpin2]+timestep*0.5
        tmptime_bkgmod = time_bkgmod[tmpin1:tmpin2]
        tmprate_bkgmod = rate_bkgmod[tmpin1:tmpin2]
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


print("Finish.")





















