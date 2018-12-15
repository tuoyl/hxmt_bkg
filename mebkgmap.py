#!/usr/bin/python
'''
Model constructed by Background Group.
Lian Jinyuan, Zhangshu, Guo Chengcheng, Jin Jing, Zhangjuan, Zhang Shu, et al.
Mail liaojinyuan@ihep.ac.cn

'''
'''
This version was written by Ge Mingyu anf Guo Chengcheng
Mail gemy@ihep.ac.cn

Add temperature effect
coe=(c[0]*Temp_fpga0_0 + c[1]*(MJD-57982) + c[2])/100 + 1

Usage:

mebkgmap lc/spec blind_det.FITS ehkfile.fits gtifile.fits deadtime.fits lcname/specname chmin chmax outnam_prefix (spec_time_arnge)
    lc/spec: lc for background lightcurve and spec for background light curve 
    blind_det_events.FITS: should only include the events for blind detecters.
    ehkfile.fits: the EHK file for the observation
    gtifile.fits: the GTI file for ME
    deadtime.fits: the Dead Time for ME
    Tempfile: input temperature file of ME
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
print "************ PRINT: mebkgmap -h for usage   *************"
print "*********************************************************"
print "HXMT background for Insight-HXMT/ME, ver-",Ver

uage_method1 = 'Method 1: mebkgmap lc/spec blind_det.FITS ehkfile.fits gtifile.fits deadtime.fits tempname lcname/specname chmin chmax outnam_prefix (spec_time_arnge)'
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
    ehkname       = sys.argv[3]
    gtifile       = sys.argv[4]
    dtname        = sys.argv[5]
    tempname      = sys.argv[6]
    sl_name       = sys.argv[7]
    chmin         = int(sys.argv[8])
    chmax         = int(sys.argv[9])
    outnam        = sys.argv[10]
    if (len(sys.argv)==12):
        slgti     = sys.argv[11]
else:
    sp_lc_select= str(raw_input("Selection(spec/lc):"))
    evtfilename = str(raw_input("Blind detector file:"))
    ehkname     = str(raw_input("EHK file:"))
    gtifile     = str(raw_input("GTI file:"))
    dtname      = str(raw_input("Dead time correction file:"))
    tempname      = str(raw_input("Temperature file:"))
    sl_name     = str(raw_input("FileName inlude spectra or lightcurve name:"))
    chmin      = int(raw_input("Minimum channel:"))
    chmax      = int(raw_input("Maximum channel:"))
    outnam      = str(raw_input("The prefix of output file name:"))
    slgti       = str(raw_input("Specific time range file(NONE):"))

#HEADAS=os.getenv('HEADAS')
HEADAS=os.getenv('REFPATH')
REFPATH=HEADAS+'/'

medetchans=1024
mefpganumber = 9

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
    print 'Non zeors index',nonzero_in[0:(fnum)]
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
    hdr['BACKSCAL'] = 0.875

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



'''Read gti extesion'''

hdulist = pf.open(gtifile)
tb = hdulist[1].data
START = tb.field(0)
STOP = tb.field(1)
hdulist.close()
GTI_num=np.size(START)
print("GTI START=",START)
print("GTI STOP=",STOP)


'''Read EHK file and position devision by 5 x 5'''

ehklist = pf.open(ehkname)
ehk_tab = ehklist[1].data
ehk_time= ehk_tab.field(0)
ehk_LON = ehk_tab.field(8)
ehk_LAT = ehk_tab.field(9)
ehklist.close()

ehk_num = np.size(ehk_time)
tflag = np.zeros(ehk_num,dtype='int')
tindex = np.zeros(ehk_num,dtype='int')
time_gtiflag(ehk_time,START,STOP,tflag)
ehk_num = flag_selection(tflag,tindex)
ehk_time = ehk_time[tindex[0:(ehk_num)]]
ehk_LON = ehk_LON[tindex[0:(ehk_num)]]
ehk_LAT = ehk_LAT[tindex[0:(ehk_num)]]

bkgmap_flag = np.zeros(ehk_num,dtype='int')
bkgmap_time(ehk_LON,ehk_LAT,bkgmap_flag)
bkgmap_flag = np.array(bkgmap_flag)
bkgmap_num = bkgmap_num(bkgmap_flag)

bkgmap_start_index = np.zeros(bkgmap_num,dtype='int')
bkgmap_stop_index = np.zeros(bkgmap_num,dtype='int')
bkgmap_interval_index(bkgmap_flag,bkgmap_start_index,bkgmap_stop_index)
bkgmap_arr_num  = bkgmap_stop_index - bkgmap_start_index + 1
bkgmap_arr_expo = bkgmap_stop_index - bkgmap_start_index + 1

bkgmap_start_time = ehk_time[bkgmap_start_index]
bkgmap_stop_time  = ehk_time[bkgmap_stop_index]
bkgmap_flag_uniq  = bkgmap_flag[bkgmap_start_index]

bkgspec_bld_bkgmap = np.zeros((bkgmap_num,medetchans),dtype='float')
bkgspec_all_bkgmap = np.zeros((bkgmap_num,medetchans+3),dtype='float')

print(bkgmap_arr_num)

'''Read Dead-Time file'''
dtimelist = pf.open(dtname)
dtime_tab = dtimelist[1].data
dtime_time = dtime_tab.field(0)
dtime_num = np.size(dtime_time)
dtime_arr = np.zeros((dtime_num,mefpganumber))

for ii in xrange(0,mefpganumber):
    tmpdt = dtime_tab.field(ii+1)
    dtime_arr[0:(dtime_num),ii] = tmpdt
dtimelist.close()

dtime_num = np.size(dtime_time)
tflag = np.zeros(dtime_num,dtype='int')
tindex = np.zeros(dtime_num,dtype='int')
time_gtiflag(dtime_time,START,STOP,tflag)
dtime_num = flag_selection(tflag,tindex)
dtime_time=dtime_time[tindex[0:(dtime_num)]]
dtime_arr=dtime_arr[tindex[0:(dtime_num)],0:mefpganumber]

dtime_cor=np.zeros(dtime_num)

for jj in xrange(0,dtime_num):
    tmpdeadtime = dtime_arr[jj,0:9]
    dtime_cor[jj] = (tmpdeadtime[1]+tmpdeadtime[4]+tmpdeadtime[7])/3.0

dtf = 1-np.mean(dtime_cor)

bkgmap_dtc = np.zeros((mefpganumber,bkgmap_num))

'''Cal the dead time correction for every time interval'''
for ii in xrange(0,mefpganumber):
    for jj in xrange(0,bkgmap_num):
        t0 = bkgmap_start_time[jj]
        t1 = bkgmap_stop_time[jj]
        tflag = np.zeros(dtime_num,dtype='int')
        tindex = np.zeros(dtime_num,dtype='int')
        time_gtiflag(dtime_time,t0,t1,tflag)
        tmpdt_num = flag_selection(tflag,tindex)
        tmpdeadtime = dtime_arr[tindex,ii]
        bkgmap_dtc[ii,jj] = np.sum(tmpdeadtime)/np.size(tmpdeadtime)

print bkgmap_dtc


'''Read Temperature file'''
Templist = pf.open(tempname)
Temp_tab = Templist[2].data
Temp_time = Temp_tab.field(0)
Temp_fpga0_0 = Temp_tab.field(1)
Temp_num = np.size(Temp_time)

tflag2 = np.zeros(Temp_num,dtype='int')
tindex2 = np.zeros(Temp_num,dtype='int')
time_gtiflag(Temp_time,START,STOP,tflag2)
Temp_num = flag_selection(tflag2,tindex2)
Temp_time=Temp_time[tindex2[0:(Temp_num)]]
Temp_fpga0_0=Temp_fpga0_0[tindex2[0:(Temp_num)]]

dt_MJD = np.mean(Temp_time)/86400.0 +55927+66.184/86400-57982
Temp_mean = np.mean(Temp_fpga0_0)
coe_t = np.loadtxt(REFPATH+'fullband_para.txt')
coe_Temp = (coe_t[0]*Temp_mean + coe_t[1]*dt_MJD + coe_t[2])/100. + 1

print(coe_Temp)

'''Read DETID==10,28,46 event data (Blind detecter)'''

evt_list  = pf.open(evtfilename)
evt_tab   = evt_list[1].data
evt_time  = evt_tab.field('Time')
evt_detid = evt_tab.field('ASIC_ID')
evt_cha   = evt_tab.field('PI')
evt_type  = evt_tab.field('Event_Type')
evt_list.close()

bld_spec_arr = np.zeros((bkgmap_num,medetchans))
#detid17_index = np.where((evt_detid == 10)|(evt_detid == 28)|(evt_detid == 46) & (evt_type == 0))
#evt_time  = evt_time[detid17_index]
#evt_cha   = evt_cha[detid17_index]

cha_ran = np.linspace(0,medetchans,medetchans+1)#@
channel = np.linspace(0,medetchans-1,medetchans)#@

'''Obtain the blind spectra for each 5x5 degrees'''
totalexpo_bld = 0
for jj in xrange(0,bkgmap_num):
    t0 = bkgmap_start_time[jj]
    t1 = bkgmap_stop_time[jj]
    ttindex = np.where((evt_time >= t0) & (evt_time < t1+1))
    tmpcha = evt_cha[ttindex]
    print 'Exposure: ', t1-t0+1, ' Photon number: ',np.size(tmpcha)
    tmpspec,bins = np.histogram(tmpcha,bins=cha_ran,range=[0,medetchans])#@
    bld_spec_arr[jj,0:medetchans] = tmpspec
    totalexpo_bld = totalexpo_bld + (t1-t0+1)
    #plt.figure()
    #plt.plot(channel,tmpspec)
    #plt.show()

'''Read background model'''
srcmapname=REFPATH+'/ME_bkgmap.fits'
srclist = pf.open(srcmapname)
srcmap_tab = srclist[1].data
srcmap_IN  = srcmap_tab.field(0)
srcmap_LON = srcmap_tab.field(1)
srcmap_LAT = srcmap_tab.field(2)
srcmap_BKG = srcmap_tab.field(3)
###det type mask


'''Cal the spectra of blind detectors from map'''
totalexpo_map = 0
for ii in xrange(0,bkgmap_num):
    tmpindex = bkgmap_flag_uniq[ii]
    tmpspec  = srcmap_BKG[tmpindex,0:2048]
    tmpexpo = bkgmap_arr_expo[ii]
    tin1 = 0
    tin2 = 1024
    if(np.sum(tmpspec[tin1:tin2]) == -1024):
        bkgspec_bld_bkgmap[ii,0:1024] = tmpspec[tin1:tin2]*0
        totalexpo_map = totalexpo_map + 0
    if(np.sum(tmpspec[tin1:tin2]) >= 0):
        bkgspec_bld_bkgmap[ii,0:1024] = tmpspec[tin1:tin2]*tmpexpo
        totalexpo_map = totalexpo_map + tmpexpo

    #print(np.sum(tmpspec[tin1:tin2]))
    #plt.figure()
    #plt.plot(channel,tmpspec[tin1:tin2]*tmpexpo)
    #plt.show()    

rr = (np.sum(bld_spec_arr,axis=0)/totalexpo_bld/dtf)/(np.sum(bkgspec_bld_bkgmap,axis=0)/totalexpo_map)

#plt.figure()
#plt.plot(channel,np.sum(bld_spec_arr,axis=0))
#plt.plot(channel,np.sum(bkgspec_bld_bkgmap,axis=0))
#plt.plot(channel,rr)
#plt.show()

rr0 = np.mean(rr[100:700])*coe_Temp * 0.875
print("Blind detector spectrum correction ratio: %f", rr0, " Temp: ", coe_Temp, "Dead time factor", dtf)
'''Cal the spectra of all detectors from map'''
tin1 = 1024
tin2 = 2048
for ii in xrange(0,bkgmap_num):
    tmpindex = bkgmap_flag_uniq[ii]
    tmpspec  = srcmap_BKG[tmpindex,0:2048]
    bkgspec_all_bkgmap[ii,0] = 0
    bkgspec_all_bkgmap[ii,1] = (bkgmap_start_time[ii]+bkgmap_start_time[ii])/2.
    #print(np.sum(tmpspec[tin1:tin2]))
    #print("ii===",ii,tmpindex,"tmpspec==",np.sum(tmpspec[0:256]))
    
    if(np.sum(tmpspec[tin1:tin2]) == -1024):
        bkgspec_bld_bkgmap[ii,0:1024] = tmpspec[tin1:tin2]*0
        bkgspec_all_bkgmap[ii,2] = 0
    if(np.sum(tmpspec[tin1:tin2]) >= 0):
        bkgspec_all_bkgmap[ii,2] = bkgmap_arr_expo[ii]
        bkgspec_all_bkgmap[ii,3:1027] = tmpspec[tin1:tin2]*bkgmap_arr_expo[ii]




'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++ Calculate time range for map ++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''

start_in_map=np.zeros(GTI_num)
stop_in_map =np.zeros(GTI_num)
time_bkgmod = bkgspec_all_bkgmap[0:bkgmap_num,1]
cnt_map = 0
for ii in xrange(0,GTI_num):
    tmpst1 = START[ii]
    tmpst2 = STOP[ii]
    tflag  =np.zeros(bkgmap_num)
    time_gtiflag(time_bkgmod,tmpst1,tmpst2,tflag)
    tmpin = np.where(tflag == 1)
    if(np.size(tmpin)==0):
        continue
    start_in_map[ii] = np.min(tmpin)
    stop_in_map[ii] = np.max(tmpin)
    cnt_map = cnt_map+1

print("Start:",start_in_map)
print("Stop:",stop_in_map)

GTI_num = cnt_map

rate_bkgmod = np.zeros(bkgmap_num)
tmpchmin = chmin
tmpchmax = chmax+1

for ii in xrange(0,bkgmap_num):
    tmpexpo_arr = bkgspec_all_bkgmap[ii,2]
    tmpspec_arr = bkgspec_all_bkgmap[ii,3:1027]
    spec_cnt2= np.sum(tmpspec_arr[tmpchmin:tmpchmax])
    tmpexpo = np.sum(tmpexpo_arr)
    rate_bkgmod[ii] = spec_cnt2/tmpexpo

#plt.figure()
#plt.plot(time_bkgmod-time_bkgmod[0],rate_bkgmod)
#plt.show()

'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''+++++++++++++++++++++++++++++++++++++++Calculate spectrum fro BLD'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''
'''++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'''



if sp_lc_select == 'spec':

    src_name=[]
    print(sl_name)
    sf = open(sl_name)
    for sline in sf:
        print sline
        src_name.append(sline)
    sf.close()
    print(src_name)
    print("Calculate background spectra now.")
    spec_ch = np.linspace(0,1023,1024)

    tmpexpo_arr = bkgspec_all_bkgmap[0:bkgmap_num,2]
    tmpspec_arr = bkgspec_all_bkgmap[0:bkgmap_num,3:1027]
    spec_cnt= np.sum(tmpspec_arr,axis=0)
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
    print("For ME src file name: ", tmpstr[0:tmppos])
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
    #rr0 = 1.
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
    print(sl_name)
    sf = open(sl_name)
    for sline in sf:
        print sline
        src_name.append(sline)
    sf.close()
    print(src_name)
    tmpstr = src_name[0]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For ME src file name: ", tmpstr[0:tmppos])
    lc_list    = pf.open(tmpstr[0:tmppos])
    lc_tab     = lc_list[1].data
    lc_time    = lc_tab.field(0)
    lc_counts  = lc_tab.field(1)
    lc_hdr     = lc_list[1].header
    lc_list.close()
    lc_num = np.size(lc_time)
    lc_bkg     = np.zeros(lc_num)-1e6
    if(chmin <0):
        print("Illigal minmum channel, which will be set to 0")
        chmin=0
    if(chmax > 1023):
        print("Illigal maximum channel, which will be set to 1023")
        chmax=1023
    print("The total light curve: ")
    for ii in xrange(0,GTI_num):
        tmpin1 = int(start_in_map[ii])
        tmpin2 = int(stop_in_map[ii])
        if(tmpin1-tmpin2==0):
            tmpin2 = tmpin1+1
        tmpst1 = time_bkgmod[tmpin1]-00
        tmpst2 = time_bkgmod[tmpin2]+00
        tmptime_bkgmod = time_bkgmod[tmpin1:tmpin2]
        tmprate_bkgmod = rate_bkgmod[tmpin1:tmpin2]
        tmpflag  =np.zeros(lc_num)
        time_gtiflag(lc_time,tmpst1,tmpst2,tmpflag)
        tmpin = np.where(tmpflag == 1)
        lc_bkg[tmpin] = np.interp(lc_time[tmpin],tmptime_bkgmod,tmprate_bkgmod)
        #print(tmpin)
    tmpin1 = np.where(lc_bkg == -1e6)
    tmpin2 = np.where(lc_bkg > -1e4)
    lc_bkg[tmpin1] = np.interp(lc_time[tmpin1],lc_time[tmpin2],lc_bkg[tmpin2])
    outname = outnam + '.lc'
    write_lcurve(outname,lc_time,lc_bkg,lc_hdr)
    #plt.figure()
    #plt.plot(lc_time-lc_time[0],lc_bkg)
    #plt.show()


print("Finish.")
















