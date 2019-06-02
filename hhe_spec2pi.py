#!/usr/bin/env python3
#
#
#

from astropy.io import fits as pf
import numpy as np
import matplotlib.pyplot as plt
import os
import commands
import sys
import time
from scipy import interpolate
import datetime

uage_method1 = 'Method 1: hhe_spec2pi src.dat bkg.dat rsp.dat src.pi bkg.pi rsp.rsp'
uage_method2 = 'Method 2: Using interactive method in prompt.'
uage_method3 = 'Method 3: hhe_spec2pi srcphafile=src.dat backfile=bkg.dat respfile=rsp.dat srcoutname=src.pi bkgoutname=bkg.pi rspoutname=rsp.rsp '

def print_usage(uage_method1,uage_method2,uage_method3):
    print(uage_method1)
    print(uage_method2)
    print(uage_method3)

def check_argument():
    srcphafile = []
    backfile   = []
    respfile   = []
    srcoutname = []
    bkgoutname = []
    rspoutname = []
    len_arg = len(sys.argv)
    if len_arg <=1:
        raise IOError("Error input argument, RUN 'hhe_spec2pi -h' for help")
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print_usage(uage_method1,uage_method2)
        return False, False, False, False, False, False
    if((len_arg>1)&(len_arg<4)):
        raise IOError("Error input argument, RUN 'hhe_spec2pi -h' for help")
        sys.exit()
    sysflag = 0
    for i in xrange(len_arg):
        arg = sys.argv[i]
        if i == 0:continue
        arg_split = arg.split('=')
        if (len(arg_split)==1):continue
        argname = arg_split[0]
        argval  = arg_split[1]
        if argname == 'srcphafile': 
            srcphafile.append(argval)
            sysflag = sysflag + 1
        if argname == 'backfile': 
            backfile.append(argval)
            sysflag = sysflag + 1
        if argname == 'respfile':
            respfile.append( argval )
            sysflag = sysflag + 1
        if argname == 'srcoutname':
            srcoutname.append( argval )
            sysflag = sysflag + 1
        if argname == 'bkgoutname':
            bkgoutname.append( argval )
            sysflag = sysflag + 1
        if argname == 'rspoutname':
            rspoutname.append( argval )
            sysflag = sysflag + 1
    if((sysflag>0)&(sysflag<7)):
        raise IOError("Error input argument, RUN 'hhe_spec2pi -h' for help")
        sys.exit()
    if (sysflag==0):
        srcphafile.append(sys.argv[1])
        backfile.append(sys.argv[2])
        respfile.append(sys.argv[3])
        srcoutname.append(sys.argv[4])
        bkgoutname.append(sys.argv[5])
        rspoutname.append(sys.argv[6])

    return srcphafile[0],backfile[0],respfile[0],srcoutname[0],bkgoutname[0],rspoutname[0]

if len(sys.argv)==2:
    if sys.argv[1]=='-h':
        print_usage(uage_method1,uage_method2,uage_method3)
    sys.exit()
elif len(sys.argv)>=2:
    srcphafile,backfile,respfile,srcoutname,bkgoutname,rspoutname=check_argument()
else:
    srcphafile= str(raw_input("Source file:"))
    backfile  = str(raw_input("Background file:"))
    respfile  = str(raw_input("Response file:"))
    srcoutname= str(raw_input("Merged PI spectral name:"))
    bkgoutname  = str(raw_input("Merged PI background name:"))
    rspoutname  = str(raw_input("Merged response file:"))

def write_bkgspec(fname,channel,counts,spec_err,expo,hdr_ext):
    spec_qua= np.zeros(256)
    spec_grp= np.ones(256)
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
    hdr['POISSERR'] = False
    hdr['GROUPING'] = 0
    hdr['QUALITY']  = 0
    hdr['MJDREFI']  = hdr_ext['MJDREFI']
    hdr['MJDREFF']  = hdr_ext['MJDREFF']
    hdu = pf.BinTableHDU.from_columns(cols,header=hdr)
    hdu.writeto(fname,overwrite=True)


def publicheader0(h,nowdate):
    h.header['DATE'] = (nowdate,'Creation date')
    h.header['TELESCOP'] = ('HXMT','Telescope used')
    h.header['INSTRUME'] = ('HE','Instrument used')
    h.header['FILTER'] = ('NONE','the instrument filter in use (if any)')
    h.header['GRATING'] = ('NONE', 'grating used, if any')
    h.header['ORIGIN'] = ('HXMTrspgenerator','source of FITS file')
    h.header['REVISION'] = (1,'Revision')
    h.header['DETNAM'] = ('HE','Revision')

def publicheader1(h,n_chan):
    h.header['CHANTYPE'] = ('PI      ','(PHA or PI) uncorrected or corrected for gain')
    h.header['DETCHANS'] = (n_chan,'raw detector channels')
    h.header['HDUCLASS'] = ('OGIP    ','file format is OGIP standard')
    h.header['HDUCLAS1'] = ('RESPONSE','extension contains response data')
    h.header['HDUVERS'] = ('1.3.0   ','version of the file format')
    h.header['TLMIN4'] = (1,'first channel in the response')
    h.header['RMFVERSN'] = ('1992a   ','obsolete')
    h.header['HDUVERS1'] = ('1.1.0   ','obsolete')
    h.header['HDUVERS2'] = ('1.2.0   ','obsolete')
    h.header['CVSD0001'] = ('2018-09-12','UTC date of beginning applicability')
    h.header['CVST0001'] = ('00:00:00','UTC time of beginning applicability')
    h.header['CDES0001'] = ('HE Matrix','brief descriptive summary of this dataset')
    h.header['CCLS0001'] = ('BCF     ','brief descriptive summary of this dataset')
    h.header['LO_THRES'] = (0,'brief descriptive summary of this dataset')
    h.header['TLMIN4'] = ('0       ','brief descriptive summary of this dataset')
    h.header['CDTP0001'] = ('DATA    ','brief descriptive summary of this dataset')
    h.header['HDUCLAS3'] = ('COUNT   ','brief descriptive summary of this dataset')

def write_matrix(outname,resp_matrix):
    now = datetime.datetime.now()
    nowdate = now.date().strftime("%Y-%m-%d")

    ######## extension 0
    hdu0 = pf.PrimaryHDU()
    hdu0.header['LONGSTRN'] = ('OGIP 1.0','The HEASARC Long String Convention may be used')
    publicheader0(hdu0,nowdate=nowdate)

    ######## extension 1
    '''Create the first extension'''
    nrows = 485
    detchans = 256
    energ_lo = np.linspace(0,nrows-1,nrows) + 15.17
    energ_hi = np.linspace(0,nrows-1,nrows) + 16.17
    n_grp  = np.ones(nrows,dtype=int)
    f_chan = np.zeros(nrows,dtype=int)
    n_chan = np.ones(nrows,dtype=int)*256
    col1 = pf.Column(name='ENERG_LO',format='1E',unit='keV',array=energ_lo)
    col2 = pf.Column(name='ENERG_HI',format='1E',unit='keV',array=energ_hi)
    col3 = pf.Column(name='N_GRP',format='1I',array=n_grp)
    col4 = pf.Column(name='F_CHAN',format='1I',array=f_chan)
    col5 = pf.Column(name='N_CHAN',format='1I',array=n_chan)
    col6 = pf.Column(name='MATRIX',format='256E',array=resp_matrix)

    hdu1 = pf.BinTableHDU.from_columns([col1, col2, col3, col4, col5, col6])
    
    hdu1.header['EXTNAME'] = ('MATRIX','Extension name')
    hdu1.header['HDUCLAS2'] = ('RSP_MATRIX','extension contains a response matrix')
    publicheader0(hdu1,nowdate=nowdate)
    publicheader1(hdu1,n_chan=detchans)
    ######## extension 2
    '''Create the second extension'''
    Emin = 15.0
    Emax = 385.0
    en_step = (Emax-Emin)/detchans
    channel= np.linspace(0,detchans-1,detchans)
    e_min_2= np.linspace(0,detchans-1,detchans)*en_step + Emin
    e_max_2= (np.linspace(0,detchans-1,detchans)+1)*en_step + Emin

    col1 = pf.Column(name='CHANNEL', format='1I', array=channel)
    col2 = pf.Column(name='E_MIN', format='1E', unit='keV', array=e_min_2)
    col3 = pf.Column(name='E_MAX', format='1E', unit='keV', array=e_max_2)
    hdu2 = pf.BinTableHDU.from_columns([col1, col2, col3])
    hdu2.header['EXTNAME'] = ('EBOUNDS','Extension name')
    hdu2.header['HDUCLAS2'] = ('EBOUNDS','extension contains a response matrix')
    publicheader0(hdu2,nowdate=nowdate)
    publicheader1(hdu2,n_chan=detchans)
    '''Write ot file'''
    temp = []
    temp.append(hdu0)
    temp.append(hdu1)
    temp.append(hdu2)
    hdulists = pf.HDUList(temp)

    if( os.path.isfile(outname)):
        os.remove(outname)

    hdulists.writeto(outname)


HEADAS=os.getenv('REFPATH')
REFPATH=HEADAS+'/he_pha_to_pi_matrix/'

src_name=[]
bkg_name=[]
rsp_name=[]
sf = open(srcphafile)
for sline in sf:
    src_name.append(sline)
sf.close()
sf = open(backfile)
for sline in sf:
    bkg_name.append(sline)
sf.close()
sf = open(respfile)
for sline in sf:
    rsp_name.append(sline)
sf.close()

srcnum = np.size(src_name)

detnum = 18;
channelnum=256
matrix_ID = [];
for ii in xrange(0,detnum):
    tmpid = 'he_pha_to_pi_matrix_' + str(ii) +'.txt'
    matrix_ID.append(tmpid)

#print(matrix_ID)

spec_id   = np.zeros(srcnum)
spec_expo = np.zeros(srcnum)
back_expo = np.zeros(srcnum)

for ii in xrange(0,srcnum):
    tmpstr = src_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    spec_list    = pf.open(tmpstr[0:tmppos])
    spec_tab     = spec_list[1].data
    spec_hdr     = spec_list[1].header
    spec_expo[ii]= spec_hdr['EXPOSURE']
    spec_tab     = spec_list[3].data
    spec_id[ii]  = spec_tab.field(0)
    spec_list.close()


    tmpstr = bkg_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    back_list    = pf.open(tmpstr[0:tmppos])
    back_tab     = back_list[1].data
    back_hdr     = back_list[1].header
    back_expo[ii]= back_hdr['EXPOSURE']

print(spec_id)
print(spec_expo)
print(back_expo)

spec_mer = np.zeros([srcnum,channelnum])
back_mer = np.zeros([srcnum,channelnum])
back_mererr = np.zeros([srcnum,channelnum])

for ii in xrange(0,srcnum):
    '''Read the convert matrix'''
    tmpmatrix = matrix_ID[int(spec_id[ii])]
    con_matrix=np.loadtxt(REFPATH+tmpmatrix)
    print(REFPATH+tmpmatrix)
    
    tmpstr = src_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For detector:", ii, " src file name: ", tmpstr[0:tmppos])
    spec_list    = pf.open(tmpstr[0:tmppos])
    spec_tab     = spec_list[1].data
    spec_channel = spec_tab.field(0)
    spec_counts  = spec_tab.field(1)
    spec_hdr     = spec_list[1].header
    spec_list.close()
    tmppispec = np.zeros(channelnum)
    tmppispecerr = np.zeros(channelnum)
    for jj in xrange(0,channelnum):
        tmppispec[jj] = np.sum(spec_counts*con_matrix[0:channelnum,jj])
        #plt.figure()
        #plt.plot(con_matrix[0:channelnum,jj])
        #plt.show()
    spec_mer[ii,0:channelnum] = tmppispec
    
    tmpstr = bkg_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    back_list    = pf.open(tmpstr[0:tmppos])
    back_tab     = back_list[1].data
    back_channel = back_tab.field(0)
    back_counts  = back_tab.field(1)
    back_countserr= back_tab.field(2)
    back_hdr     = back_list[1].header
    back_list.close()
    for jj in xrange(0,channelnum):
        tmppispec[jj] = np.sum(back_counts*con_matrix[0:channelnum,jj])
        tmppispecerr[jj]=np.sum(back_countserr*con_matrix[0:channelnum,jj])
    back_mer[ii,0:channelnum] = tmppispec
    back_mererr[ii,0:channelnum] = tmppispecerr

weight = np.ones(srcnum)
spec_pi = np.zeros(channelnum)
back_pi = np.zeros(channelnum)
back_pierr = np.zeros(channelnum)

for ii in xrange(0,srcnum):
    weight[ii] = back_expo[ii]/spec_expo[ii]

for ii in xrange(0,srcnum):
    print("ii===",weight[ii])
    spec_pi = spec_pi + weight[ii]*spec_mer[ii,0:channelnum]
    back_pi = back_pi + back_mer[ii,0:channelnum]
    back_pierr = back_pierr + back_mererr[ii,0:channelnum]

'''
for ii in xrange(0,srcnum):
    weight[ii] = spec_expo[ii]/back_expo[ii]

for ii in xrange(0,srcnum):
    print("ii===",weight[ii])
    spec_pi = spec_pi + spec_mer[ii,0:channelnum]
    back_pi = back_pi + weight[ii]*back_mer[ii,0:channelnum]
'''

channel = np.linspace(0,channelnum-1,channelnum)
expo = np.mean(back_expo)
#expo = np.mean(spec_expo)
#plt.figure()
#plt.plot(spec_pi-back_pi)
#plt.show()

write_bkgspec(srcoutname,channel,spec_pi,np.sqrt(spec_pi),expo,spec_hdr)
write_bkgspec(bkgoutname,channel,back_pi,back_pierr,expo,spec_hdr)

'''Merge respfile'''
nrows = 485
detchans = 256
resp_matrix_pi = np.zeros([nrows,detchans])
resp_arf = np.zeros(nrows)

rmf_ID = [];
for ii in xrange(0,detnum):
    tmpid = 'hePI_NaI_detID_' + str(ii) +'.rmf'
    rmf_ID.append(tmpid)

for ii in xrange(0,srcnum):
    tmpmatrix = matrix_ID[int(spec_id[ii])]
    con_matrix=np.loadtxt(REFPATH+tmpmatrix)
    print(np.size(con_matrix))
    
    tmpstr = rsp_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For detector:", ii, " src file name: ", tmpstr[0:tmppos])
    resp_list    = pf.open(tmpstr[0:tmppos])
    resp_tab     = resp_list[1].data
    tresp_matrix  = resp_tab.field(5)
    #resp_hdr     = spec_list[1].header
    resp_list.close()
    tmppispec = np.zeros(channelnum)
    tresp_arf = np.zeros(nrows)
    tresp_matrix2 = tresp_matrix
    for ll in xrange(0,nrows):
        tmprsp = tresp_matrix[ll,0:detchans]
        tresp_arf[ll] = np.sum(tmprsp)
    #Read basic rmf
    rmf_list    = pf.open(REFPATH+rmf_ID[int(spec_id[ii])])
    rmf_tab     = rmf_list[1].data
    rmf_matrix  = rmf_tab.field(5)
    #resp_hdr     = spec_list[1].header
    rmf_list.close()

    for ll in xrange(0,nrows):
        resp_matrix_pi[ll,0:detchans] = resp_matrix_pi[ll,0:detchans] + rmf_matrix[ll,0:detchans]*tresp_arf[ll]


'''
for ii in xrange(0,srcnum):
    tmpmatrix = matrix_ID[int(spec_id[ii])]
    con_matrix=np.loadtxt(REFPATH+tmpmatrix)
    print(np.size(con_matrix))
    
    tmpstr = rsp_name[ii]
    tmppos = tmpstr.find('\n')
    if (len(tmpstr[0:tmppos]) == 0):
        print("Input file name error:")
        sys.exit()
    print("For detector:", ii, " src file name: ", tmpstr[0:tmppos])
    resp_list    = pf.open(tmpstr[0:tmppos])
    resp_tab     = resp_list[1].data
    tresp_matrix  = resp_tab.field(5)
    #resp_hdr     = spec_list[1].header
    resp_list.close()
    tmppispec = np.zeros(channelnum)
    tresp_arf = np.zeros(nrows)
    tresp_matrix2 = tresp_matrix
    for ll in xrange(0,nrows):
        tmprsp = tresp_matrix[ll,0:detchans]
        tresp_arf[ll] = np.sum(tmprsp)
        tresp_matrix2[ll,0:detchans] = tresp_matrix[ll,0:detchans]/np.sum(tmprsp)
    plt.figure()
    plt.plot(tresp_arf)
    plt.show
    for ll in xrange(0,nrows):
        tmprsp = tresp_matrix2[ll,0:detchans]
        for jj in xrange(0,channelnum):
            tmppispec[jj] = np.sum(tmprsp*con_matrix[0:channelnum,jj])
        tmppispec = tmppispec/np.sum(tmppispec)
        resp_matrix_pi[ll,0:detchans] = resp_matrix_pi[ll,0:detchans] + tmppispec*tresp_arf[ll]
'''
'''
for ii in xrange(0,485):
    plt.figure()
    plt.plot(resp_matrix_pi[ii,0:256])
    plt.show()
'''
write_matrix(rspoutname,resp_matrix_pi)


