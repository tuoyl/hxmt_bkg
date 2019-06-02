#!coding=utf-8
"""
#Version: 1.0                            Jun. 3rd, 2018
Version: 2.0                            Sep. 13th, 2018


USAGE:
    python newgti.py legradename oldgtiname newgtiname


DESCRIPTION:
    This task regenerates a new gtifile that can eliminate the 'peaks' in 
    HXMT-LE lightcure which can not be deleted automatically by 
    the HXMTDAS task legtigen.
    The data choosing method put forward by Jinyuan LIAO which comparing 
    the difference of count between LE big FoV and LE small FoV. Because 
    each CCD of these two FoV has same area, the source should creates 
    similar count in the pointing mode of HXMT, but the count of backgroud
    shows a positive correlation with the FoV. So we can use this 
    criteria to choose good data. This criterion can be presented by,

                   D-M(D) < Fsigma*error(D)

    where the value of [Fsigma] can be set in the main function [finggti] 
    by parameter (sigma).

    The first step is choosing the legrade data by the former version of 
    gtifile, if there are more than one GTI extension in the gtifile, 
    it will choose the intersection of all extensions.
    The second step is calculate the difference of 18 CCDs of big FoV and 
    60 CCDs of small FoV. Then, applying the above criteria.
    The final step is check the good time interval and delete the 
    short one (dt<100s, or can be set in the main function [findgti]
    by parameter (dtime), in units of second). Then the new gtifile can be 
    generated.

    This task also can be imported as a module in a PYTHON script. You need
    to put this file in a proper path which is included in your python 
    enviroment variables. And you can import the main function [findgti] in 
    your script.


PARAMETERS
    (sigma) [float]
        Please read the description. You can change this value in the end 
        of this script.

    (dtime) [float]
        Please read the description. You can change this value in the end 
        of this script.

Update:
factor = 5.35



QUESTIONS & BUGS:
    The author of this script is Yue ZHANG and MingYu GE, the idea comes 
    from JingYuan LIAO.
    Please report problems to zhangyue@ihep.ac.cn or gemy@ihep.ac.cn or
    liaojinyuan@ihep.ac.cn.
    And if you have any questions about this task, also contact with 
    these three emals.
"""
from __future__ import division
from astropy.io import fits
import numpy as np
from collections import OrderedDict
from sys import argv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

try:
    # Python 2
    xrange
except NameError:
    # Python 3, xrange is now named range
    xrange = range

print( "-"*10,"import end","-"*10 )
########################
def comp4(a1,a2,b1,b2):
    """两个区间交集,a1<a2; b1<b2"""
    if a2<b1 or b2<a1:#'空集'
        gtii = []
    else:
        lst1 = sorted([a1,a2,b1,b2])
        gtii = [lst1[1], lst1[2]]
    return gtii
def comp4lst(a1,a2,lstb1,lstb2):
    """一个区间与一个区间列表的交集,a1<a2; lstb1<lstb2"""
    gtilst = []
    for i1,f1 in enumerate(lstb1):
        gtii = comp4(a1,a2,f1,lstb2[i1])
        if len(gtii) != 0:
            gtilst.append(gtii)
    return gtilst
def comp4lstlst(lsta1,lsta2,lstb1,lstb2):
    """两个列表区间的交集,lsta1<lsta2; lstb1<lstb2"""
    lst1 = []; lst2 = []
    for i1,f1 in enumerate(lsta1):
        lstii = comp4lst(f1, lsta2[i1], lstb1, lstb2)
        for i2 in lstii:
            lst1.append(i2[0])
            lst2.append(i2[1])
    return lst1,lst2
def samegti1(fitsname, gtiornot):
    """单个fits文件中多个GTIextension的交集
    gti文件的gtiornot=1,event文件gtiornot=2"""
    if gtiornot == 1:
        keystart = 'start'
        keystop = 'stop'
    elif gtiornot == 2:
        keystart='tstart'
        keystop='tstop'
    fd = fits.open(fitsname)
    lst0 = fd[gtiornot].data.field(keystart)
    lst1 = fd[gtiornot].data.field(keystop)
    gtiexpendn = len(fd)-3
    for i1 in range(gtiexpendn):
        nn = i1+gtiornot
        startlsti = fd[nn].data.field(keystart)
        stoplsti = fd[nn].data.field(keystop)
        lst0, lst1 = comp4lstlst(lst0, lst1, startlsti, stoplsti)
    fd.close()
    return lst0, lst1
def Two_gti_1(fits1, fits2, type1, type2):
    """两个fits文件的gti交集"""
    gti11, gti12 = samegti1(fits1, type1)
    gti21, gti22 = samegti1(fits2, type2)
    gti1, gti2 = comp4lstlst(gti11, gti12, gti21, gti22)
    #[start0, start1, ..., startn]
    #[stop0, stop1, ..., stopn]
    return gti1, gti2
########################
def togti(boollst, tlst):
    """得到bool值列表和与其对应的时间列表的时间段"""
    boolcha = boollst[1:]^boollst[:-1] #寻找边界,TF翻转时为T
    boolcha[0] = boollst[0]^boolcha[0]
    boolcha[-1] = boollst[-1]^boolcha[-1]
    timeint = tlst[:-1][boolcha]
    time0 = timeint[::2]; time1 = timeint[1::2]
    #timestart, timestop
    return time0, time1
   #^^^^^^^END
def binnsecond(gtiStart, gtiStop, n):
    """若相邻好时间段间隔<=ns，则合并"""
    subtract = gtiStart[1:]-gtiStop[:-1]
    subtractB = subtract<=n
    gtiStartNew = []; gtiStopNew = []
    gtiStartNew.append(gtiStart[0])
    for i1 in range(len(subtractB)):
        if subtractB[i1]==False:
            gtiStartNew.append(gtiStart[i1+1])
            gtiStopNew.append(gtiStop[i1])
    gtiStopNew.append(gtiStop[-1])
#    ########
#    subtractNew = np.array(gtiStartNew[1:])-np.array(gtiStopNew[:-1])
#    subtract = np.concatenate([subtract, [999]])
#    subtractNew = np.concatenate([subtractNew, [999]])
#    for i1 in range(len(gtiStart)):
#        print gtiStart[i1], gtiStop[i1], subtract[i1]
#    print 'xxx'
#    for i1 in range(len(gtiStartNew)):
#        print gtiStartNew[i1], gtiStopNew[i1], subtractNew[i1]
    return np.array(gtiStartNew), np.array(gtiStopNew)
   #^^^^^^^END
def creatnewgti(gtiold, gtiout, tstart, tstop):
    #copy imformation from old gtifile
    old = fits.open(gtiold)
    oldheader0 = old[0].header
    oldheader1 = old[1].header
    oldheader2 = old[-1].header
    olddata2 = old[-1].data
    old.close()
    ######## 产生新fits
    ######## extension 0
    primarydata = fits.PrimaryHDU(header=oldheader0)
    ######## extension 1 2
    c11 = fits.Column(name='START', array=tstart, format='1D')
    c12 = fits.Column(name='STOP', array=tstop, format='1D')
    data1 = fits.BinTableHDU.from_columns([c11,c12])
    c21=fits.Column(name='a', array=[0,0],format='1D')
    c22=fits.Column(name='b', array=[0,0],format='1D')
    data2 = fits.BinTableHDU.from_columns([c21,c22])
    ########
    newfits = fits.HDUList([primarydata,data1,data2])
    newfits.writeto(gtiout,overwrite=True) #for server
#    newfits.writeto(gtiout,overwrite=True)
    ######## fits写入新数据
    naxis2 = int(len(tstart))
    oldheader1['NAXIS2'] = naxis2
    newfits = fits.open(gtiout,mode='update')
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
    print( 'done')
#    newfits.writeto(gtiout,overwrite=True)
   #^^^^^^^END
def plot_check(binlst0, ctBig0, ctSmall0,
        binlsta, ctBiga, ctSmalla,
        binlst0i, ctBig0i, ctSmall0i,
        jji, badtimeB, cha, cri, chamedian, dtime):
    
    badstarti, badstopi = togti(badtimeB, binlst0i)
    ########
    gd0, gd1 = togti(np.array(jji.values()), np.array(jji.keys()))
    gd00=[];gd11=[]
    for i1,f1 in enumerate(gd0):
        if gd1[i1]-f1>dtime:
            gd00.append(f1)
            gd11.append(gd1[i1])
    #^^^得到每次循环的好时间段
    ######## 每次初始值
    rangeB0 = np.zeros((binlst0i.shape))
    for i1,ta in enumerate(gd00):
        tb = gd11[i1]
        rangeBi = np.logical_and(binlst0i>ta,binlst0i<tb)
        allrangeB = np.logical_or(rangeB0,rangeBi)
        rangeB0 = allrangeB
    binlsta1 = binlst0i[allrangeB]
    ctBiga1 = ctBig0i[allrangeB]
    ctSmalla1 = ctSmall0i[allrangeB]
    ########################
    plt.figure(figsize=(12,8))
    ########
    plt.subplot(211)
    plt.plot(binlst0, ctSmall0, '.', label='small')#fits原始
    plt.plot(binlst0, ctBig0, '.', label='big')
    plt.plot(binlsta1, ctSmalla1,'.',label='small2')#每次初始值
    plt.plot(binlsta1, ctBiga1, '.', label='big2')
    for i2 in range(len(gd00)):#每次好时间段
        plt.axvspan(gd00[i2], gd11[i2], facecolor='#cccccc',
                alpha=0.5)
    plt.legend(loc=9)
    ########
    plt.subplot(212)
    plt.plot(binlst0i, cha, '.', label='cha')#本次差值
    plt.plot(binlsta, (ctBiga-ctSmalla), 's', markerfacecolor='none')#符合条件的
    plt.axhline(cri,color='red', label=str(cri))
    plt.axhline(-cri,color='red')
    for i1 in range(len(badstarti)):
        plt.axvspan(badstarti[i1], badstopi[i1], facecolor='#cccccc',
                alpha=0.5)
    plt.plot(binlst0i, badtimeB*100, color='black')
    plt.legend(loc=9)
    ########
    plt.show()
   #^^^^^^^END
def newgti(sigma, dtime, binlst0i, ctBig0i, ctSmall0i, jji,
        binlst0=0, ctBig0=0, ctSmall0=0):#forPlot
    cha = ctBig0i-ctSmall0i
    chamedian = np.median(cha)
    print('median==',chamedian)
    abscha = np.abs(cha-chamedian)
    cri = sigma*np.mean((ctBig0i+ctSmall0i)**0.5)
    badtimeB = abscha>cri  #找出不符合条件的时间点
#    badtimeB = cha>cri     #找出不符合条件的时间点
    for ii,ti in enumerate(binlst0i):
        if badtimeB[ii]:jji[ti]=False#将不符合的点赋值false
    ######## #为下次做准备，符合条件的点继续判断
#    gdtimeB = (abscha<=cri) #使用gdtime产生好数据
    gdtimeB = (cha<=cri) #使用gdtime产生好数据
    binlsta = binlst0i[gdtimeB]
    ctBiga = ctBig0i[gdtimeB]
    ctSmalla = ctSmall0i[gdtimeB]
    #
    #plot_check(binlst0, ctBig0, ctSmall0, #fits原始值
    #        binlsta, ctBiga, ctSmalla,#每次好事例
    #        binlst0i,ctBig0i, ctSmall0i,#每次初始值
    #        jji, badtimeB, cha, cri, chamedian, dtime)
    del badtimeB
    del binlst0i,ctBig0i,ctSmall0i
    return binlsta, ctBiga, ctSmalla, chamedian, jji
   #^^^^^^^END

def write_baddec(fname,detid_bad,t1,t2,type_bad,status):
    spec_col1 = fits.Column(name='DetID', format='I', array=detid_bad)
    spec_col2 = fits.Column(name='TIMERANGE', format='20A', array=t1)
    spec_col3 = fits.Column(name='TIMERANGE2', format='20A', array=t2)
    spec_col4 = fits.Column(name='TYPE', format='20A', array=type_bad)
    spec_col5 = fits.Column(name='STATUS', format='B', array=status)
    cols = fits.ColDefs([spec_col1, spec_col2, spec_col3, spec_col4, spec_col5])
    hdr = fits.Header()
    hdr['EXTNAME']  = "detectorStatus"
    hdr['TELESCOP'] = 'HXMT'
    hdr['INSTRUME'] = 'ME'
    hdr['VERSION'] = '1.00'
    hdr['CREATOR'] = 'HXMTDAS'
    hdr['UNITNUM'] = '1728'
    hdr['TSTART'] = '0'
    hdr['DATE'] = '2018-01-01T00:00:00'
    hdu = fits.BinTableHDU.from_columns(cols,header=hdr)
    hdu.writeto(fname,overwrite=True)

def findgti(piname, gtiname, gtioutname, sigma=5, dtime=30):
    fpi = fits.open(piname)
    timelst = fpi[1].data.field('time')
#    Detid = fpi[1].data.field('det_id')
    Asicid0_5 = fpi[1].data.field('asic_id')
    Fpgaid0_8 = fpi[1].data.field('fpga_id')
    EvtPI = fpi[1].data.field('PI')
    fpi.close()
    #选择数据1,gti文件和pi文件所有gti扩展的交集
    #gti01, gti02 = Two_gti_1(piname, gtiname, 2, 1) #function
    #选择数据2,使用gti交集筛选PI文件数据
    fgti = fits.open(gtiname)
    gti01_st = fgti[1].data.field('START')
    gti01_sp = fgti[1].data.field('STOP')
    fgti.close()

    rangeB0 = np.zeros((timelst.shape))
    print(gti01_st,gti01_sp)
    rangeB0 = np.zeros((timelst.shape))
    for i1, ta in enumerate(gti01_st):
        ta = np.floor(ta)+1
        tb = np.floor(gti01_sp[i1])
        dt = tb-ta
        if dt<dtime:continue
        rangeBi = np.logical_and(timelst>=ta,timelst<=tb)
        allrangeB = np.logical_or(rangeB0,rangeBi)
        rangeB0 = allrangeB
    del rangeB0
    timelst = timelst[allrangeB]
#    Detid = Detid[allrangeB]
    Asicid0_5 = Asicid0_5[allrangeB]
    Fpgaid0_8 = Fpgaid0_8[allrangeB]
    EvtPI = EvtPI[allrangeB]
    print( 'chooseing the gradedata in the input gtifile')
    #^^^^^^^END
    #1得到大视场与小视场事例,2得到大小视场计数率
    boxid = Fpgaid0_8//3
    asic0_53 = Asicid0_5+6*Fpgaid0_8
    asic0_17 = asic0_53-boxid*18
    del boxid, asic0_53, Asicid0_5, Fpgaid0_8
    FovSmallOnlyB = (((asic0_17>=0)&(asic0_17<=5)) | (asic0_17==7) | ((asic0_17>=12)&(asic0_17<=17))) 
    FovBigOnlyB = ((asic0_17==8)|(asic0_17==9))
    print(len(FovSmallOnlyB),len(EvtPI >=100))
    FovSmallOnlyB = FovSmallOnlyB & (EvtPI >=100) & (EvtPI <=600)
    FovBigOnlyB = FovBigOnlyB & (EvtPI >=100) & (EvtPI <=600)
    timesmall = timelst[FovSmallOnlyB]
    timebig = timelst[FovBigOnlyB]
    del FovSmallOnlyB, FovBigOnlyB
    ########
    t0 = timelst[0]; te = timelst[-1]
    binlst0 = np.arange(t0, te, 1)
    ctSmall0 = np.histogram(timesmall, bins=binlst0)[0]
    ctBig0 = np.histogram(timebig, bins=binlst0)[0]
    gooddataB = np.logical_and(ctSmall0>0, ctBig0>0) #排除0时间段
    ctSmall0 = ctSmall0[gooddataB]
    ctBig0 = ctBig0[gooddataB]
    binlst0 = binlst0[:-1][gooddataB]

    #tmpfac = np.median(ctSmall0/ctBig0)
    #print( "factor===",tmpfac)
    #factor = tmpfac
    factor = 5.35


    ctBig0 = factor*ctBig0
    ctSmall0 = smooth(ctSmall0, 20)
    ctBig0 = smooth(ctBig0, 20)
    jj = OrderedDict(zip(np.arange(t0, te-1, 1), gooddataB)) #每个时间点的好坏
    del timelst, timesmall, timebig
    del gooddataB
    #^^^^^^^END
    ######## 筛选好数据
    #plt.plot(binlst0, ctSmall0, '.', label='small')
    #plt.plot(binlst0, ctBig0,'.', label='small')
    #plt.plot(binlst0, ctBig0-ctSmall0,'.', label='small')
    #plt.plot(binlst0, ctBig0/ctSmall0,'.', label='small')
    #plt.show()
    binlsti = binlst0; ctBigi = ctBig0; ctSmalli = ctSmall0
    Dmedian = 100; chamedian0 = 1e6
    loopi=1
    while Dmedian > 0.01:
        binlsti, ctBigi, ctSmalli, chamediani, jj = \
                newgti(sigma, dtime, binlsti,ctBigi, ctSmalli, jj, 
                        binlst0, ctBig0, ctSmall0)#forPlot
        Dmedian = chamedian0 - chamediani
        chamedian0 = chamediani
        print( 'loop'+str(loopi),'median=',chamediani,'Dmedian=',Dmedian)
        loopi+=1
#    plt.show()
    ######## 完成筛选，产生gti数据
    xx1 = np.array([i1 for i1 in jj.values()])
    xx2 = np.array([i1 for i1 in jj.keys()])
    gtistart, gtistop = togti(xx1, xx2)
    #gtistart, gtistop = togti(np.array(jj.values()), np.array(jj.keys()))
    gtistart, gtistop = binnsecond(gtistart, gtistop, 5)
    gtistart1=[]; gtistop1=[]
#    print(gtistop-gtistart)
    for i1,ta in enumerate(gtistart):
        tb = gtistop[i1]
        if tb-ta>dtime:
            gtistart1.append(round(ta))
            gtistop1.append(round(tb))
#    print 'create new gtifile'
    print(np.array(gtistop1)-np.array(gtistart1))
    creatnewgti(gtiname, gtioutname, gtistart1, gtistop1)
   #^^^^^^^END
def findgti_ver2(piname, gtiname, gtioutname, baddecfile, baddec_new, sigma=3, dtime=30,ratio_cut=1.5):
    fpi = fits.open(piname)
    timelst = fpi[1].data.field('time')
    Detid = fpi[1].data.field('det_id')
    Asicid0_5 = fpi[1].data.field('asic_id')
    Fpgaid0_8 = fpi[1].data.field('fpga_id')
    EvtPI = fpi[1].data.field('PI')
    fpi.close()
    #选择数据1,gti文件和pi文件所有gti扩展的交集
    #gti01, gti02 = Two_gti_1(piname, gtiname, 2, 1) #function
    #选择数据2,使用gti交集筛选PI文件数据
    fgti = fits.open(gtiname)
    gti01_st = fgti[1].data.field('START')
    gti01_sp = fgti[1].data.field('STOP')
    fgti.close()

    rangeB0 = np.zeros((timelst.shape))
    print(gti01_st,gti01_sp)
    rangeB0 = np.zeros((timelst.shape))
    for i1, ta in enumerate(gti01_st):
        ta = np.floor(ta)+1
        tb = np.floor(gti01_sp[i1])
        dt = tb-ta
        if dt<dtime:continue
        rangeBi = np.logical_and(timelst>=ta,timelst<=tb)
        allrangeB = np.logical_or(rangeB0,rangeBi)
        rangeB0 = allrangeB
    del rangeB0
    '''Select the events by OLD gti'''
    timelst = timelst[allrangeB]
    Detid = Detid[allrangeB]
    Asicid0_5 = Asicid0_5[allrangeB]
    Fpgaid0_8 = Fpgaid0_8[allrangeB]
    EvtPI = EvtPI[allrangeB]
    print( 'chooseing the gradedata in the input gtifile')
    #^^^^^^^END
    #1得到大视场与小视场事例,2得到大小视场计数率
    boxid = Fpgaid0_8//3
    asic0_53 = Asicid0_5+6*Fpgaid0_8
    asic0_17 = asic0_53-boxid*18
    del boxid, asic0_53, Asicid0_5, Fpgaid0_8
    FovSmallOnlyB = (((asic0_17>=0)&(asic0_17<=5)) | (asic0_17==7) | ((asic0_17>=12)&(asic0_17<=17))) 
    FovBigOnlyB = ((asic0_17==8)|(asic0_17==9))
    print(len(FovSmallOnlyB),len(EvtPI >=100))
    FovSmallOnlyB = FovSmallOnlyB & (EvtPI >=100) & (EvtPI <=600)
    FovBigOnlyB = FovBigOnlyB & (EvtPI >=100) & (EvtPI <=600)
    timesmall = timelst[FovSmallOnlyB]
    timebig = timelst[FovBigOnlyB]
    del FovSmallOnlyB, FovBigOnlyB
    ########
    t0 = timelst[0]; te = timelst[-1]
    binlst0 = np.arange(t0, te, 1)
    ctSmall0 = np.histogram(timesmall, bins=binlst0)[0]
    ctBig0 = np.histogram(timebig, bins=binlst0)[0]
    gooddataB = np.logical_and(ctSmall0>0, ctBig0>0) #排除0时间段
    ctSmall0 = ctSmall0[gooddataB]
    ctBig0 = ctBig0[gooddataB]
    binlst0 = binlst0[:-1][gooddataB]

    #tmpfac = np.median(ctSmall0/ctBig0)
    #print( "factor===",tmpfac)
    #factor = tmpfac
    factor = 5.35


    ctBig0 = factor*ctBig0
    ctSmall0 = smooth(ctSmall0, 20)
    ctBig0 = smooth(ctBig0, 20)
    jj = OrderedDict(zip(np.arange(t0, te-1, 1), gooddataB)) #每个时间点的好坏
    del timesmall, timebig
    del gooddataB
    #^^^^^^^END
    ######## 筛选好数据
    #plt.plot(binlst0, ctSmall0, '.', label='small')
    #plt.plot(binlst0, ctBig0,'.', label='small')
    #plt.plot(binlst0, ctBig0-ctSmall0,'.', label='small')
    #plt.plot(binlst0, ctBig0/ctSmall0,'.', label='small')
    #plt.show()
    binlsti = binlst0; ctBigi = ctBig0; ctSmalli = ctSmall0
    Dmedian = 100; chamedian0 = 1e6
    loopi=1
    while Dmedian > 0.01:
        binlsti, ctBigi, ctSmalli, chamediani, jj = \
                newgti(sigma, dtime, binlsti,ctBigi, ctSmalli, jj, 
                        binlst0, ctBig0, ctSmall0)#forPlot
        Dmedian = chamedian0 - chamediani
        chamedian0 = chamediani
        print( 'loop'+str(loopi),'median=',chamediani,'Dmedian=',Dmedian)
        loopi+=1
#    plt.show()
    ######## 完成筛选，产生gti数据
    xx1 = np.array([i1 for i1 in jj.values()])
    xx2 = np.array([i1 for i1 in jj.keys()])
    gtistart, gtistop = togti(xx1, xx2)
    #gtistart, gtistop = togti(np.array(jj.values()), np.array(jj.keys()))
    gtistart, gtistop = binnsecond(gtistart, gtistop, 5)
    gtistart1=[]; gtistop1=[]
#    print(gtistop-gtistart)
    for i1,ta in enumerate(gtistart):
        tb = gtistop[i1]
        if tb-ta>dtime:
            gtistart1.append(round(ta))
            gtistop1.append(round(tb))
#    print 'create new gtifile'
    print(np.array(gtistop1)-np.array(gtistart1))
    creatnewgti(gtiname, gtioutname, gtistart1, gtistop1)
    if(baddecfile!=''):
        rangeB0 = np.zeros((timelst.shape))
        for i1, ta in enumerate(gtistart1):
            tb = gtistop1[i1]
            dt = tb-ta
            rangeBi = np.logical_and(timelst>ta,timelst<tb)
            allrangeB = np.logical_or(rangeB0,rangeBi)
            rangeB0 = allrangeB
        del rangeB0
        '''Select the events by NEW gti'''
        timelst = timelst[allrangeB]
        Detid   = Detid[allrangeB]
        EvtPI   = EvtPI[allrangeB]
        '''Read the bad detector information'''
        badf = fits.open(baddecfile)
        detid_bad = badf[1].data.field('DetID')
        t1_bad = badf[1].data.field('TIMERANGE')
        t2_bad = badf[1].data.field('TIMERANGE2')
        typ_bad = badf[1].data.field('TYPE')
        sta_bad = badf[1].data.field('STATUS')
        badf.close()
        '''Spectrum calcultion'''
        medetchans=1024
        det_num   =1728
        spec_ch = np.linspace(0,medetchans-1,medetchans)
        cha_ran = np.linspace(0,medetchans,medetchans+1)
        did_ran = np.linspace(0,det_num,det_num+1)
        spec_arr = np.zeros((det_num,medetchans))
        spec_ratio=np.zeros((det_num,1))
        spec_cnts=np.zeros((det_num,1))
        spec_detid=np.zeros((det_num,1))
        spec_flag =np.zeros((det_num,1))
        spec_arr2d,xd,yd = np.histogram2d(Detid,EvtPI,bins=(did_ran,cha_ran))
        print(np.size(spec_arr2d),np.size(xd),np.size(yd))
        print(xd,yd)
        #Xd,Yd=np.meshgrid(yd[0:1024],xd[0:1728])
        #fig = plt.figure()
        #ax0 = Axes3D(fig)
        #ax0.plot_surface(Xd,Yd,spec_arr2d,rstride=1, cstride=1)
        #plt.figure()
        #plt.plot(np.sum(spec_arr2d,axis=0))
        #plt.show()
        #plt.figure()
        #plt.plot(np.sum(spec_arr2d,axis=1))
        #plt.show()

        for ii in xrange(0,det_num):
            spec_detid[ii] = ii
            is_bad = np.where(detid_bad==ii)
            #print(ii,is_bad,np.size(is_bad),10*32,11*32,28*32,29*32,46*32,47*32)
            if (np.size(is_bad)==1): continue
            if (ii>=10*32) & (ii < 11*32): continue
            if (ii>=28*32) & (ii < 29*32): continue
            if (ii>=46*32) & (ii < 47*32): continue
            spec_flag[ii] = 1
            tmpspec = spec_arr2d[ii,0:1024]
            #rangeBi = np.where(Detid==ii)
            #tmpspec,bins = np.histogram(EvtPI[rangeBi],bins=cha_ran,range=[0,medetchans])
            #print(ii, ' Total number: ',np.sum(tmpspec))
            #spec_arr[ii,0:medetchans] = tmpspec
            tmpcnt1 = float(np.sum(tmpspec[0:200]))
            tmpcnt2 = float(np.sum(tmpspec[200:400]))+1
            spec_ratio[ii] = tmpcnt1/tmpcnt2
            spec_cnts[ii]  = np.sum(tmpspec[200:400])
            #print(ii,spec_ratio[ii],tmpcnt1,tmpcnt2)
            #plt.figure()
            #plt.plot(tmpspec)
            #plt.plot(spec_arr2d[ii,0:1024])
            #plt.show()
        '''Search bad detectors'''
        #print(spec_ratio)
        #plt.figure()
        #plt.plot(spec_ratio)
        #plt.show()
        good_index  = np.where(spec_flag==1)
        spec_ratio2 = spec_ratio[good_index]
        spec_detid2 = spec_detid[good_index]
        spec_cnts2  = spec_cnts[good_index]
        mspec_ratio = np.mean(spec_ratio2[np.where(spec_ratio2<=10)])
        mspec_ratio_norm = spec_ratio2/mspec_ratio
        ratio_bins = np.linspace(0,10,100)
        tmprr,bins = np.histogram(mspec_ratio_norm,bins=ratio_bins)
        #plt.figure()
        #plt.plot(bins[0:99],tmprr)
        #plt.show()
        baddet_index2 = np.where(mspec_ratio_norm>ratio_cut)
        detid_bad_new  =  spec_detid2[baddet_index2]
        new_num1 = np.size(detid_bad)
        new_num2 = np.size(detid_bad_new)
        new_num = new_num1 + new_num2
        print(new_num1,new_num2,detid_bad_new)
        baddet_arr = np.zeros(new_num)
        baddet_arr[0:(new_num1)] = detid_bad
        if new_num2 > 0:
            baddet_arr[new_num1:(new_num2+new_num1)] = detid_bad_new
        baddet_arr=np.unique(baddet_arr)
        new_num = np.size(baddet_arr)
        t1_newbad = np.zeros(new_num,'|S20')
        t2_newbad = np.zeros(new_num,'|S20')
        typ_newbad = np.zeros(new_num,'|S20')
        sta_newbad = np.zeros(new_num,'B')
        for ii in xrange(0,new_num):
            t1_newbad[ii]= '0'
            t2_newbad[ii]= 'INDEF'
            typ_newbad[ii]= 'bad'
            sta_newbad[ii]= '0'
        write_baddec(baddec_new,baddet_arr,t1_newbad,t2_newbad,typ_bad,sta_bad)
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
   #^^^^^^^END
########################
if __name__ == '__main__':
    gradename = argv[1]
    gti0name = argv[2]
    gtioutname = argv[3]
    if(len(argv)==4):
        findgti(gradename, gti0name, gtioutname, sigma=1.5, dtime=30)
    if(len(argv)>=5):
        baddecfile = argv[4]
        baddec_new = argv[5]
        ratio_cut = 1.5
        if(len(argv)==7):
            ratio_cut = float(argv[6])
        findgti_ver2(gradename, gti0name, gtioutname, baddecfile,baddec_new, sigma=1.5, dtime=30, ratio_cut=1.5)
    
#    gradename = '/home/shlocal/zhangyueDATA/task_newgti/gtinew/ME/megrade_101.fits'
#    gti0name = '/home/shlocal/zhangyueDATA/task_newgti/gtinew/ME/megti_101.fits'
#    gtioutname='megtinew.fits'
