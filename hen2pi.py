#!/usr/bin/env python
'''
Convert Energy to PI for ME & LE

Formula for LE: PI=1536*(E-0.1)/13.
Formula for ME: PI=1024*(E-3)/60.

'''
import sys
import os

if len(sys.argv)==2:
    help  = sys.argv[1]
    if(help=='-h'):
        print("hen2pi LE energy")
        print("hen2pi ME energy")
        sys.exit()

if len(sys.argv)==3:
    det_select = sys.argv[1]
    E          = float(sys.argv[2])

if det_select == 'LE':
    print("LE")
    if((E<0.1) | (E>13.1)):
        print("Energy should be input from 0.1 to 13.1")
        sys.exit()
    print("PI is ",int((E-0.1)/13.0*1536),  " for Energy: ", E)
    print("Finish")
if det_select == 'ME':
    print("ME")
    if((E<3.) | (E>63.)):
        print("Energy should be input from 3 to 63")
        sys.exit()
    print("PI is ",int((E-3)/60.0*1024),  " for Energy: ", E)
    print("Finish")

if (det_select != 'ME')&(det_select != 'LE'):
    print("Detector input error")





