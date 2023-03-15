''' This comes with RAT, so source rat.env first! '''
import os
import couchdb
import sys
import glob
from rat import *
#from rat import PMTInfoUtil
from ROOT import *
import numpy as np
from numpy import *
import math
from math import sqrt
from datetime import datetime, time as datetime_time, timedelta
from array import array

##21604, 75803 seconds = 0.877349 days
##18831, 103383 seconds = 1.19657 days
##22387, 73903.6 seconds = 0.855366 days
##22394, 74506.8 seconds = 0.862347 days
### for 21604, LAr
lo = 15333.9-2*245.15;
hi = 15333.9+3*245.15;
### for 18831
lo_old = 15500
hi_old = 16500

ffroi = TFile("saveSideBandROI.root")
#### saveSideBandROI.root")
####roi_802_days_11March2020_nsc_rp60.root")
## from loose cut to tight cut
roicut_top05 = ffroi.Get("top05")
roicut_top05.SetVarX("nscb")

roicut_top30 = ffroi.Get("top30")
roicut_top30.SetVarX("nscb")

roicut_top55 = ffroi.Get("top55")
roicut_top55.SetVarX("nscb")


fLar   = TFile("Combined_Sina_Lar_u232_run21604.root");
f18831 = TFile("Combined_Sina_Lar_u232_run18831.root");
f22387 = TFile("Combined_Sina_Lar_u232_run22387.root");
f22394 = TFile("Combined_Sina_Lar_u232_run22394.root");

fileList = [fLar, f18831, f22387, f22394]

#### livetime sec
list_time = [75803, 103383, 73903.6, 74506.8]
list_ntp = []
ii = 0

ntupleTF2 = TNtuple("ntupleTF2","save TF2 info","nscb:rprompt60Bayes:mbX:mbY:mbZ:pulseGar")



for ff in fileList:
    ntp = ff.Get("ntupleTF2")
    list_ntp.append(ntp)
    ii = ii+1

##ntupleTF2 = TNtuple("ntupleTF2","save TF2 info","nscb:rprompt60Bayes:mbX:mbY:mbZ:pulseGar")
runs_countGamma_L10 = []
runs_countGamma_L11 = []
runs_countGamma_L12 = []
ListHist_L10 = []
ListHist_L11 = []
ListHist_L12 = []

fname = ["H_21604", "H_18831","H_22387","H_22394"]
ffnew = TFile("processed.root","recreate")
ffnew.cd()
ii = 0
for ntp in list_ntp:
   ### roi: top05 
   htemp0 = TH1F("htemp0", "", 25000, 0, 25000) ## nscb
   ntp.Project("htemp0", "nscb", "")
   htemp0.SetName(fname[ii]+"level10")
   htemp0.Write(fname[ii]+"level10")
   countGamma_L10 = htemp0.Integral(htemp0.FindBin(lo), htemp0.FindBin(hi))
   if fname[ii] == "H_18831":
       countGamma_L10 = htemp0.Integral(htemp0.FindBin(lo_old), htemp0.FindBin(hi_old))
   htemp1 = TH1F("htemp1", "", 25000, 0, 25000) ## nscb
   ntp.Project("htemp1", "nscb", "sqrt(mbX*mbX+mbY*mbY+mbZ*mbZ)<750")
   htemp1.SetName(fname[ii]+"level11")
   htemp1.Write(fname[ii]+"level11")
   countGamma_L11 = htemp1.Integral(htemp1.FindBin(lo), htemp1.FindBin(hi))
   if fname[ii] == "H_18831":
       countGamma_L11 = htemp1.Integral(htemp1.FindBin(lo_old), htemp1.FindBin(hi_old))
   htemp2 = TH1F("htemp2", "", 25000, 0, 25000)
   ntp.Project("htemp2", "nscb", "sqrt(mbX*mbX+mbY*mbY+mbZ*mbZ)<750 && pulseGar>2")
   htemp2.SetName(fname[ii]+"level12")
   htemp2.Write(fname[ii]+"level12")
   countGamma_L12 = htemp2.Integral(htemp2.FindBin(lo), htemp2.FindBin(hi))
   if fname[ii] == "H_18831":
       countGamma_L12 = htemp2.Integral(htemp2.FindBin(lo_old), htemp2.FindBin(hi_old))

   runs_countGamma_L10.append(countGamma_L10)
   runs_countGamma_L11.append(countGamma_L11) 
   runs_countGamma_L12.append(countGamma_L12)
   ## htemp2.Draw()
   ## raw_input()
   ListHist_L10.append(htemp0)
   ListHist_L10[ii].SetDirectory(0)
   ListHist_L11.append(htemp1)
   ListHist_L11[ii].SetDirectory(0)
   ListHist_L12.append(htemp2)
   ListHist_L12[ii].SetDirectory(0)
   htemp0.Reset()
   htemp1.Reset()
   htemp2.Reset()
   ii += 1

for ii in range(4):
   print "level 10", fname[ii], runs_countGamma_L10[ii]/list_time[ii]*1000, "mHz"
   print "level 11", fname[ii], runs_countGamma_L11[ii]/list_time[ii]*1000, "mHz"
   print "level 12", fname[ii], runs_countGamma_L12[ii]/list_time[ii]*1000, "mHz"

#ListHist_L11[0].Draw()
