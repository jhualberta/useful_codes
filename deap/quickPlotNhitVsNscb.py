import sys
import glob
from rat import *
from ROOT import *
import numpy as np
from numpy import sqrt
import time
from time import sleep
import struct
#from run_numbers import runs
from datetime import datetime, time as datetime_time, timedelta
## u232 source
fcalF = TFile("save_calF_nhitNSCB.root")
fcalC = TFile("save_calC_nhitNSCB.root")
fvac = TFile("save_pureVac_nhitNSCB.root")

level = 7;

h2d_calF = fcalF.Get("h2d_level"+str(level))
h2d_calC = fcalC.Get("h2d_level"+str(level))
h2d_vac = fvac.Get("h2d_level"+str(level))

h2d_calF.GetXaxis().SetRangeUser(0,400)
h2d_calC.GetXaxis().SetRangeUser(0,400)
h2d_vac.GetXaxis().SetRangeUser(0,400)

h2d_calF.GetYaxis().SetLabelOffset(1.2)
h2d_calC.GetYaxis().SetLabelOffset(1.2)
h2d_vac.GetYaxis().SetLabelOffset(1.2)

h2d_calF.GetXaxis().SetLabelOffset(0.8)
h2d_calC.GetXaxis().SetLabelOffset(0.8)
h2d_vac.GetXaxis().SetLabelOffset(0.8)

x1 = 60; x2 = 400
y1 = (-1.132e-3*x1) + 0.75
y2 = (-1.132e-3*x2) + 0.75
greenline = TLine(x1, y1, x2, y2)

greenline.SetLineColor(kGreen+2)
greenline.SetLineWidth(3)

## Red line: for separating neck-like/LG like Cherenkov
## NHit/PE = (-1.132e-3*PE) + 0.471692:
x1r = 60; x2r = 400
y1r = (-1.132e-3*x1r) + 0.471692
y2r = (-1.132e-3*x2r) + 0.471692
redline = TLine(x1r, y1r, x2r, y2r)
redline.SetLineColor(kRed)
redline.SetLineWidth(3)

c = TCanvas("c","",1400,600)
c.Divide(3,1)
c.cd(1)
gPad.SetLogz()
h2d_calC.Draw("colz")
redline.Draw("same")
greenline.Draw("same")
c.cd(2)
gPad.SetLogz()
h2d_calF.Draw("colz")
redline.Draw("same")
greenline.Draw("same")

c.cd(3)
gPad.SetLogz()
h2d_vac.Draw("colz")
redline.Draw("same")
greenline.Draw("same")

raw_input()
