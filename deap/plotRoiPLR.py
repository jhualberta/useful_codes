#/usr/bin/python
import operator
import sys,os
import collections
import subprocess
from subprocess import *
from ROOT import *
from rat import *
import numpy as np
from numpy import array
#### roi 55 data ############
ff = TFile("NewROI_0.1NRAccUpper_SecondPaperLower_190821.root")
fsb = TFile("saveSideBandROI.root")
roi55 = fsb.Get("top55")
roi30 = fsb.Get("top30")
roi05 = fsb.Get("top05")

roi55.SetName("roi_top55_gr")
roi55.SetVarX("nSCBayes");
roi55.SetVarY("rprompt60Bayes");

roi30.SetName("roi_top30_gr")
roi30.SetVarX("nSCBayes");
roi30.SetVarY("rprompt60Bayes");

roi05.SetName("roi_top05_gr")
roi05.SetVarX("nSCBayes");
roi05.SetVarY("rprompt60Bayes");

gr_roi_low = ff.Get("roi_low")
hroi_upper = ff.Get("hUpper") # histogram!!

npoints_low = gr_roi_low.GetN()-1 ## 106-1, x from 95 to 199
### there are 212 entires, don't use this!!
npoints_upper = int(hroi_upper.GetEntries())-1


nsteps1 = 20 # steps for vertical lines
nsteps2 = 40
npoints = nsteps1 + npoints_low + nsteps2 + npoints_low
print "!!!", npoints_low, npoints_low, npoints

roiPLR = TCutG("roi_plr_gr", npoints);
roiPLR.SetVarX("nSCBayes");
roiPLR.SetVarY("rprompt60Bayes");

## 212
## x = 95 to 200

## TGraph is the lower bound
xmin = gr_roi_low.GetX()[0] ## xmin = 95
yleftLower = gr_roi_low.GetY()[0] ## to draw vertical line, the left one
yleftUpper = hroi_upper.GetBinContent(hroi_upper.FindBin(xmin))

### alternative, to use xmax = 199
xmax = gr_roi_low.GetX()[npoints_low - 1]
yrightLower = gr_roi_low.GetY()[npoints_low - 1]
yrightUpper = hroi_upper.GetBinContent(hroi_upper.FindBin(xmax))

### build the TCutG by counter-clockwise (clockwise will give area<0!!), first draw left vertical line
### first draw the top curve, from right to left 
zip_upper = []
xx = xmin
print xmin, xmax, yleftLower, yleftUpper, yrightLower, yrightUpper

for i in range( npoints_low ):
    yy = hroi_upper.GetBinContent(i+1)
    #print xx, yy
    zip_upper.append( (xx, yy) )
    ##skip the same point
    #if i == 0:
    #    continue
    roiPLR.SetPoint(npoints_low-1-i, xx, yy)
    xx += 1

### 2nd, draw left vertical line, y from top to bottom
leftVerticalLine_x = [round(xmin,0) for i in np.linspace(yleftLower, yleftUpper, nsteps1)]
leftVerticalLine_y = [round(yy,6) for yy in np.linspace(yleftUpper, yleftLower, nsteps1)]
for i in range(nsteps1):
    xx = leftVerticalLine_x[i]
    yy = leftVerticalLine_y[i]
    #if i == 0:
    #    yy = yleftLower
    #if i == nsteps1-1:
    #    yy = yleftUpper
    roiPLR.SetPoint(i+npoints_low, xx, yy)


### 3rd, draw the bottom curve, from left to right
print "!!!!!!!!!", npoints_low
for i in range( npoints_low ):
    #if i == 0:
    #    continue
    xx = gr_roi_low.GetX()[i]
    yy = gr_roi_low.GetY()[i]
    ## print xx, yy
    roiPLR.SetPoint( i+npoints_low+nsteps1, xx, round(yy,6) )

### last, draw right vertical line, from bottom to top
rightVerticalLine_x = [round(xmax,0) for i in np.linspace(yrightLower, yrightUpper, nsteps2)]
rightVerticalLine_y = [round(yy,6) for yy in np.linspace(yrightLower, yrightUpper, nsteps2)]

for i in range(nsteps2):
    xx = rightVerticalLine_x[i]
    yy = rightVerticalLine_y[i]
    #if i == 0:
    #    continue
    #    yy = yrightLower
    #if i == nsteps-1:
    #    yy = yrightUpper
    #print xx, yy
    roiPLR.SetPoint( i+nsteps1+npoints_low+npoints_low, xx, round(yy,6) )

#
#roiPLR.Print()
roiPLR.Draw()
raw_input()


plr_area = roiPLR.Area()
top55_area = roi55.Area()
top30_area = roi30.Area()
top05_area = roi05.Area()

print "areas (count counter-clockwise): plr, top55, top30, top05", plr_area, top55_area, top30_area, top05_area
fnew = TFile("save_roi_PLR.root","recreate")
fnew.cd()
roiPLR.Write()
roi55.Write()
roi30.Write()
roi05.Write()
