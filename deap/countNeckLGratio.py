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
from array import array
#time1 = 48.08*3600*24 ## cal C neck runs!!!
choosefile = str(sys.argv[1])

flist = ["Combined_calC_neck_run28103to28487.root", "Combined_calF_neck_run28512to28849.root", "Combined_neckStudy_358runs_run30627to32149.root"]

timelist = [48.08*3600*24, 35.46*3600*24, 241.05*3600*24]

f = TFile(flist[0]); time1 = timelist[0]

if choosefile == 'calc':
    f = TFile(flist[0]); time1 = timelist[0]     
elif choosefile == 'calf': 
    f = TFile(flist[1]); time1 = timelist[1]
elif choosefile == 'vac':
    f = TFile(flist[2]); time1 = timelist[2]

RCUT = 720

froi = TFile("save_roi_PLR.root")
top05 = froi.Get("roi_top05_gr")
top05.SetVarX("nscb")
top05.SetVarY("rprompt")

top30 = froi.Get("roi_top30_gr")
top30.SetVarX("nscb")
top30.SetVarY("rprompt")

top55 = froi.Get("roi_top55_gr")
top55.SetVarX("nscb")
top55.SetVarY("rprompt")

plr = froi.Get("roi_plr_gr")
plr.SetVarX("nscb")
plr.SetVarY("rprompt")

sb05 = froi.Get("sideband05")
sb05.SetVarX("nscb")
sb05.SetVarY("rprompt")

sb30 = froi.Get("sideband30")
sb30.SetVarX("nscb")
sb30.SetVarY("rprompt")

sb55 = froi.Get("sideband55")
sb55.SetVarX("nscb")
sb55.SetVarY("rprompt")

#sb55.SetVarX("nscb")
zContF = TFile("tf2mb_nSCBayes_deltaz_contours.root","READ");
zCut = zContF.Get("cont90_cut;1");
rContF = TFile("tf2mb_nSCBayes_dist_after_deltaz90_contours.root","READ");
rCut = rContF.Get("cont85_cut;1");

level = 8 ## start from level 8 fmaxpe, 15-8+1 = 8

#cutVal_LG = 0.471692 ## cut value of neck-like/LG-like
#cutVal_neck = 0.75 ## cut value of two top neck-like populations!!

cutVal_LG = 0.43 ## cut value of neck-like/LG-like
cutVal_neck = 0.72 ## cut value of two top neck-like populations!!

## count[0]: fmaxpe
## count[1]: neckVetoN == 1
## count[2]: mbR
## count[3]
## count[4] 12 CFT2R + MBZ
## count[5] 13 CFB3R
## count[6] 14 TF2-MB(z)
## count[7] 15 TF2-MB(r)
### count neck-like
countNeckROI05 = [0 for i in range(level)]
countNeckROI30 = [0 for i in range(level)]
countNeckROI55 = [0 for i in range(level)]
countNeckPLR = [0 for i in range(level)]
countNeckSB05 = [0 for i in range(level)]
countNeckSB30 = [0 for i in range(level)]
countNeckSB55 = [0 for i in range(level)]

### two neck populations!! 
countNeck1ROI05 = [0 for i in range(level)]
countNeck1ROI30 = [0 for i in range(level)]
countNeck1ROI55 = [0 for i in range(level)]
countNeck1PLR = [0 for i in range(level)]
countNeck1SB05 = [0 for i in range(level)]
countNeck1SB30 = [0 for i in range(level)]
countNeck1SB55 = [0 for i in range(level)]

countNeck2ROI05 = [0 for i in range(level)]
countNeck2ROI30 = [0 for i in range(level)]
countNeck2ROI55 = [0 for i in range(level)]
countNeck2PLR = [0 for i in range(level)]
countNeck2SB05 = [0 for i in range(level)]
countNeck2SB30 = [0 for i in range(level)]
countNeck2SB55 = [0 for i in range(level)]

### count LG-like
countLgROI05 = [0 for i in range(level)]
countLgROI30 = [0 for i in range(level)]
countLgROI55 = [0 for i in range(level)]
countLgPLR = [0 for i in range(level)]
countLgSB05 = [0 for i in range(level)]
countLgSB30 = [0 for i in range(level)]
countLgSB55 = [0 for i in range(level)]

tree = f.Get("T")
entries = tree.GetEntries()
for i in range(entries):
  tree.GetEntry(i)
  neckVetoN = tree.neckVeto
  nscb = tree.nSCBayes
  rprompt = tree.rprompt60Bayes
  nhit = tree.nhit
  ### check neck/LG (red line)
  yy = nhit/nscb
  pulseG = tree.pulseG
  cft2r = tree.cft2r
  cfb3r = tree.cfb3r
  evtx = tree.evtx
  evty = tree.evty
  evtz = tree.evtz

  tf2x = tree.tf2x
  tf2y = tree.tf2y
  tf2z = tree.tf2z

  runID = tree.runID
  subrunID = tree.subrunID
  eventID = tree.eventID

  mbr = sqrt(evtx*evtx+evty*evty+evtz*evtz)
  promptPE = nscb*rprompt;
  mbPos = TVector3(evtx, evty, evtz);
  tf2Pos = TVector3(tf2x, tf2y, tf2z);

  if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
      if top05.IsInside(nscb, rprompt):countNeckROI05[0] += 1
      if top30.IsInside(nscb, rprompt):countNeckROI30[0] += 1
      if top55.IsInside(nscb, rprompt):countNeckROI55[0] += 1
      if plr.IsInside(nscb, rprompt):countNeckPLR[0] += 1
      if sb05.IsInside(nscb, rprompt): countNeckSB05[0] += 1
      if sb30.IsInside(nscb, rprompt): countNeckSB30[0] += 1
      if sb55.IsInside(nscb, rprompt): countNeckSB55[0] += 1
      ## neck population 1 (green line, above) 
      if yy>=(-1.132e-3*nscb) + cutVal_neck:
          if top05.IsInside(nscb, rprompt):countNeck1ROI05[0] += 1
          if top30.IsInside(nscb, rprompt):countNeck1ROI30[0] += 1
          if top55.IsInside(nscb, rprompt):countNeck1ROI55[0] += 1
          if plr.IsInside(nscb, rprompt): countNeck1PLR[0] += 1
          if sb05.IsInside(nscb, rprompt): countNeck1SB05[0] += 1
          if sb30.IsInside(nscb, rprompt): countNeck1SB30[0] += 1
          if sb55.IsInside(nscb, rprompt): countNeck1SB55[0] += 1
      else: 
          if top05.IsInside(nscb, rprompt):countNeck2ROI05[0] += 1
          if top30.IsInside(nscb, rprompt):countNeck2ROI30[0] += 1
          if top55.IsInside(nscb, rprompt): countNeck2ROI55[0] += 1
          if plr.IsInside(nscb, rprompt):countNeck2PLR[0] += 1
          if sb05.IsInside(nscb, rprompt): countNeck2SB05[0] += 1
          if sb30.IsInside(nscb, rprompt): countNeck1SB30[0] += 1
          if sb55.IsInside(nscb, rprompt): countNeck1SB55[0] += 1
  else: ## LG-like
      ## neck population 1 (green line, above)
      if top05.IsInside(nscb, rprompt):countLgROI05[0] += 1
      if top30.IsInside(nscb, rprompt):countLgROI30[0] += 1
      if top55.IsInside(nscb, rprompt):countLgROI55[0] += 1
      if plr.IsInside(nscb, rprompt):countLgPLR[0] += 1
      if sb05.IsInside(nscb, rprompt): countLgSB05[0] += 1
      if sb30.IsInside(nscb, rprompt): countLgSB30[0] += 1 
      if sb55.IsInside(nscb, rprompt): countLgSB55[0] += 1
  ### cut level 9
  if neckVetoN == 0:
      if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
           if top05.IsInside(nscb, rprompt):countNeckROI05[1] += 1
           if top30.IsInside(nscb, rprompt):countNeckROI30[1] += 1
           if top55.IsInside(nscb, rprompt):countNeckROI55[1] += 1
           if plr.IsInside(nscb, rprompt):countNeckPLR[1] += 1
           if sb05.IsInside(nscb, rprompt): countNeckSB05[1] += 1
           if sb30.IsInside(nscb, rprompt): countNeckSB30[1] += 1
           if sb55.IsInside(nscb, rprompt): countNeckSB55[1] += 1
           ## neck population 1 (green line, above) 
           if yy>=(-1.132e-3*nscb) + cutVal_neck:
               if top05.IsInside(nscb, rprompt):countNeck1ROI05[1] += 1
               if top30.IsInside(nscb, rprompt):countNeck1ROI30[1] += 1
               if top55.IsInside(nscb, rprompt):countNeck1ROI55[1] += 1
               if plr.IsInside(nscb, rprompt):countNeck1PLR[1] += 1
               if sb05.IsInside(nscb, rprompt): countNeck1SB05[1] += 1
               if sb30.IsInside(nscb, rprompt): countNeck1SB30[1] += 1
               if sb55.IsInside(nscb, rprompt): countNeck1SB55[1] += 1
           else:
               if top05.IsInside(nscb, rprompt):countNeck2ROI05[1] += 1
               if top30.IsInside(nscb, rprompt):countNeck2ROI30[1] += 1
               if top55.IsInside(nscb, rprompt):countNeck2ROI55[1] += 1
               if plr.IsInside(nscb, rprompt):countNeck2PLR[1] += 1
               if sb05.IsInside(nscb, rprompt): countNeck2SB05[1] += 1
               if sb30.IsInside(nscb, rprompt): countNeck2SB30[1] += 1
               if sb55.IsInside(nscb, rprompt): countNeck2SB55[1] += 1
      else: ## LG-like
          ## neck population 1 (green line, above)
          if top05.IsInside(nscb, rprompt):countLgROI05[1] += 1
          if top30.IsInside(nscb, rprompt):countLgROI30[1] += 1
          if top55.IsInside(nscb, rprompt):countLgROI55[1] += 1
          if plr.IsInside(nscb, rprompt):countLgPLR[1] += 1
          if sb05.IsInside(nscb, rprompt): countLgSB05[1] += 1
          if sb30.IsInside(nscb, rprompt): countLgSB30[1] += 1
          if sb55.IsInside(nscb, rprompt): countLgSB55[1] += 1

      ### cut level 10 
      if mbr<RCUT: ### moved from 800 to 720 (PLR analysis)
          if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
               if top05.IsInside(nscb, rprompt):countNeckROI05[2] += 1
               if top30.IsInside(nscb, rprompt):countNeckROI30[2] += 1
               if top55.IsInside(nscb, rprompt):countNeckROI55[2] += 1
               if plr.IsInside(nscb, rprompt):countNeckPLR[2] += 1
               if sb05.IsInside(nscb, rprompt): countNeckSB05[2] += 1
               if sb30.IsInside(nscb, rprompt): countNeckSB30[2] += 1
               if sb55.IsInside(nscb, rprompt): countNeckSB55[2] += 1
               ## neck population 1 (green line, above) 
               if yy>=(-1.132e-3*nscb) + cutVal_neck:
                   if top05.IsInside(nscb, rprompt):countNeck1ROI05[2] += 1
                   if top30.IsInside(nscb, rprompt):countNeck1ROI30[2] += 1
                   if top55.IsInside(nscb, rprompt):countNeck1ROI55[2] += 1
                   if plr.IsInside(nscb, rprompt):countNeck1PLR[2] += 1
                   if sb05.IsInside(nscb, rprompt): countNeck1SB05[2] += 1
                   if sb30.IsInside(nscb, rprompt): countNeck1SB30[2] += 1
                   if sb55.IsInside(nscb, rprompt): countNeck1SB55[2] += 1
               else:
                   if top05.IsInside(nscb, rprompt):countNeck2ROI05[2] += 1
                   if top30.IsInside(nscb, rprompt):countNeck2ROI30[2] += 1
                   if top55.IsInside(nscb, rprompt):countNeck2ROI55[2] += 1
                   if plr.IsInside(nscb, rprompt):countNeck2PLR[2] += 1
                   if sb05.IsInside(nscb, rprompt): countNeck2SB05[2] += 1
                   if sb30.IsInside(nscb, rprompt): countNeck2SB30[2] += 1
                   if sb55.IsInside(nscb, rprompt): countNeck2SB55[2] += 1
          else: ## LG-like
              ## neck population 1 (green line, above)
              if top05.IsInside(nscb, rprompt):countLgROI05[2] += 1
              if top30.IsInside(nscb, rprompt):countLgROI30[2] += 1
              if top55.IsInside(nscb, rprompt):countLgROI55[2] += 1
              if plr.IsInside(nscb, rprompt):countLgPLR[2] += 1
              if sb05.IsInside(nscb, rprompt): countLgSB05[2] += 1
              if sb30.IsInside(nscb, rprompt): countLgSB30[2] += 1
              if sb55.IsInside(nscb, rprompt): countLgSB55[2] += 1
          ### cut level 11
          if pulseG>2:
              if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
                   if top05.IsInside(nscb, rprompt):countNeckROI05[3] += 1
                   if top30.IsInside(nscb, rprompt):countNeckROI30[3] += 1
                   if top55.IsInside(nscb, rprompt):countNeckROI55[3] += 1
                   if plr.IsInside(nscb, rprompt):countNeckPLR[3] += 1
                   if sb05.IsInside(nscb, rprompt): countNeckSB05[3] += 1
                   if sb30.IsInside(nscb, rprompt): countNeckSB30[3] += 1
                   if sb55.IsInside(nscb, rprompt): countNeckSB55[3] += 1
                   ## neck population 1 (green line, above) 
                   if yy>=(-1.132e-3*nscb) + cutVal_neck:
                       if top05.IsInside(nscb, rprompt):countNeck1ROI05[3] += 1
                       if top30.IsInside(nscb, rprompt):countNeck1ROI30[3] += 1
                       if top55.IsInside(nscb, rprompt):countNeck1ROI55[3] += 1
                       if plr.IsInside(nscb, rprompt):countNeck1PLR[3] += 1
                       if sb05.IsInside(nscb, rprompt): countNeck1SB05[3] += 1
                       if sb30.IsInside(nscb, rprompt): countNeck1SB30[3] += 1
                       if sb55.IsInside(nscb, rprompt): countNeck1SB55[3] += 1
                   else:
                       if top05.IsInside(nscb, rprompt):countNeck2ROI05[3] += 1
                       if top30.IsInside(nscb, rprompt):countNeck2ROI30[3] += 1
                       if top55.IsInside(nscb, rprompt):countNeck2ROI55[3] += 1
                       if plr.IsInside(nscb, rprompt):countNeck2PLR[3] += 1
                       if sb05.IsInside(nscb, rprompt): countNeck2SB05[3] += 1
                       if sb30.IsInside(nscb, rprompt): countNeck2SB30[3] += 1
                       if sb55.IsInside(nscb, rprompt): countNeck2SB55[3] += 1
              else: ## LG-like
                  ## neck population 1 (green line, above)
                  if top05.IsInside(nscb, rprompt):countLgROI05[3] += 1
                  if top30.IsInside(nscb, rprompt):countLgROI30[3] += 1
                  if top55.IsInside(nscb, rprompt):countLgROI55[3] += 1
                  if plr.IsInside(nscb, rprompt):countLgPLR[3] += 1
                  if sb05.IsInside(nscb, rprompt): countLgSB05[3] += 1
                  if sb30.IsInside(nscb, rprompt): countLgSB30[3] += 1
                  if sb55.IsInside(nscb, rprompt): countLgSB55[3] += 1
              #### cut level 12
              if cft2r<0.04 and evtz<550: 
                  if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
                       if top05.IsInside(nscb, rprompt):countNeckROI05[4] += 1
                       if top30.IsInside(nscb, rprompt):countNeckROI30[4] += 1
                       if top55.IsInside(nscb, rprompt):countNeckROI55[4] += 1 
                       if plr.IsInside(nscb, rprompt):countNeckPLR[4] += 1
                       if sb05.IsInside(nscb, rprompt): countNeckSB05[4] += 1
                       if sb30.IsInside(nscb, rprompt): countNeckSB30[4] += 1
                       if sb55.IsInside(nscb, rprompt): countNeckSB55[4] += 1
                       ## neck population 1 (green line, above) 
                       if yy>=(-1.132e-3*nscb) + cutVal_neck:
                           if top05.IsInside(nscb, rprompt):countNeck1ROI05[4] += 1
                           if top30.IsInside(nscb, rprompt):countNeck1ROI30[4] += 1
                           if top55.IsInside(nscb, rprompt):countNeck1ROI55[4] += 1
                           if plr.IsInside(nscb, rprompt):countNeck1PLR[4] += 1
                           if sb05.IsInside(nscb, rprompt): countNeck1SB05[4] += 1
                           if sb30.IsInside(nscb, rprompt): countNeck1SB30[4] += 1
                           if sb55.IsInside(nscb, rprompt): countNeck1SB55[4] += 1
                       else:
                           if top05.IsInside(nscb, rprompt):countNeck2ROI05[4] += 1
                           if top30.IsInside(nscb, rprompt):countNeck2ROI30[4] += 1
                           if top55.IsInside(nscb, rprompt):countNeck2ROI55[4] += 1
                           if plr.IsInside(nscb, rprompt):countNeck2PLR[4] += 1
                           if sb05.IsInside(nscb, rprompt): countNeck2SB05[4] += 1
                           if sb30.IsInside(nscb, rprompt): countNeck2SB30[4] += 1
                           if sb55.IsInside(nscb, rprompt): countNeck2SB55[4] += 1
                  else: ## LG-like
                      ## neck population 1 (green line, above)
                      if top05.IsInside(nscb, rprompt):countLgROI05[4] += 1
                      if top30.IsInside(nscb, rprompt):countLgROI30[4] += 1
                      if top55.IsInside(nscb, rprompt):countLgROI55[4] += 1
                      if plr.IsInside(nscb, rprompt):countLgPLR[4] += 1
                      if sb05.IsInside(nscb, rprompt): countLgSB05[4] += 1
                      if sb30.IsInside(nscb, rprompt): countLgSB30[4] += 1
                      if sb55.IsInside(nscb, rprompt): countLgSB55[4] += 1
                  
                  #### cut level 13
                  if cfb3r<0.1:
                      if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
                           if top05.IsInside(nscb, rprompt):countNeckROI05[5] += 1
                           if top30.IsInside(nscb, rprompt):countNeckROI30[5] += 1
                           if top55.IsInside(nscb, rprompt):countNeckROI55[5] += 1
                           if plr.IsInside(nscb, rprompt):countNeckPLR[5] += 1
                           if sb05.IsInside(nscb, rprompt): countNeckSB05[5] += 1
                           if sb30.IsInside(nscb, rprompt): countNeckSB30[5] += 1
                           if sb55.IsInside(nscb, rprompt): countNeckSB55[5] += 1
                           ## neck population 1 (green line, above) 
                           if yy>=(-1.132e-3*nscb) + cutVal_neck:
                               if top05.IsInside(nscb, rprompt):countNeck1ROI05[5] += 1
                               if top30.IsInside(nscb, rprompt):countNeck1ROI30[5] += 1
                               if top55.IsInside(nscb, rprompt):countNeck1ROI55[5] += 1
                               if plr.IsInside(nscb, rprompt):countNeck1PLR[5] += 1
                               if sb05.IsInside(nscb, rprompt): countNeck1SB05[5] += 1
                               if sb30.IsInside(nscb, rprompt): countNeck1SB30[5] += 1
                               if sb55.IsInside(nscb, rprompt): countNeck1SB55[5] += 1
                           else:
                               if top05.IsInside(nscb, rprompt):countNeck2ROI05[5] += 1
                               if top30.IsInside(nscb, rprompt):countNeck2ROI30[5] += 1
                               if top55.IsInside(nscb, rprompt):countNeck2ROI55[5] += 1
                               if plr.IsInside(nscb, rprompt):countNeck2PLR[5] += 1
                               if sb05.IsInside(nscb, rprompt): countNeck2SB05[5] += 1
                               if sb30.IsInside(nscb, rprompt): countNeck2SB30[5] += 1
                               if sb55.IsInside(nscb, rprompt): countNeck2SB55[5] += 1
                      else: ## LG-like
                          ## neck population 1 (green line, above)
                          if top05.IsInside(nscb, rprompt):countLgROI05[5] += 1
                          if top30.IsInside(nscb, rprompt):countLgROI30[5] += 1
                          if top55.IsInside(nscb, rprompt):countLgROI55[5] += 1
                          if plr.IsInside(nscb, rprompt):countLgPLR[5] += 1
                          if sb05.IsInside(nscb, rprompt): countLgSB05[5] += 1
                          if sb30.IsInside(nscb, rprompt): countLgSB30[5] += 1
                          if sb55.IsInside(nscb, rprompt): countLgSB55[5] += 1

                      ## cut level 14 
                      if (zCut.IsInside(promptPE,(tf2Pos.Z()-mbPos.Z()))): 
                          if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like 
                               if top05.IsInside(nscb, rprompt):countNeckROI05[6] += 1
                               if top30.IsInside(nscb, rprompt):countNeckROI30[6] += 1
                               if top55.IsInside(nscb, rprompt):countNeckROI55[6] += 1
                               if plr.IsInside(nscb, rprompt):
                                   print "!!!! level 14 in ROI", runID, subrunID, eventID
                                   countNeckPLR[6] += 1
                               if sb05.IsInside(nscb, rprompt): countNeckSB05[6] += 1
                               if sb30.IsInside(nscb, rprompt): countNeckSB30[6] += 1
                               if sb55.IsInside(nscb, rprompt): countNeckSB55[6] += 1
                               ## neck population 1 (green line, above) 
                               if yy>=(-1.132e-3*nscb) + cutVal_neck:
                                   if top05.IsInside(nscb, rprompt):countNeck1ROI05[6] += 1
                                   if top30.IsInside(nscb, rprompt):countNeck1ROI30[6] += 1
                                   if top55.IsInside(nscb, rprompt):countNeck1ROI55[6] += 1
                                   if plr.IsInside(nscb, rprompt):countNeck1PLR[6] += 1
                                   if sb05.IsInside(nscb, rprompt): countNeck1SB05[6] += 1
                                   if sb30.IsInside(nscb, rprompt): countNeck1SB30[6] += 1
                                   if sb55.IsInside(nscb, rprompt): countNeck1SB55[6] += 1 
                               else:
                                   if top05.IsInside(nscb, rprompt):countNeck2ROI05[6] += 1
                                   if top30.IsInside(nscb, rprompt):countNeck2ROI30[6] += 1
                                   if top55.IsInside(nscb, rprompt):countNeck2ROI55[6] += 1
                                   if plr.IsInside(nscb, rprompt):countNeck2PLR[6] += 1
                                   if sb05.IsInside(nscb, rprompt): countNeck2SB05[6] += 1
                                   if sb30.IsInside(nscb, rprompt): countNeck2SB30[6] += 1
                                   if sb55.IsInside(nscb, rprompt): countNeck2SB55[6] += 1

                          else: ## LG-like
                              ## neck population 1 (green line, above)
                              if top05.IsInside(nscb, rprompt):countLgROI05[6] += 1
                              if top30.IsInside(nscb, rprompt):countLgROI30[6] += 1
                              if top55.IsInside(nscb, rprompt):countLgROI55[6] += 1
                              if plr.IsInside(nscb, rprompt):countLgPLR[6] += 1
                              if sb05.IsInside(nscb, rprompt): countLgSB05[6] += 1
                              if sb30.IsInside(nscb, rprompt): countLgSB30[6] += 1
                              if sb55.IsInside(nscb, rprompt): countLgSB55[6] += 1

                          ## cut level 15
                          if (rCut.IsInside(promptPE,(tf2Pos-mbPos).Mag())):
                              if yy>=(-1.132e-3*nscb) + cutVal_LG: ## neck-like
                                   if top05.IsInside(nscb, rprompt):countNeckROI05[7] += 1
                                   if top30.IsInside(nscb, rprompt):countNeckROI30[7] += 1
                                   if top55.IsInside(nscb, rprompt):countNeckROI55[7] += 1
                                   if plr.IsInside(nscb, rprompt):countNeckPLR[7] += 1
                                   if sb05.IsInside(nscb, rprompt): countNeckSB05[7] += 1
                                   if sb30.IsInside(nscb, rprompt): countNeckSB30[7] += 1
                                   if sb55.IsInside(nscb, rprompt): countNeckSB55[7] += 1
                                   ## neck population 1 (green line, above)
                                   if yy>=(-1.132e-3*nscb) + cutVal_neck:
                                       if top05.IsInside(nscb, rprompt):countNeck1ROI05[7] += 1
                                       if top30.IsInside(nscb, rprompt):countNeck1ROI30[7] += 1
                                       if top55.IsInside(nscb, rprompt):countNeck1ROI30[7] += 1
                                       if plr.IsInside(nscb, rprompt):
                                           print "!!!! seriously level 15 in neck-like ROI", runID, subrunID, eventID
                                           countNeck1PLR[7] += 1
                                       if sb05.IsInside(nscb, rprompt): countNeck1SB05[7] += 1
                                       if sb30.IsInside(nscb, rprompt): countNeck1SB30[7] += 1    
                                       if sb55.IsInside(nscb, rprompt): countNeck1SB55[7] += 1
                                   else:
                                       if top05.IsInside(nscb, rprompt):countNeck2ROI05[7] += 1
                                       if top30.IsInside(nscb, rprompt):countNeck2ROI30[7] += 1
                                       if top55.IsInside(nscb, rprompt):countNeck2ROI55[7] += 1
                                       if plr.IsInside(nscb, rprompt):countNeck2PLR[7] += 1
                                       if sb05.IsInside(nscb, rprompt): countNeck2SB05[7] += 1
                                       if sb30.IsInside(nscb, rprompt): countNeck2SB30[7] += 1
                                       if sb55.IsInside(nscb, rprompt): countNeck2SB55[7] += 1

                              else: ## LG-like
                                  ## neck population 1 (green line, above)
                                  if top05.IsInside(nscb, rprompt):countLgROI05[7] += 1
                                  if top30.IsInside(nscb, rprompt):countLgROI30[7] += 1
                                  if top55.IsInside(nscb, rprompt):countLgROI55[7] += 1
                                  if plr.IsInside(nscb, rprompt):countLgPLR[7] += 1
                                  if sb05.IsInside(nscb, rprompt): countLgSB05[7] += 1
                                  if sb30.IsInside(nscb, rprompt): countLgSB30[7] += 1
                                  if sb55.IsInside(nscb, rprompt): countLgSB55[7] += 1


print "neck-roi, neck-sideband, LG-roi, LG-sb, neck-roi30, neck-roi05, neck-roi30, neck-roi05"
print "ROI 55"

print "PLR analysis"
for  i in range(level):
    print countNeckPLR[i]/time1*1000
print "-------"
for  i in range(level):
    print countLgPLR[i]/time1*1000


print "ROI 55"
for  i in range(level):
    #print countNeckSB55[i]/time1*1000, countNeckROI55[i]/time1*1000, countLgSB55[i]/time1*1000, countLgROI55[i]/time1*1000 
    #print "level ", i+8, countNeckROI55[i]/time1*1000, countNeckSB55[i]/time1*1000, countLgROI55[i]/time1*1000, countLgSB05[i]/time1*1000
    print countNeckROI55[i]/time1*1000
print "-------"
for  i in range(level):
    print countLgROI55[i]/time1*1000 



print "ROI 30"
for  i in range(level):
    #print countNeckSB30[i]/time1*1000, countNeckROI30[i]/time1*1000, countLgSB30[i]/time1*1000, countLgROI30[i]/time1*1000 
    print countNeckROI30[i]/time1*1000
print "-------"
for  i in range(level):
    print countLgROI30[i]/time1*1000

print "ROI 05"
for  i in range(level):
    #print countNeckSB05[i]/time1*1000, countNeckROI05[i]/time1*1000, countLgSB05[i]/time1*1000, countLgROI05[i]/time1*1000 
    print countNeckROI05[i]/time1*1000
print "-------"
for  i in range(level):
    print countLgROI05[i]/time1*1000 

