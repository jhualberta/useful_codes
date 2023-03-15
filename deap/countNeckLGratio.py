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
f = TFile("Combined_calC_neck_run28103to28487.root"); time1 = 48.08*3600*24;
#f = TFile("Combined_calF_neck_run28512to28849.root"); time1 = 35.46*3600*24 ## cal F equator runs!!!
#f = TFile("Combined_neckStudy_358runs_run30627to32149.root"); time1 = 241.05*3600*24 ## vacuum

froi = TFile("saveSideBandROI.root")
top05 = froi.Get("top05")
#top05.SetVarX("nscb")
top30 = froi.Get("top30")
#top30.SetVarX("nscb")
top55 = froi.Get("top55")
#top55.SetVarX("nscb")
sb = froi.Get("cutSideBand")
#sb.SetVarX("nscb")
zContF = TFile("tf2mb_nSCBayes_deltaz_contours.root","READ");
zCut = zContF.Get("cont90_cut;1");
rContF = TFile("tf2mb_nSCBayes_dist_after_deltaz90_contours.root","READ");
rCut = rContF.Get("cont85_cut;1");

level = 8 ## start from level 8 fmaxpe, 15-8+1 = 8
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
countNeckSideband = [0 for i in range(level)]
### two neck populations!! 
countNeck1ROI05 = [0 for i in range(level)]
countNeck1ROI30 = [0 for i in range(level)]
countNeck1ROI55 = [0 for i in range(level)]
countNeck1Sideband = [0 for i in range(level)]

countNeck2ROI05 = [0 for i in range(level)]
countNeck2ROI30 = [0 for i in range(level)]
countNeck2ROI55 = [0 for i in range(level)]
countNeck2Sideband = [0 for i in range(level)]
### count LG-like
countLgROI05 = [0 for i in range(level)]
countLgROI30 = [0 for i in range(level)]
countLgROI55 = [0 for i in range(level)]
countLgSideband = [0 for i in range(level)]

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

  if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
      if top05.IsInside(nscb, rprompt):countNeckROI05[0] += 1
      if top30.IsInside(nscb, rprompt):countNeckROI30[0] += 1
      if top55.IsInside(nscb, rprompt):countNeckROI55[0] += 1
      if sb.IsInside(nscb, rprompt): countNeckSideband[0] += 1
      ## neck population 1 (green line, above) 
      if yy>=(-1.132e-3*nscb) + 0.75:
          if top05.IsInside(nscb, rprompt):countNeck1ROI05[0] += 1
          if top30.IsInside(nscb, rprompt):countNeck1ROI30[0] += 1
          if top55.IsInside(nscb, rprompt):countNeck1ROI55[0] += 1
          if sb.IsInside(nscb, rprompt): countNeck1Sideband[0] += 1
      else: 
          if top05.IsInside(nscb, rprompt):countNeck2ROI05[0] += 1
          if top30.IsInside(nscb, rprompt):countNeck2ROI30[0] += 1
          if top55.IsInside(nscb, rprompt):countNeck2ROI55[0] += 1
          if sb.IsInside(nscb, rprompt): countNeck2Sideband[0] += 1
    
  else: ## LG-like
      ## neck population 1 (green line, above)
      if top05.IsInside(nscb, rprompt):countLgROI05[0] += 1
      if top30.IsInside(nscb, rprompt):countLgROI30[0] += 1
      if top55.IsInside(nscb, rprompt):countLgROI55[0] += 1
      if sb.IsInside(nscb, rprompt): countLgSideband[0] += 1

  ### cut level 9
  if neckVetoN == 0:
      if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
           if top05.IsInside(nscb, rprompt):countNeckROI05[1] += 1
           if top30.IsInside(nscb, rprompt):countNeckROI30[1] += 1
           if top55.IsInside(nscb, rprompt):countNeckROI55[1] += 1
           if sb.IsInside(nscb, rprompt): countNeckSideband[1] += 1
           ## neck population 1 (green line, above) 
           if yy>=(-1.132e-3*nscb) + 0.75:
               if top05.IsInside(nscb, rprompt):countNeck1ROI05[1] += 1
               if top30.IsInside(nscb, rprompt):countNeck1ROI30[1] += 1
               if top55.IsInside(nscb, rprompt):countNeck1ROI55[1] += 1
               if sb.IsInside(nscb, rprompt): countNeck1Sideband[1] += 1
           else:
               if top05.IsInside(nscb, rprompt):countNeck2ROI05[1] += 1
               if top30.IsInside(nscb, rprompt):countNeck2ROI30[1] += 1
               if top55.IsInside(nscb, rprompt):countNeck2ROI55[1] += 1
               if sb.IsInside(nscb, rprompt): countNeck2Sideband[1] += 1

      else: ## LG-like
          ## neck population 1 (green line, above)
          if top05.IsInside(nscb, rprompt):countLgROI05[1] += 1
          if top30.IsInside(nscb, rprompt):countLgROI30[1] += 1
          if top55.IsInside(nscb, rprompt):countLgROI55[1] += 1
          if sb.IsInside(nscb, rprompt): countLgSideband[1] += 1

      ### cut level 10 
      if mbr<800:
          if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
               if top05.IsInside(nscb, rprompt):countNeckROI05[2] += 1
               if top30.IsInside(nscb, rprompt):countNeckROI30[2] += 1
               if top55.IsInside(nscb, rprompt):countNeckROI55[2] += 1
               if sb.IsInside(nscb, rprompt): countNeckSideband[2] += 1
               ## neck population 1 (green line, above) 
               if yy>=(-1.132e-3*nscb) + 0.75:
                   if top05.IsInside(nscb, rprompt):countNeck1ROI05[2] += 1
                   if top30.IsInside(nscb, rprompt):countNeck1ROI30[2] += 1
                   if top55.IsInside(nscb, rprompt):countNeck1ROI55[2] += 1
                   if sb.IsInside(nscb, rprompt): countNeck1Sideband[2] += 1
               else:
                   if top05.IsInside(nscb, rprompt):countNeck2ROI05[2] += 1
                   if top30.IsInside(nscb, rprompt):countNeck2ROI30[2] += 1
                   if top55.IsInside(nscb, rprompt):countNeck2ROI55[2] += 1
                   if sb.IsInside(nscb, rprompt): countNeck2Sideband[2] += 1
  
          else: ## LG-like
              ## neck population 1 (green line, above)
              if top05.IsInside(nscb, rprompt):countLgROI05[2] += 1
              if top30.IsInside(nscb, rprompt):countLgROI30[2] += 1
              if top55.IsInside(nscb, rprompt):countLgROI55[2] += 1
              if sb.IsInside(nscb, rprompt): countLgSideband[2] += 1

          ### cut level 11
          if pulseG>2:
              if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
                   if top05.IsInside(nscb, rprompt):countNeckROI05[3] += 1
                   if top30.IsInside(nscb, rprompt):countNeckROI30[3] += 1
                   if top55.IsInside(nscb, rprompt):countNeckROI55[3] += 1
                   if sb.IsInside(nscb, rprompt): countNeckSideband[3] += 1
                   ## neck population 1 (green line, above) 
                   if yy>=(-1.132e-3*nscb) + 0.75:
                       if top05.IsInside(nscb, rprompt):countNeck1ROI05[3] += 1
                       if top30.IsInside(nscb, rprompt):countNeck1ROI30[3] += 1
                       if top55.IsInside(nscb, rprompt):countNeck1ROI55[3] += 1
                       if sb.IsInside(nscb, rprompt): countNeck1Sideband[3] += 1
                   else:
                       if top05.IsInside(nscb, rprompt):countNeck2ROI05[3] += 1
                       if top30.IsInside(nscb, rprompt):countNeck2ROI30[3] += 1
                       if top55.IsInside(nscb, rprompt):countNeck2ROI55[3] += 1
                       if sb.IsInside(nscb, rprompt): countNeck2Sideband[3] += 1

              else: ## LG-like
                  ## neck population 1 (green line, above)
                  if top05.IsInside(nscb, rprompt):countLgROI05[3] += 1
                  if top30.IsInside(nscb, rprompt):countLgROI30[3] += 1
                  if top55.IsInside(nscb, rprompt):countLgROI55[3] += 1
                  if sb.IsInside(nscb, rprompt): countLgSideband[3] += 1

              #### cut level 12
              if cft2r<0.04 and evtz<550: 
                  if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
                       if top05.IsInside(nscb, rprompt):countNeckROI05[4] += 1
                       if top30.IsInside(nscb, rprompt):countNeckROI30[4] += 1
                       if top55.IsInside(nscb, rprompt):countNeckROI55[4] += 1
                       if sb.IsInside(nscb, rprompt): countNeckSideband[4] += 1
                       ## neck population 1 (green line, above) 
                       if yy>=(-1.132e-3*nscb) + 0.75:
                           if top05.IsInside(nscb, rprompt):countNeck1ROI05[4] += 1
                           if top30.IsInside(nscb, rprompt):countNeck1ROI30[4] += 1
                           if top55.IsInside(nscb, rprompt):countNeck1ROI55[4] += 1
                           if sb.IsInside(nscb, rprompt): countNeck1Sideband[4] += 1
                       else:
                           if top05.IsInside(nscb, rprompt):countNeck2ROI05[4] += 1
                           if top30.IsInside(nscb, rprompt):countNeck2ROI30[4] += 1
                           if top55.IsInside(nscb, rprompt):countNeck2ROI55[4] += 1
                           if sb.IsInside(nscb, rprompt): countNeck2Sideband[4] += 1

                  else: ## LG-like
                      ## neck population 1 (green line, above)
                      if top05.IsInside(nscb, rprompt):countLgROI05[4] += 1
                      if top30.IsInside(nscb, rprompt):countLgROI30[4] += 1
                      if top55.IsInside(nscb, rprompt):countLgROI55[4] += 1
                      if sb.IsInside(nscb, rprompt): countLgSideband[4] += 1
                  
                  #### cut level 13
                  if cfb3r<0.1:
                      if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
                           if top05.IsInside(nscb, rprompt):countNeckROI05[5] += 1
                           if top30.IsInside(nscb, rprompt):countNeckROI30[5] += 1
                           if top55.IsInside(nscb, rprompt):countNeckROI55[5] += 1
                           if sb.IsInside(nscb, rprompt): countNeckSideband[5] += 1
                           ## neck population 1 (green line, above) 
                           if yy>=(-1.132e-3*nscb) + 0.75:
                               if top05.IsInside(nscb, rprompt):countNeck1ROI05[5] += 1
                               if top30.IsInside(nscb, rprompt):countNeck1ROI30[5] += 1
                               if top55.IsInside(nscb, rprompt):countNeck1ROI55[5] += 1
                               if sb.IsInside(nscb, rprompt): countNeck1Sideband[5] += 1
                           else:
                               if top05.IsInside(nscb, rprompt):countNeck2ROI05[5] += 1
                               if top30.IsInside(nscb, rprompt):countNeck2ROI30[5] += 1
                               if top55.IsInside(nscb, rprompt):countNeck2ROI55[5] += 1
                               if sb.IsInside(nscb, rprompt): countNeck2Sideband[5] += 1

                      else: ## LG-like
                          ## neck population 1 (green line, above)
                          if top05.IsInside(nscb, rprompt):countLgROI05[5] += 1
                          if top30.IsInside(nscb, rprompt):countLgROI30[5] += 1
                          if top55.IsInside(nscb, rprompt):countLgROI55[5] += 1
                          if sb.IsInside(nscb, rprompt): countLgSideband[5] += 1

                      ## cut level 14 
                      if (zCut.IsInside(promptPE,(tf2Pos.Z()-mbPos.Z()))): 
                          if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like 
                               if top05.IsInside(nscb, rprompt):countNeckROI05[6] += 1
                               if top30.IsInside(nscb, rprompt):countNeckROI30[6] += 1
                               if top55.IsInside(nscb, rprompt):
                                   print "!!!! level 14 in ROI", runID, subrunID, eventID
                                   countNeckROI55[6] += 1
                               if sb.IsInside(nscb, rprompt): countNeckSideband[6] += 1
                               ## neck population 1 (green line, above) 
                               if yy>=(-1.132e-3*nscb) + 0.75:
                                   if top05.IsInside(nscb, rprompt):countNeck1ROI05[6] += 1
                                   if top30.IsInside(nscb, rprompt):countNeck1ROI30[6] += 1
                                   if top55.IsInside(nscb, rprompt):countNeck1ROI55[6] += 1
                                   if sb.IsInside(nscb, rprompt): countNeck1Sideband[6] += 1
                               else:
                                   if top05.IsInside(nscb, rprompt):countNeck2ROI05[6] += 1
                                   if top30.IsInside(nscb, rprompt):countNeck2ROI30[6] += 1
                                   if top55.IsInside(nscb, rprompt):countNeck2ROI55[6] += 1
                                   if sb.IsInside(nscb, rprompt): countNeck2Sideband[6] += 1

                          else: ## LG-like
                              ## neck population 1 (green line, above)
                              if top05.IsInside(nscb, rprompt):countLgROI05[6] += 1
                              if top30.IsInside(nscb, rprompt):countLgROI30[6] += 1
                              if top55.IsInside(nscb, rprompt):countLgROI55[6] += 1
                              if sb.IsInside(nscb, rprompt): countLgSideband[6] += 1
                          
                          ## cut level 15
                          if (rCut.IsInside(promptPE,(tf2Pos-mbPos).Mag())):
                              if yy>=(-1.132e-3*nscb) + 0.471692: ## neck-like
                                   if top05.IsInside(nscb, rprompt):countNeckROI05[7] += 1
                                   if top30.IsInside(nscb, rprompt):countNeckROI30[7] += 1
                                   if top55.IsInside(nscb, rprompt):countNeckROI55[7] += 1
                                   if sb.IsInside(nscb, rprompt): countNeckSideband[7] += 1
                                   ## neck population 1 (green line, above)
                                   if yy>=(-1.132e-3*nscb) + 0.75:
                                       if top05.IsInside(nscb, rprompt):countNeck1ROI05[7] += 1
                                       if top30.IsInside(nscb, rprompt):countNeck1ROI30[7] += 1
                                       if top55.IsInside(nscb, rprompt):
                                           print "!!!! level 15 in ROI", runID, subrunID, eventID
                                           countNeck1ROI55[7] += 1
                                       if sb.IsInside(nscb, rprompt): countNeck1Sideband[7] += 1
                                   else:
                                       if top05.IsInside(nscb, rprompt):countNeck2ROI05[7] += 1
                                       if top30.IsInside(nscb, rprompt):countNeck2ROI30[7] += 1
                                       if top55.IsInside(nscb, rprompt):countNeck2ROI55[7] += 1
                                       if sb.IsInside(nscb, rprompt): countNeck2Sideband[7] += 1

                              else: ## LG-like
                                  ## neck population 1 (green line, above)
                                  if top05.IsInside(nscb, rprompt):countLgROI05[7] += 1
                                  if top30.IsInside(nscb, rprompt):countLgROI30[7] += 1
                                  if top55.IsInside(nscb, rprompt):countLgROI55[7] += 1
                                  if sb.IsInside(nscb, rprompt): countLgSideband[7] += 1


print "neck-roi, neck-sideband, LG-roi, LG-sb"
for  i in range(level):

    if i>
    print countNeckROI55[i]/time1*1000, countNeckSideband[i]/time1*1000, countLgROI55[i]/time1*1000, countLgSideband[i]/time1*1000
    #print "level ", i+8, countNeckROI55[i]/time1*1000, countNeckSideband[i]/time1*1000, countLgROI55[i]/time1*1000, countLgSideband[i]/time1*1000


