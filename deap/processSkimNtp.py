# contact: Jie Hu (jhu9@ualberta.ca)
# date: 10 July 2023
# use PyROOT, apply to the ntuple files. example:
#######################
import os
import sys
import glob
#from rat import *
from ROOT import *
import numpy as np
import struct
#from run_numbers import runs
from datetime import datetime, time as datetime_time, timedelta

plrNSCBcut1 = 90  ## additional for PLR roi
plrNSCBcut2 = 200 ## additional for PLR roi

#Run            = sys.argv[1]
file0 = str(sys.argv[1])
fileName = os.path.basename(file0)

ffile = TFile(file0)
data = ffile.Get("data_mc")
nentries = data.GetEntries()
update = 1000000

fout = TFile("skimmedAr39_"+fileName,"recreate")

dstreeclone = data.CloneTree(0);
dstreeclone.SetDirectory(fout);

for event in range (nentries):
    if (event+1)%update == 0:
        print (event+1), "Analyzed..."

    data.GetEntry(event)

    dtmTrigSrc = data.dtmTrigSrc
    calcut = data.calcut
    deltat = data.deltat
    numEarlyPulses = data.numEarlyPulses
    subeventN = data.subeventN
    eventTime = data.eventTime

    tf2x = data.timefit2X
    tf2y = data.timefit2Y
    tf2z = data.timefit2Z
    mbx = data.mblikelihoodX
    mby = data.mblikelihoodY
    mbz = data.mblikelihoodZ
    mbr = data.mblikelihoodR

    mbPos = TVector3(mbx, mby, mbz)
    tf2Pos = TVector3(tf2x, tf2y, tf2z)
    qpe = data.qPE
    fprompt = data.fprompt
    neckVetoN = data.neckVetoN
    nscb = data.nSCBayes
    rprompt60Bayes = data.rprompt60Bayes
    fmaxpe = data.fmaxpe

    promptPE =  nscb*rprompt60Bayes;

    ### Contency based cut, reconstructed Z (90% Ar39 acceptance)
    zContF = TFile("tf2mb_nSCBayes_deltaz_contours.root","READ");
    zCut = zContF.Get("cont90_cut;1");

    rContF = TFile("tf2mb_nSCBayes_dist_after_deltaz90_contours.root","READ");
    rCut = rContF.Get("cont85_cut;1");

    chargetopring = data.chargetopring
    chargesecondring = data.chargesecondring
    chargebottomring = data.chargebottomring
    chargesecondbottomring = data.chargesecondbottomring
    chargethirdbottomring = data.chargethirdbottomring

    if not(qpe>0): continue
    cft2r = (chargetopring + chargesecondring)/qpe
    cfb3r = (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qpe
 
    pmtMaxpe = data.dbCherenkov_pmtMaxpe
    dist1 = data.dbCherenkov_dist1
    dist2 = data.dbCherenkov_dist2
    dist3 = data.dbCherenkov_dist3
    dist4 = data.dbCherenkov_dist4
    dist5 = data.dbCherenkov_dist5

    checkAr39 = nscb>plrNSCBcut1 and nscb<plrNSCBcut2
    #checkLevel1 = (dtmTrigSrc&0x82 == 0)
    #checkLevel2 = checkLevel1 and (calcut&0x31f8 == 0)
    ## for MC, cut begins here
    checkLevel3 = (deltat>20000)
    #checkLevel3 = checkLevel2 and (deltat>20000)
    checkLevel4 = checkLevel3 and (numEarlyPulses <= 3)
    checkLevel5 = checkLevel4 and (subeventN == 1)
    checkLevel6 = checkLevel5 and (2250<eventTime and eventTime<2700)
    checkLevel7 = checkLevel6 and (qpe>0) #Jie used (qpe>60), and Spencer not used this cut
    checkLevel8 = checkLevel7 and (fmaxpe<0.4)
    checkLevel9 = checkLevel8 and (neckVetoN == 0)
    #checkLevel10 = checkLevel9 and (mbr<720)
    checkLevel10 = checkLevel9 and (cft2r<0.04 and mbz<550)
    checkLevel11 = checkLevel10 and (cfb3r<0.1)
    checkLevel12 = checkLevel11 and zCut.IsInside(promptPE,(tf2Pos.Z()-mbPos.Z()))
    checkLevel13 = checkLevel12 and rCut.IsInside(promptPE,(tf2Pos-mbPos).Mag())

    checkDCC = (not( pmtMaxpe > 0.06 and (dist1 > 600. or dist2 > 600. or dist3 > 600. or dist4 > 600. or dist5 > 600.) ))

    if checkLevel10 and checkAr39 and checkDCC:
        dstreeclone.Fill()

fout.cd()
dstreeclone.Write()
