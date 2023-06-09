from ROOT import *
from rat import *
from scipy import interpolate
import numpy as np
import sys, os
import matplotlib.pyplot as plt
from numpy import array, pi, random, sqrt

label_data = ['tpb', 'Gar0p01', 'Gar0p015', 'Gar0p02', 'Gar1']
labelGar = label_data[1] 
print "Now we process ", labelGar
# x0 = TVector3(0,0,-850+3./1000) ## on TPB surface
x0 = TVector3(0,0,-850) ## on TPB surface
DEBUG = False
rav = 850 ## inner AV radius
rtpb = 850 - 3./1000
rAVout = 850 + 35./1000 # SRIM says only travelling 35.23 um in AV

path=os.getcwd()
file_tpb = open('helium_in_tpb.dat')
file_gar_0p01 = open('helium_in_gar_0p01bar.dat')
file_gar_0p015 = open('helium_in_gar_0p015bar.dat')
file_gar_0p02 = open('helium_in_gar_0p02bar.dat')
file_gar_1 = open('helium_in_gar_1bar.dat')

fnew = TFile("saveSrim_"+labelGar+".root", "recreate")

processfile1 = []
processfile2 = []
processfile3 = []
processfile4 = []
processfile5 = []

count = 0
for i in file_tpb:
    processfile1.append(i.rstrip())

for i in file_gar_0p01:
    processfile2.append(i.rstrip())

for i in file_gar_0p015:
    processfile3.append(i.rstrip())

for i in file_gar_0p02:
    processfile4.append(i.rstrip())

for i in file_gar_1:
    processfile5.append(i.rstrip())

deltaE = []; deltaL = [];
data1 = {'dE': deltaE, 'dL': deltaL}
data2 = {'dE': deltaE, 'dL': deltaL}
data3 = {'dE': deltaE, 'dL': deltaL}
data4 = {'dE': deltaE, 'dL': deltaL}
data5 = {'dE': deltaE, 'dL': deltaL}

wholeDataSet = {'tpb': data1, 'Gar0p01': data2, 'Gar0p015': data3, 'Gar0p02': data4, 'Gar1': data5}

def getSrimInfo(process0, label):
    E = []
    L_Eloss= []
    for line in process0:
        data = line.split(' ')
        if data[1] == 'MeV':
             E.append(float(data[0]))
        elif data[1] == 'keV':
             E.append(float(data[0])/1000)
        if data[3] == 'mm':
             L_Eloss.append(float(data[2]))
        elif data[3] == 'm':
             L_Eloss.append(float(data[2])*1000)
        elif data[3] == 'um':
             L_Eloss.append(float(data[2])/1000)
        elif data[3] == 'A':
             L_Eloss.append(float(data[2])/1e7)

    lenE = len(E)
    lenL = len(L_Eloss)
    if lenE != lenL:
        print "wrong!!!"
        return [-1,-1]

    #labels = wholeDataSet.keys() # why wrong orders??
    print "-----", label_data[label], "-----------------"
    print "Emin", E[0], "Emax", E[len(E)-1], "MeV"
    print "Lmin", L_Eloss[0], "Lmax", L_Eloss[len(L_Eloss)-1], "mm"

    deltaE = []
    deltaL = []
    for i in range(lenE):
        dE = E[lenE-1] - E[i]
        dL = L_Eloss[lenL-1] - L_Eloss[i]
        deltaE.append(round(dE, 6))
        deltaL.append(round(dL, 6))

    # save = [deltaE, deltaL]
    save = [E, L_Eloss]
    return save

save1 = getSrimInfo(processfile1, 0)
deltaE1 = save1[0]
#print "Look here", deltaE1
deltaL1 = save1[1]
data1['dE'] = np.array(deltaE1)
data1['dL'] = np.array(deltaL1)

save2 = getSrimInfo(processfile2, 1)
deltaE2 = save2[0]
deltaL2 = save2[1]
data2['dE'] = np.array(deltaE2)
data2['dL'] = np.array(deltaL2)

save3 = getSrimInfo(processfile3, 2)
deltaE3 = save3[0]
deltaL3 = save3[1]
data3['dE'] = np.array(deltaE3)
data3['dL'] = np.array(deltaL3)

save4 = getSrimInfo(processfile4, 3)
deltaE4 = save4[0]
deltaL4 = save4[1]
data4['dE'] = np.array(deltaE4)
data4['dL'] = np.array(deltaL4)

save5 = getSrimInfo(processfile5, 4)
deltaE5 = save5[0]
deltaL5 = save5[1]
data5['dE'] = np.array(deltaE5)
data5['dL'] = np.array(deltaL5)


pressure = 6.0484e+04
## load SRIM table

# direction
# randomized direction
NN = 1000
Phi = [] 
Theta = []

def CalScintPath(rSphere, direction):
    scintpathInAV = 0
    fx = x0.X()
    fy = x0.Y()
    fz = x0.Z()
    dx = direction.X()
    dy = direction.Y()
    dz = direction.Z()
    vtx2 =  x0*x0 # (x0 - oAV)*(x0 - oAV)
    rVertex = sqrt( vtx2 ); # Roffset = |X0 - oAV|
    udDotVtx = x0*direction
    sqrVal = udDotVtx*udDotVtx - vtx2 + rSphere*rSphere
    a1 = 0; a2 = 0; aplus = 0; abig = 0; asmall = 0;
    ## print sqrVal
    if sqrVal<0:
         print "ray-sphere no interception"
    else: # line passes AV sphere
         # find the line-sphere intersect points; a1, a2 are the path lengths btw vertex and intersection points
         a1 = -udDotVtx + sqrt(sqrVal);
         a2 = -udDotVtx - sqrt(sqrVal); # always a2<a1
         if DEBUG:
             print "ray-sphere a1 =", a1, "a2=", a2
         if( rVertex<rSphere): # vertex inside the AV, equivalent to a1*a2<0
             if ( a1*a2<0 ):
                 aplus = a1; # a1>0>a2
                 scintpathInAV = aplus;
         ## vertex in AV
         else:  # rVertex>=rSphere, vertex in external
             if( a1>=0 and a2>=0):
                 abig = a1; # far intersection point
                 asmall = a2; # near intersection point
                 scintpathInAV = abig - asmall;
           # ensure a1 and a2 are positive
         # vertex in external
         # pass through AV
    return scintpathInAV

for i in range(NN):
    phi = random.uniform(-pi, pi) 
    theta = random.uniform(0, pi) 
    Phi.append(phi)
    Theta.append(theta)

hDistInScint = TH1F("hDistInScint", "geometric path in scint(GAr)", 1700,0,1700)
hDistInAv = TH1F("hDistInAv", "geometric path in Av sphere", 1700,0,1700)
hDistInTpb = TH1F("hDistInTpb", "geometric path in TPB bulk", 2000,0,2)

hTotalDist1 = TH1F("hTotalDist1", "tpb surf vertex, path with Eloss, scintPath+tpbPath", 1700,0,1700)
hTotalDist2 = TH1F("hTotalDist2", "AV surf vertex, path with Eloss, tpbPath+scintPath", 1700,0,1700)

hEtpb = TH1F("hEtpb", "energy deposited in TPB sphere", 100,0,5.3)
## only for pure tpb path!!
hEtpbVsGeoDistInTpb = TH2F("hEtpbVsGeoDistInTpb", "geometric path vs energy deposited in TPB sphere", 1000,0,1, 100,0,5.3)
hEgar = TH1F("hEgar", "energy deposited in tpb sphere", 100,0,5.3)
hEgarVsGeoDistInGar= TH2F("hEgarVsGeoDistInGar", "geometric path vs energy deposited in"+labelGar, 1700,0,1700, 100,0,5.3)

hEvsTotalDist_tpbSurf = TH2F("hEvsTotalDist_tpbSurf", "tpb surf, total distance (with Eloss) vs energy deposited in"+labelGar, 1700,0,1700, 100,0,5.3)
hEvsTotalDist_avSurf = TH2F("hEvsTotalDist_avSurf", "av surf, total distance (with Eloss) vs energy deposited in"+labelGar, 1700,0,1700, 100,0,5.3)

hEdepositTpb0 = TH1F("hEdepositTpb0", "energy deposited for TPB only track, all near-side", 100,0,5.3)
hEdepositTpb0VsDist = TH2F("hEdepositTpb0VsDist", "energy deposited totally in TPB vs. distance, for TPB only track, all near-side", 1700, 0, 1700, 100, 0, 5.3)

### only for Av surf events
hEdepositTpb1 = TH1F("hEdepositTpb1", "energy deposited in near-side TPB", 100, 0, 5.3)
hEdepositTpb2 = TH1F("hEdepositTpb2", "energy deposited in far-side TPB", 100, 0, 5.3)
hEdepositGAr = TH1F("hEdepositGAr", "energy deposited in GAr", 100, 0, 5.3)

hEdepositTpb1VsDist = TH2F("hEdepositTpb1VsDist", "energy deposited in near-side TPB vs distance", 1700, 0, 1700, 100, 0, 5.3)
hEdepositTpb2VsDist = TH2F("hEdepositTpb2VsDist", "energy deposited in far-side TPB vs distance", 1700, 0, 1700, 100, 0, 5.3)
hEdepositGArVsDist = TH2F("hEdepositGArVsDist", "energy deposited in GAr vs distance", 1700, 0, 1700, 100, 0, 5.3)

for i in range(NN):
    phi = Phi[i]
    theta = Theta[i]
    uu = TVector3( cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta) )

    print "--- ---- eventID ", i, "X0 = (", x0.X(), ",", x0.Y(), ",", x0.Z(), "), u=(", round(uu.X(),4), ",", round(uu.Y(),4), ",", round(uu.Z(),4),")"
    lengthInScint = CalScintPath(rtpb, uu) # to the sphere of TPB surface, Lp_scint
    lengthInAv = CalScintPath(rav, uu) # to the sphere of AV surface, Lp_Av, then Lp_tpb = Lp_Av - Lp_scint

    print "LAv=", lengthInAv, ", Lscint=", lengthInScint, "Ltpb =", lengthInAv - lengthInScint

    if lengthInAv <= 0: # never care about the rays not passing through the inner AV sphere
        continue

    if lengthInScint>0: 
        hDistInScint.Fill(lengthInScint)
    #print "geometrical travel length in GAr, r<r_tpb, path length", lengthInScint

    lengthInTpb = 0
    if lengthInAv>0:
        lengthInTpb = lengthInAv - lengthInScint

    if lengthInTpb>0:
        hDistInTpb.Fill(lengthInTpb)

    #print lengthInScint, lengthInAv, lengthInTpb
    ## extrapolate energies!!!
    if x0.Mag()<rav: ## for the vertex at tpb surf
        if lengthInAv>0:
            hTotalDist1.Fill(lengthInAv)
    else: ## for the vertex at tpb surf
        if lengthInAv>0:
            hTotalDist2.Fill(lengthInAv)

    ## !!! Now calculate energy, energy loss and paths with energy loss
    Egar = -9999
    Etpb = -9999
    actualLengthInScint = 0 # path length with energy loss
    actualLengthInTpb = 0 
    Ear2tpb = -9999 # for tpb surf case, the energy step when alpha goes from agron to tpb
    Etpb2ar = -9999 # for Av surf case, the energy step when alpha goes from tpb to argon
    Eactual = 5.3 # actual energy at the last step

    if x0.Mag()<rav: # for tpb surface vertex
        if lengthInAv>0:
           ddL = wholeDataSet[labelGar]['dL']
           ddE = wholeDataSet[labelGar]['dE']
           #print ddL
           #print ddE
           if lengthInScint>ddL[0]: ## the shortest length in GAr when E at 10 keV
               try:
                  ff = interpolate.interp1d(ddL, ddE, kind = 'linear')
                  Egar = round( ff(lengthInScint),6 )
                  print "interpolate E in Gar", Egar, "MeV", "Eactual = 5.3-", Egar, "=", 5.3-Egar, "MeV"
                  if Egar>=0 and Egar<=5.3:
                       Ear2tpb = 5.3 - Egar
                       Eactual = 5.3 - Egar
                       ff1 = interpolate.interp1d(ddE, ddL, kind = 'linear')
                       actualLengthInScint = round( ff1(Eactual),6 ) ## with this energy, how much long the alpha can continue to travel ignoring the geometric path
                       if actualLengthInScint>lengthInScint:
                           print "geo length in scint", lengthInScint, "<actualLengthInAv", round(actualLengthInScint,6), "set as", round(lengthInScint,6)
                           actualLengthInScint = lengthInScint
                       print "possible length in Gar due to energy", actualLengthInScint, "mm"
                  # actualLength = ff(Etpb)
               except:
                  print "failed to find interpolation"
                  Egar = 0
                  actualLengthInScint = 0
           else:
               Egar = 0
               actualLengthInScint = 0
               #print "path too long, stop. all absorbed. E = 0"
           #print "length in GAr", lengthInScint, "mm. E=", Egar, "MeV; actual length", actualLength

           if Egar>0 and Egar<=5.3:
               Eactual = 5.3 - Egar
               hEgar.Fill(Eactual)
               hEgarVsGeoDistInGar.Fill(lengthInScint, Eactual)
               if lengthInTpb<0: # rarely happens!!
                   print "are you serious? why the ray crosses the scint but not tpb? strange geometry!!"
                   hEvsTotalDist_tpbSurf.Fill(actualLengthInScint, Eactual)

        if lengthInTpb>0:
           ddL = wholeDataSet['tpb']['dL']
           ddE = wholeDataSet['tpb']['dE']
           if lengthInScint <= 1e-12:# for the tpb surface case, if there is no scint paths but only tpb paths
               print "Possible TPB only path, LAv=", lengthInAv, "Lscint =", lengthInScint
               if lengthInTpb>ddL[0]:
                   try:
                      ff = interpolate.interp1d(ddL, ddE, kind = 'linear')
                      Etpb = round( ff(lengthInTpb),6 )
                      Eactual = 5.3 - Etpb
                      print "interpolate Etpb", Etpb, "MeV, Eactual=5.3-", Etpb, "=", Eactual
                      ff1 = interpolate.interp1d(ddE, ddL, kind = 'linear')
                      actualLengthInTpb = round( ff1(Eactual), 6)
                      if actualLengthInTpb>lengthInTpb:
                         print "need to trunct the actualLength:", actualLengthInTpb, ">", lengthInTpb
                         actualLengthInTpb = lengthInTpb
                      print "Eactual=", Eactual, "possible length in tpb due to energy", actualLengthInTpb, "mm"
                   except:
                      print "TPB only track. Etpb failed to find interpolation"
               else:
                   print "Tpb length is less than SRIM length for 10keV", lengthInTpb 
                   Etpb = 0
                   actualLengthInTpb = 0
                   #print "path too long, stop. all absorbed. E = 0"
               # print "length in tpb", lengthInTpb*1000, "um. E=", Etpb, "MeV; actual length", actualLength
               if Etpb>0:
                  Eactual = 5.3 - Etpb
                  hEtpb.Fill(Eactual)
                  hEtpbVsGeoDistInTpb.Fill(lengthInTpb, Eactual)
                  hEvsTotalDist_tpbSurf.Fill(actualLengthInTpb, Eactual)
           else: # lengthInScint>0, alpha crosses the Argon, need to check how much energy lost in Argon
               print "Tpb surf vertex. path crosses scint and then tpb"
               if Eactual<0.01:#<10 keV
                    print "E argon to tpb is too small"
                    actualLengthInTpb = 0
                    hEtpbVsGeoDistInTpb.Fill(lengthInTpb, Eactual) ## E lost only in GAr
                    hEvsTotalDist_tpbSurf.Fill(actualLengthInScint, Eactual)
               else:
                    if lengthInTpb>ddL[0]:
                        try:
                            ff = interpolate.interp1d(ddL, ddE, kind = 'linear')
                            Etpb = round( ff(lengthInTpb),6 )
                            Eactual = 5.3 - Egar - Etpb
                            print "interpolate Etpb loss here", Etpb, "MeV", "Eactual = 5.3-Egar-Etpb=5.3-",Egar,"-", Etpb, "=", Eactual, "MeV"
                            ff1 = interpolate.interp1d(ddE, ddL, kind = 'linear')
                            actualLengthInTpb = round( ff1(Eactual), 6)
                            if actualLengthInTpb>lengthInTpb:
                                actualLengthInTpb = lengthInTpb
                            print "possible length in tpb due to energy", actualLengthInTpb, "mm"
                        except:
                            print "Etpb failed to find interpolation"
                    else:
                        Etpb = 0
                        lengthInTpb = 0

                    if Eactual>=0:
                        hEtpb.Fill(Etpb)
                        hEtpbVsGeoDistInTpb.Fill(lengthInScint+lengthInTpb, Eactual)
                        hEvsTotalDist_tpbSurf.Fill(actualLengthInScint+actualLengthInTpb, Eactual)

    else: # for the Av surf vertex
        print "vertex at av surface"
        if lengthInAv>0: # if not, the ray can go into AV bulk, no scintillation, don't care!!
           # load GAr data 
           ddL = wholeDataSet[labelGar]['dL']
           ddE = wholeDataSet[labelGar]['dE']
           # load tpb data
           ddL1 = wholeDataSet['tpb']['dL']
           ddE1 = wholeDataSet['tpb']['dE']
           Eactual = 5.3
           actualLengthInTpb   = 0 # tpb only track
           actualLengthInTpb1  = 0 # near tpb
           actualLengthInTpb2  = 0 # far tpb
           actualLengthInScint = 0 # track in scint
           ### Track type 1: paths only in TPB, single-clusters
           if lengthInScint<=0 and lengthInTpb>ddL1[0]: #1.5 um
               print "!!!!!!!!!!!!! tpb track only? that's rare !!!", lengthInAv, lengthInScint, lengthInTpb
               try:
                  ff = interpolate.interp1d(ddL1, ddE1, kind = 'linear')
                  Etpb = round( ff(lengthInTpb),6 )
                  if Etpb >=0:
                       Eactual = 5.3 - Etpb
                       ff1 = interpolate.interp1d(ddE1, ddL1, kind = 'linear')
                       actualLengthInTpb= round( ff1(Eactual),6 ) ## with this energy, how much long the alpha can continue to travel ignoring the geometric path
                       if actualLengthInTpb>lengthInTpb:
                           print "geo length in tpb", lengthInTpb, "<actualLengthInTpb", actualLengthInTpb
                           actualLengthInTpb = lengthInTpb
                       print "possible length in Gar due to energy", round(actualLengthInScint,6), "mm with E=", Eactual,"MeV" 
                  # actualLength = ff(Etpb)
               except:
                  print "failed to find interpolation"
                  Eactual = 0 
                  actualLengthInTpb = 0 
               if Eactual>0 and actualLengthInTpb>0:
                  hEdepositTpb0.Fill(Eactual)
                  hEdepositTpb0VsDist.Fill(actualLengthInTpb, Eactual) 

           ### Track type 2: paths go through tpb and scint volume, double-clusters
           if lengthInScint>ddL[0] and lengthInTpb>ddL1[0]:
               ### tpb length divided into a half, the near path and the far path!!
               ### first/near tpb track!
               print "now path goes through the near tpb"
               l1 = lengthInTpb/2; l2 = lengthInTpb/2
               Eactual1 = 0 # near tpb
               Eactual2 = 0 # far tpb
               Escint = 0
               if l1>ddL1[0]:
                   try:
                      ff = interpolate.interp1d(ddL1, ddE1, kind = 'linear')
                      Etpb = round( ff(l1),6 )
                      if Etpb >=0:
                           Eactual1 = 5.3 - Etpb
                           ff1 = interpolate.interp1d(ddE1, ddL1, kind = 'linear')
                           actualLengthInTpb1 = round( ff1(Eactual1),6 ) ## with this energy, how much long the alpha can continue to travel ignoring the geometric path
                           if actualLengthInTpb1 > l1:
                               print "geo length in TPB: ", round(l1,6), "< actualLengthInTpb1:", round(actualLengthInTpb1,6), "set to", round(l1,6)
                               actualLengthInTpb1 = l1
                           print "type 2, possible length in Gar due to energy ", round(actualLengthInTpb1,6), "mm. E1 =", Eactual1, "MeV"
                      # actualLength = ff(Etpb)
                   except:
                      print "failed to find interpolation"
                      Eactual1 = 0 
                      actualLengthInTpb1 = 0 
                   if Eactual1>0 and actualLengthInTpb1>0:
                      ### save the far tpb points 
                      hEdepositTpb1.Fill(Eactual1)
                      hEdepositTpb1VsDist.Fill(actualLengthInTpb1, Eactual1)
                      ### now checking the path in GAr, the energy begin at Eactual1
                      try:
                         ### evaluate the path in LAr
                         ff = interpolate.interp1d(ddE, ddL, kind = 'linear')
                         actualLengthInScint = round( ff(Eactual1), 6 )
                         if actualLengthInScint>lengthInScint:
                             actualLengthInScint = lengthInScint
                      except:
                          print "failed to find interpolation"
                          Escint = 0 
                          actualLengthInScint = 0
                      ## now evaluate the energy loss in the GAr region
                      print "Now pass the near TPB and enter the GAr. E=", Eactual2
                      if actualLengthInScint>ddL[0]:
                          try:
                              ff = interpolate.interp1d(ddL, ddE, kind = 'linear')
                              E0 = round( ff(actualLengthInScint), 6)
                              Escint = Eactual1 - E0
                          except:
                              Escint = 0
                          if Escint>0:
                              hEdepositGAr.Fill(Escint)
                              hEdepositGArVsDist.Fill(actualLengthInScint, Escint)
                              ### now check whether alpha goes into far-side TPB
                              try:
                                  ### evaluate the path in TPB
                                  ff1 = interpolate.interp1d(ddE1, ddL1, kind = 'linear')
                                  actualLengthInTpb2 = round( ff1(Escint), 6)
                              except:
                                  actualLengthInTpb2 = 0
                              print "Now pass the GAr and enter the far TPB"
                              if actualLengthInTpb2>ddL1[0]:
                                  try:
                                      ff = interpolate.interp1d(ddL1, ddE1, kind = 'linear')
                                      Etpb = round( ff(actualLengthInTpb2), 6)
                                      Eactual2 = Escint - Etpb 
                                  except:
                                      Eactual2 = 0

                                  if Eactual2>0:
                                      hEdepositTpb1.Fill(Eactual2)
                                      hEdepositTpb1VsDist.Fill(actualLengthInTpb2, Eactual2)

fnew.cd()
hDistInScint.Write()
hTotalDist1.Write()
hTotalDist2.Write()
hDistInTpb.Write()
hEtpb.Write()
hEtpbVsGeoDistInTpb.Write()

hEgar.Write()
hEgarVsGeoDistInGar.Write()

hEvsTotalDist_tpbSurf.Write()
hEvsTotalDist_avSurf.Write()

hEdepositTpb0.Write()
hEdepositTpb0VsDist.Write()
hEdepositGAr.Write()
hEdepositGArVsDist.Write()
hEdepositTpb1.Write()
hEdepositTpb1VsDist.Write()
hEdepositTpb2.Write()
hEdepositTpb2VsDist.Write()


c1 = TCanvas("c1","dist in scint",600,400)
c1.cd()
hDistInScint.Draw()

c2 = TCanvas("c2","dist in Av",600,400)
c2.cd()
hDistInTpb.Draw()
#hDistInAv.Draw()

c3 = TCanvas("c3","TPB surf events, actual E vs actual distance",600,400)
c3.cd()
hEvsTotalDist_tpbSurf.Draw("colz")

c4 = TCanvas("c4","dist in TPB",600,400)
c4.cd()
hEtpbVsGeoDistInTpb.Draw("colz")

c5 = TCanvas("c5","dist in "+labelGar,600,400)
c5.cd()
hEgarVsGeoDistInGar.Draw("colz")

raw_input()
