import ROOT
import sys, os, getopt
from ROOT import TCanvas, TTree, TFile, TH1F, AddressOf
import numpy as np
import csv
from array import array
import random
import matplotlib.pyplot as plt
from random import choice
fname = 'saveNCherenPhoton.root'
binN = 240
ff = TFile(fname,"READ")
h1 = ff.Get("hNCherPhoton2MeV")
h1.Scale(1./h1.Integral())
h2 = ff.Get("hNCherPhoton3MeV")
h2.Scale(1./h2.Integral())
h3 = ff.Get("hNCherPhoton4MeV")
h3.Scale(1./h3.Integral())
h4 = ff.Get("hNCherPhoton5MeV")
h4.Scale(1./h4.Integral())
h5 = ff.Get("hNCherPhoton6MeV")
h5.Scale(1./h5.Integral())
h6 = ff.Get("hNCherPhoton7MeV")
h6.Scale(1./h6.Integral())
h7 = ff.Get("hNCherPhoton8MeV")
h7.Scale(1./h7.Integral())
h8 = ff.Get("hNCherPhoton9MeV")
h8.Scale(1./h8.Integral())
h9 = ff.Get("hNCherPhoton10MeV")
h9.Scale(1./h9.Integral())
list_ni = []
list_ni_2MeV = [];list_pi_2MeV = []
list_ni_3MeV = [];list_pi_3MeV = []
list_ni_4MeV = [];list_pi_4MeV = []
list_ni_5MeV = [];list_pi_5MeV = []
list_ni_6MeV = [];list_pi_6MeV = []
list_ni_7MeV = [];list_pi_7MeV = []
list_ni_8MeV = [];list_pi_8MeV = []
list_ni_9MeV = [];list_pi_9MeV = []
list_ni_10MeV = [];list_pi_10MeV = []

nphoton = 0
xbin = 6000./240
for i in range(binN):
    n_i = nphoton + i*xbin
    list_ni.append(n_i)
    p_i = h1.GetBinContent(i+1)
    list_pi_2MeV.append(p_i)
    p_i = h2.GetBinContent(i+1)
    list_pi_3MeV.append(p_i)
    p_i = h3.GetBinContent(i+1)
    list_pi_4MeV.append(p_i)
    p_i = h4.GetBinContent(i+1)
    list_pi_5MeV.append(p_i)
    p_i = h5.GetBinContent(i+1)
    list_pi_6MeV.append(p_i)
    p_i = h6.GetBinContent(i+1)
    list_pi_7MeV.append(p_i)
    p_i = h7.GetBinContent(i+1)
    list_pi_8MeV.append(p_i)
    p_i = h8.GetBinContent(i+1)
    list_pi_9MeV.append(p_i)
    p_i = h9.GetBinContent(i+1)
    list_pi_10MeV.append(p_i)

#hNew2MeV = TH1F("hNew2MeV","",binN,0,6000)
binNN = 100
hProbE = TH1F("hProbE","",binNN,0,10)
#for seed in range(10000):
histList = []
hTot = TH1F("hTot","",binNN,0,10)
seedN = 10000
for seed in range(seedN):
 histList.append(hProbE)
 xx = 0
 kk = 0
 for loop in range(binNN):
   # bin = []
   if(xx<=3 and xx>=2):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_2MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
       kk = randomN
   if(xx<=4 and xx>3):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_3MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
   if(xx<=5 and xx>4):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_4MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   if(xx<=6 and xx>5):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_5MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   if(xx<=7 and xx>6):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_6MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   if(xx<=8 and xx>7):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_7MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   if(xx<=9 and xx>8):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_8MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   if(xx<=10 and xx>9):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_9MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   if(xx>10):
       randomN = np.random.choice(np.array(list_ni), p=list_pi_10MeV)
       histList[seed].SetBinContent(loop+1,randomN)
       #hProbE.SetBinContent(loop,randomN)
       #hProbE.Fill(xx,randomN)
   xx = xx + 10./binNN

for seed in range(seedN):
  hTot.Add(histList[seed])


hTot.Scale(1./seedN)
#random_number = np.random.choice(np.array(list_ni), p=list_pi_2MeV)
#hNew2MeV.Fill(random_number)

listTable = []

for i in range(binNN):
  listTable.append(int(hTot.GetBinContent(i+1)))

print listTable

#c1 = TCanvas("c1","",500,400)
#hTot.Draw()
##hNew2MeV.Draw()   
#c1.Draw()
#raw_input()
