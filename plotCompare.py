import ROOT
#from __future__ import print_function
from ROOT import TCanvas, TGraph, TH1F
from ROOT import gROOT
from array import array
from ROOT import *
fmpw = []
fmpw.append(ROOT.TFile("ResolMCmpw_FitMultiWater_1p40_Water_6MeV_10000evt_ISOFill.root"))
fmpw.append(ROOT.TFile("ResolMCmpw_FitMPW_mode75_TresCut_Water_6MeV_10000evt_ISOFill.root"))
fmpw.append(ROOT.TFile("ResolMCmpw_FitMPW_mode75_TresCut3_Water_6MeV_10000evt_ISOFill.root"))
cDeltaRvsTime = TCanvas( 'cDeltaRvsTime', 'Delta R vs Time', 200, 10, 1000, 1000 )
cDeltaRvsTime.SetGrid()
gpad = ROOT.gPad
gstyle = ROOT.gStyle
gstyle.SetOptFit(1011)
hDeltaX = [fmpw[i].Get("hDeltaX") for i in range(3)]
hDeltaY = [fmpw[i].Get("hDeltaY") for i in range(3)]
hDeltaZ = [fmpw[i].Get("hDeltaZ") for i in range(3)]
hT = [fmpw[i].Get("htRes") for i in range(3)]

cDeltaRvsTime.Divide(2,2)
cDeltaRvsTime.cd(1)
gpad.SetLogy()
gpad.SetGrid()
#hDeltaX[0].GetXaxis().SetTitle("|fitPos-mcPos|")
hDeltaX[0].Draw()
hDeltaX[1].SetLineColor(4)
hDeltaX[1].Draw("sames")
hDeltaX[2].SetLineColor(2)
hDeltaX[2].Draw("sames")

cDeltaRvsTime.cd(2)
gpad.SetLogy()
gpad.SetGrid()
#hDeltaY[0].GetXaxis().SetTitle("|fitPos-mcPos|")
hDeltaY[0].Draw()
hDeltaY[1].SetLineColor(4)
hDeltaY[1].Draw("sames")
hDeltaY[2].SetLineColor(2)
hDeltaY[2].Draw("sames")

cDeltaRvsTime.cd(3)
gpad.SetLogy()
gpad.SetGrid()
#hDeltaZ[0].GetXaxis().SetTitle("|fitPos-mcPos|")
hDeltaZ[0].Draw()
hDeltaZ[1].SetLineColor(4)
hDeltaZ[1].Draw("sames")
hDeltaZ[2].SetLineColor(2)
hDeltaZ[2].Draw("sames")

cDeltaRvsTime.cd(4)
gpad.SetLogy()
gpad.SetGrid()
#hT[0].GetXaxis().SetTitle("|fitPos-mcPos|")
hT[0].Draw()
hT[1].SetLineColor(4)
hT[1].Draw("sames")
hT[2].SetLineColor(2)
hT[2].Draw("sames")



cDeltaRvsTime.Draw()
text = raw_input()
