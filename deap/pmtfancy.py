#!/usr/bin/python
import os,sys
from math import pi, acos
import numpy as np
import ROOT

pmt_radius = 85.1
lg_radius = 9.5

pmtpos = './Aksel_PMTpos_260.ratdb'
if not os.path.isfile(pmtpos):
    print 'ERROR: need Aksel_PMTpos_260.ratdb to work'
    print 'INFO: inside TH2DPMTfancy.py, edit the variable pmtpos to the right place'
    sys.exit(1)

class pmtpositions:
    def __init__(self, fn=pmtpos):
        self.fn = fn
        self.xpos, self.ypos, self.zpos = None, None, None
        self.pmtvec = []
    def Load(self):
#        print 'debug',self.fn
        self.f = file(self.fn,'r')
        for l in self.f:
            if l.startswith('x'): self.xpos = [float(val) for val in l.split('[')[1][:-4].split(',')]
            if l.startswith('y'): self.ypos = [float(val) for val in l.split('[')[1][:-4].split(',')]
            if l.startswith('z'): self.zpos = [float(val) for val in l.split('[')[1][:-4].split(',')]
        if None in [self.xpos, self.ypos, self.zpos]:
            print 'ERROR: could not read valid PMT position data from', self.fn
            sys.exit(1)

        for i,(x,y,z) in enumerate(zip(self.xpos, self.ypos, self.zpos)):
            v = ROOT.TVector3(x,y,z)
            self.pmtvec.append(pmt_radius*v.Unit())

        return True

pp = pmtpositions()
pp.Load()

class TH2DPMTfancy:
    def __init__(self, histoname, histotitle, nbinsx=150):
        self.h = ROOT.TH2D(histoname, histotitle+';#phi [rad];cos#theta', nbinsx, -pi, pi, nbinsx, -1., 1.)
        self.h.SetStats(0)
        # self.h.GetXaxis().CenterTitle()
    def Fill(self, userdata):
        if type(userdata)==type([]):
            data_of_length255 = userdata
        elif type(userdata)==type({}):
            data_of_length255 = 255*[0.]
            for i in range(255):
                if i in userdata.keys():
                    data_of_length255[i] = userdata[i]

        vecbin = ROOT.TVector3()
        hptr = self.h
        counter = 0
        for iphi in range(0, hptr.GetNbinsX()):
            if not iphi%10: print 'iphi=',iphi,'of',hptr.GetNbinsX()
            phibin = hptr.GetXaxis().GetBinUpEdge(iphi)
            for itheta in range(0, hptr.GetNbinsY()):
                thetabin = acos(hptr.GetYaxis().GetBinUpEdge(itheta))
                vecbin.SetMagThetaPhi(pmt_radius,thetabin,phibin)
                for pmtid in range(255):
                    if data_of_length255[pmtid]==0: continue
                    # if (pp.pmtvec[pmtid]-vecbin).Mag()<lg_radius:
                    if (pp.pmtvec[pmtid]-vecbin).Mag()<8.51:
                        hptr.SetBinContent(iphi, itheta, data_of_length255[pmtid]*1000)
                        # hptr.SetBinContent(iphi, itheta, 2)
                        counter+=1
        print(counter)
        fout = ROOT.TFile("test.root", "RECREATE")
        self.h.Write(self.h.GetName())
    
    def Draw(self, dopt):
        # self.h.SetMarkerStyle(1)
        self.h.SetLineColor(15)
        self.h.SetLineWidth(1)
        self.h.SetLineStyle(4)
        # self.h.SetMarkerSize(10)
        # self.h.SetMarkerColor(2)
        self.h.Draw(dopt)

## Jie: NOTE: run the code below first to create test.root, then run the next __main__
#if __name__=='__main__':
#    h = TH2DPMTfancy('h','test', 500)
#    # h = TH2DPMTfancy('h','test', 150)
#    myda = [x+1 for x in range(255)]
#    #uncomment to load 1d histogram with 255 bins from input file
#    #ipf = ROOT.TFile(ipfn)
#    #h_inputdata = ipf.Get('histogram_name')
#    #myda = [h_inputdata.GetBinContent(ibin+1) for ibin in range(255)]
#    Canvas = ROOT.TCanvas ("Canvas", "Canvas", 1200, 1000);           Canvas.SetTickx(1);          Canvas.SetTicky(1);     Canvas.SetLogz(1);    Canvas.SetLogy(0);   Canvas.SetGrid(0,0)
#    h.Fill(myda)
#    h.Draw('CONT0');     Canvas.SaveAs("test0.png"); Canvas.Clear()
#    h.Draw('CONT1');     Canvas.SaveAs("test1.png"); Canvas.Clear()
#    h.Draw('CONT2');     Canvas.SaveAs("test2.png"); Canvas.Clear()
#    h.Draw('CONT3');     Canvas.SaveAs("test3.png"); Canvas.Clear()
#    h.Draw('CONT4');     Canvas.SaveAs("test4.png"); Canvas.Clear()
#    # h.Draw('SURF2');     Canvas.SaveAs("test.png"); Canvas.Clear()
#    # au = raw_input('>')
#    Canvas.Close()

if __name__=='__main__':

    fin = ROOT.TFile("test.root")
    h = fin.Get('h'); h.SetDirectory(0);
    fin.Close()
    h.SetLineColor(15)
    h.SetLineWidth(1)
    h.SetLineStyle(4)
    h.GetXaxis().CenterTitle(0)
    h.SetTitle('')

    Canvas = ROOT.TCanvas ("Canvas", "Canvas", 1400, 1000);Canvas.SetLogz(1);Canvas.SetLogy(0);Canvas.SetGrid(0,0)
    Canvas.SetTickx(1);Canvas.SetTicky(0); 
    # for pmtid in range(255):
    # h.Draw('CONT0');     Canvas.SaveAs("test0.png"); Canvas.Clear()
    # h.Draw('CONT1');     Canvas.SaveAs("test1.png"); Canvas.Clear()
    # h.Draw('CONT2');     Canvas.SaveAs("test2.png"); Canvas.Clear()
    h.Draw('CONT3')
    for pmtid in range(255):
        x,y,z = pp.pmtvec[pmtid][0], pp.pmtvec[pmtid][1], pp.pmtvec[pmtid][2]
        Tl = ROOT.TLatex ()
        Tl.SetTextAlign(22)
        Tl.SetTextSize(0.02)

        phi = np.arctan(y/x)
        if (x<0 and y<0):
            phi -=np.pi
        if (x<0 and y>0):
            phi +=np.pi
        costheta = (z/np.sqrt(x**2+y**2+z**2)).tolist()

        if phi>3.10:
            phi=3.10
        if phi<-3.10:
            phi = -3.10
        # Tl.DrawLatex(phi, costheta, '#color[4]{'+str(pmtid)+"}")
        Tl.DrawLatex(phi, costheta, str(pmtid))

    phis, costhetas = [-1.9, -0.63, 0.62, 1.9, 3.10, -2.5, -1.25, 0, 1.25, 2.5, -0.5], [0.45]*5+[-0.45]*5 + [-1.05]
    for i in range(11):
        phi, costheta = phis[i], costhetas[i]
        Tl = ROOT.TLatex()
        Tl.SetTextAlign(22)
        Tl.SetTextSize(0.04)
        Tl.DrawLatex(phi, costheta, '#color[2]{'+str(i)+'}')

    Canvas.SaveAs("PMTs.png"); Canvas.Clear()

    Canvas.Close()
