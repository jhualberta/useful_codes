#!/usr/bin/python
import os,sys
from math import pi, acos
import numpy as np
import ROOT
from ROOT import *
from array import array
pmt_radius = 85.1
lg_radius = 9.5
#checkEventID = 37

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

class TH2DPMTfancy: ### draw fancy PMTs
    def __init__(self, histoname, histotitle, nbinsx=150):
        self.h = ROOT.TH2D(histoname, histotitle+';#phi [rad];cos#theta', nbinsx, -pi, pi, nbinsx, -1., 1.)
        self.h.SetStats(0)
        # self.h.GetXaxis().CenterTitle()
    def Fill(self, userdata):
        if type(userdata)==type([]): ## userdata is list
            data_of_length255 = userdata
        elif type(userdata)==type({}): ## userdata is dictionary
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
                    if (pp.pmtvec[pmtid]-vecbin).Mag()<8.51:#lg_radius:
                    #if (pp.pmtvec[pmtid]-vecbin).Mag()<8.51: ### very small one?
                        #hptr.SetBinContent(iphi, itheta, data_of_length255[pmtid]*1000)
                        hptr.SetBinContent(iphi, itheta, data_of_length255[pmtid]*10)
                        # hptr.SetBinContent(iphi, itheta, 2)
                        counter+=1
        print(counter)
        fout = ROOT.TFile("test_fancyPMT.root", "RECREATE")
        self.h.Write(self.h.GetName())

    def Draw(self, dopt):
        # self.h.SetMarkerStyle(1)
        self.h.SetLineColor(15)
        self.h.SetLineWidth(1)
        self.h.SetLineStyle(4)
        # self.h.SetMarkerSize(10)
        # self.h.SetMarkerColor(2)
        self.h.Draw(dopt)

## Jie: NOTE: run the code below first to create test.root, then run the __main__ below
if __name__=='__main__':
    h = TH2DPMTfancy('h','test', 500)
    # h = TH2DPMTfancy('h','test', 150)

    print "now processing %s\n"%(sys.argv[1])
    file0 = str(sys.argv[1])
    fileName = os.path.basename(file0)
    checkEventID = int(sys.argv[2])
    fin = TFile(file0)
    tree1 = fin.Get("T1")

    runIDval = array('l',[0])
    subrunIDval = array('l',[0])
    evtID = array('l',[0])
    mbx = array('f',[0])
    mby = array('f',[0])
    mbz = array('f',[0])
    
    tf2x = array('f',[0])
    tf2y = array('f',[0])
    tf2z = array('f',[0])
    
    mcx = array('f', [0])
    mcy = array('f', [0])
    mcz = array('f', [0])
    
    qpeVal = array('f',[0]) #unsigned double 
    fpromptVal = array('f',[0])
    nSCBayesVal = array('f',[0])
    rprompt60BayesVal = array('f',[0])
    fmaxpeVal = array('f',[0])
    nhitVal = array('i',[0])
    pulseGar = array('f',[0])
    cft2r = array('f',[0]) ## (chargetopring + chargesecondring)/qpe < 0.04
    cfb3r = array('f',[0]) ## (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qpe<0.1
    neckVetoVal = array('f',[0])
    nPMTs1 = array('I',[0])
    
    Nmax = 255
    pmtPhi = array('f',Nmax*[0])
    pmtCosTheta = array('f',Nmax*[0])
    pmtCharge = array('f',Nmax*[0])
    pmtChargeNSCB = array('f',Nmax*[0])
    pmtChargePrompt = array('f',Nmax*[0])
    listPmtID = array('I',Nmax*[0])
    myda = {}
    nentries = tree1.GetEntries()
    for event in range(nentries):
        if (event%20000 == 0): print "processed", event
        tree1.GetEntry(event)
        eventID = tree1.eventID
        if eventID != checkEventID: continue
        listPMTid = tree1.listPmtID
        mcx = tree1.mcx; mcy = tree1.mcy; mcz = tree1.mcz

        mbx = tree1.mbx; mby = tree1.mby; mbz = tree1.mbz

        tf2x = tree1.tf2x; tf2y = tree1.tf2y; tf2z = tree1.tf2z;
        
        posMC = TVector3(mcx, mcy, mcz)
        posMBL = TVector3(mbx, mby, mbz)
        posTF2 = TVector3(tf2x, tf2y, tf2z)
        pmtCharge = tree1.pmtCharge
        pmtChargeNSCB = tree1.pmtChargeNSCB
        pmtChargePrompt = tree1.pmtChargePrompt

        #print [x for x in pmtCharge]
        #print [x for x in pmtChargeNSCB]
        #print [x for x in pmtChargePrompt]
        #print "charge", pmtCharge[10]
        checkPMTids = [pmtid for pmtid in listPMTid]
        #for pmtid in listPMTid:
        #    index = checkPMTids.index(pmtid)
        #    #cosAngle_MC = pmtPos.Unit()*posMC.Unit()
        #    #cosAngle_MB = pmtPos.Unit()*posMBL.Unit()
        #    #cosAngle_TF2 = pmtPos.Unit()*posTF2.Unit()

        #    #print index, pmtid
        #    pmtqpe = pmtCharge[index]
        #    pmtnscb = pmtChargeNSCB[index]
        #    pmtprompt = pmtChargePrompt[index]
        #    
        #    #myda[pmtid] = pmtqpe
        #    myda[pmtid] = pmtnscb
        #    #myda[pmtid] = pmtprompt
        
        # ### loop all PMTids, also fill 0 charges!!
        for pmtid in range(255):
            if pmtid in listPMTid:
               index = checkPMTids.index(pmtid)
               #print index, pmtid
               pmtqpe = pmtCharge[index]
               pmtnscb = pmtChargeNSCB[index]
               pmtprompt = pmtChargePrompt[index]
               myda[pmtid] = pmtqpe
               #myda[pmtid] = pmtnscb
               #myda[pmtid] = pmtprompt
            else:
               myda[pmtid] = 0

    #myda = [x+1 for x in range(255)]### put my data here!!!
    #uncomment to load 1d histogram with 255 bins from input file
    #ipf = ROOT.TFile(ipfn)
    #h_inputdata = ipf.Get('histogram_name')
    #myda = [h_inputdata.GetBinContent(ibin+1) for ibin in range(255)]
    Canvas = ROOT.TCanvas ("Canvas", "Canvas", 1200, 1000);           Canvas.SetTickx(1);          Canvas.SetTicky(1);     Canvas.SetLogz(1);    Canvas.SetLogy(0);   Canvas.SetGrid(0,0)
    h.Fill(myda)
    h.Draw('CONT0');     Canvas.SaveAs("test0.png"); Canvas.Clear()
    #h.Draw('CONT1');     Canvas.SaveAs("test1.png"); Canvas.Clear()
    #h.Draw('CONT2');     Canvas.SaveAs("test2.png"); Canvas.Clear()
    #h.Draw('CONT3');     Canvas.SaveAs("test3.png"); Canvas.Clear()
    #h.Draw('CONT4');     Canvas.SaveAs("test4.png"); Canvas.Clear()
    # h.Draw('SURF2');     Canvas.SaveAs("test.png"); Canvas.Clear()
    # au = raw_input('>')
    Canvas.Close()

if __name__=='__main__':

    fin = ROOT.TFile("test_fancyPMT.root")
    h = fin.Get('h'); h.SetDirectory(0);
    fin.Close()
    h.SetLineColor(15)
    h.SetLineWidth(1)
    h.SetLineStyle(4)
    h.GetXaxis().CenterTitle(0)
    h.SetTitle('')
    h.Draw('CONT0')
    Canvas = ROOT.TCanvas ("Canvas", "Canvas", 1400, 1000);Canvas.SetLogz(1);Canvas.SetLogy(0);Canvas.SetGrid(0,0)
    Canvas.SetTickx(1);Canvas.SetTicky(0); 
    # for pmtid in range(255):
    # h.Draw('CONT0');     Canvas.SaveAs("test0.png"); Canvas.Clear()
    # h.Draw('CONT1');     Canvas.SaveAs("test1.png"); Canvas.Clear()
    # h.Draw('CONT2');     Canvas.SaveAs("test2.png"); Canvas.Clear()
    h.Draw('CONT0')
    ### draw PMT ids
    for pmtid in range(255):
        x,y,z = pp.pmtvec[pmtid][0], pp.pmtvec[pmtid][1], pp.pmtvec[pmtid][2]
        Tl = ROOT.TLatex ()
        Tl.SetTextAlign(22)
        Tl.SetTextSize(0.02)

        phi  = np.arctan2(y, x) # Jie: switch to arctan2
        #phi = np.arctan(y/x)
        #if (x<0 and y<0):
        #    phi -=np.pi
        #if (x<0 and y>0):
        #    phi +=np.pi
        costheta = (z/np.sqrt(x**2+y**2+z**2)).tolist()

        if phi>3.10:
            phi=3.10
        if phi<-3.10:
            phi = -3.10
        # Tl.DrawLatex(phi, costheta, '#color[4]{'+str(pmtid)+"}")
        Tl.DrawLatex(phi, costheta, str(pmtid))
    ### draw gp ids
    phis, costhetas = [-1.9, -0.63, 0.62, 1.9, 3.10, -2.5, -1.25, 0, 1.25, 2.5, -0.5], [0.45]*5+[-0.45]*5 + [-1.05]
    for i in range(11):
        phi, costheta = phis[i], costhetas[i]
        Tl = ROOT.TLatex()
        Tl.SetTextAlign(22)
        Tl.SetTextSize(0.04)
        Tl.DrawLatex(phi, costheta, '#color[2]{'+str(i)+'}')

    Canvas.SaveAs("PMTs_"+str(checkEventID)+".png");
    raw_input()
    #Canvas.Clear()
    #Canvas.Close()