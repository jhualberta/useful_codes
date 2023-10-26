#from ratproc.base import Processor
import sys
import os
import operator
import glob
#from rat import *
#from rat import PMTInfoUtil
from ROOT import *
import numpy as np
from numpy import * 
import math
from array import array
from scipy.special import gammaln

pmtpos = './Aksel_PMTpos_260.ratdb'
if not os.path.isfile(pmtpos):
    print 'ERROR: need Aksel_PMTpos_260.ratdb to work'
    print 'INFO: inside TH2DPMTfancy.py, edit the variable pmtpos to the right place'
    sys.exit(1)

filePDF = TFile("surface_tagging_pdfs.root")

print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
fileName = os.path.basename(file0)

fout = TFile("dump_"+fileName,"recreate")

# load all the required PDFs
## MC position for alpha surface and ar40 PDF
Hsurf_truth_angle = filePDF.Get("surf_truth_angle_qPE")
Har40_truth_angle = filePDF.Get("ar40_truth_angle_qPE")
### !!!! updated PDF
#fileNew = TFile("saveQPE_checkHistosOFF_MBR800.root")
#Har40_truth_angle = fileNew.Get("HpmtID_costhetaMC")

### dead pmts
Hsurf_PMTs_angle  = filePDF.Get("surf_PMTs_angle_qPE")
Hsurf_gap10_angle = filePDF.Get("surf_gap10_angle_qPE")
Hsurf_gaps_angle  = filePDF.Get("surf_gaps_angle_qPE")

### build new PDFs
Hcostheta_TF2 = TH1F("Hcostheta_TF2","costheta TF2",255,-1,1)
Hcostheta_MBL = TH1F("Hcostheta_MBL","costheta MBL", 255,-1,1)
Hcostheta_MC  = TH1F("Hcostheta_MC","costheta MC truth", 255,-1,1) 
Hcostheta_MC_scale  = TH1F("Hcostheta_MC_scale","costheta MC truth, scaled by dead PMTs", 255,-1,1) 

### building actual PDFs, ar40_truth_angle
Hsort_costheta_MC = TH1F("Hsort_costheta_MC","sorted costheta MC truth", 255,0,255)
Hsort_costheta_MBL = TH1F("Hsort_costheta_MBL","sorted costheta MBL truth", 255,0,255)
Hsort_costheta_TF2 = TH1F("Hsort_costheta_TF2","sorted costheta TF2 truth", 255,0,255)
Hsort_costheta_MC_scale = TH1F("Hsort_costheta_MC_scale","sorted costheta MC truth, scaled by dead PMTs", 255,0,255)
Hcostheta_MC_badPMT149 = TH1F("Hcostheta_MC_badPMT149","sorted costheta MC truth, iPMT == 149", 255,0,255)
Hcostheta_MC_badPMT204 = TH1F("Hcostheta_MC_badPMT204","sorted costheta MC truth, iPMT == 204", 255,0,255)
#### checkings
HpmtID_pmtCharges = TH1F("HpmtID_pmtCharges", "pmtID vs pmtCharges", 255,0,255)
Hsort_costheta_ar40MBL  = TH1F("Hsort_costheta_ar40MBL", "sorted costheta ar40 MBL", 255, 0, 255)
HpmtID_costheta_ar40MBL = TH1F("HpmtID_costheta_ar40MBL", "pmtID costheta ar40 MBL", 255, 0, 255)
HpmtID_pdf_ar40MBL      = TH1F("HpmtID_pdf_ar40MBL", "pmtID qpe*PDF ar40 MBL", 255, 0, 255)
Hsort_costheta_surfMBL  = TH1F("Hsort_costheta_surfMBL", "sorted costheta surf MBL", 255, 0, 255)
HpmtID_costheta_surfMBL = TH1F("HpmtID_costheta_surfMBL", "pmtID costheta surf MBL", 255, 0, 255)
HpmtID_pdf_surfMBL = TH1F("HpmtID_pdf_surfMBL", "pmtID qpe*PDF surf MBL", 255, 0, 255)

Hsort_costheta_ar40TF2  = TH1F("Hsort_costheta_ar40TF2", "sorted costheta ar40 TF2", 255, 0, 255)
HpmtID_costheta_ar40TF2 = TH1F("HpmtID_costheta_ar40TF2", "pmtID costheta ar40 TF2", 255, 0, 255)
HpmtID_pdf_ar40TF2      = TH1F("HpmtID_pdf_ar40TF2", "pmtID qpe*PDF ar40 TF2", 255, 0, 255)
Hsort_costheta_surfTF2  = TH1F("Hsort_costheta_surfTF2", "sorted costheta surf TF2", 255, 0, 255)
HpmtID_costheta_surfTF2 = TH1F("HpmtID_costheta_surfTF2", "pmtID costheta surf TF2", 255, 0, 255)
HpmtID_pdf_surfTF2 = TH1F("HpmtID_pdf_surfTF2", "pmtID qpe*PDF surf TF2", 255, 0, 255)

Hsave_pmtsPDF = TH1F("Hsave_pmtsPDF","min(NLL_pmtsPDF)", 2000,-1000,1000)
Hsave_gapsPDF = TH1F("Hsave_gapsPDF","min(NLL_gapsPDF)", 2000,-1000,1000)

Hsave_ar40PDF_MBL = TH1F("Hsave_ar40PDF_MBL","ar40PDF",2000,-1000,1000)
Hsave_surfPDF_MBL = TH1F("Hsave_surfPDF_MBL","surfPDF",2000,-1000,1000)
Hsave_ar40_min_MBL = TH1F("Hsave_ar40_min_MBL", "ar40_min_MBL", 2000,-1000,1000)

Hsave_ar40PDF_TF2 = TH1F("Hsave_ar40PDF_TF2","ar40PDF",2000,-1000,1000)
Hsave_surfPDF_TF2 = TH1F("Hsave_surfPDF_TF2","surfPDF",2000,-1000,1000)
Hsave_ar40_min_TF2 = TH1F("Hsave_ar40min_TF2","ar40_min_TF2",2000,-1000,1000) 

# load offline pmt positions
ff = file(pmtpos,'r')
pmtpos_offline = []
xpos = []
ypos = []
zpos = []

dead_pmt_id = [149,204]
pmt_radius = 851
for l in ff:
    if l.startswith('x'): xpos = [float(val) for val in l.split('[')[1][:-4].split(',')]
    if l.startswith('y'): ypos = [float(val) for val in l.split('[')[1][:-4].split(',')]
    if l.startswith('z'): zpos = [float(val) for val in l.split('[')[1][:-4].split(',')]

    for i,(x,y,z) in enumerate(zip(xpos, ypos, zpos)):
        v = TVector3(x,y,z)
        pmtpos_offline.append(pmt_radius*v.Unit())

## 2 bad pmts
#bad_pmtpos_offline = []
#for pmtid in dead_pmt_id:
#   bad_pmtpos_offline.append(pmtpos_offline[pmtid])

#############################################################
## dump likelihoods
checkEventID = 37
##Fix for 2 dead pmts, 11 gaps
HDump_NLL = TH1F("HDump_NLL","",255,0,255)
ListHDump_NLL_deadPMT = []
ListHDump_Q_deadPMT = []
ListHDump_pdf_deadPMT = []
ListHDumpPMTid_NLL_deadPMT = []
ListHDumpPMTid_Q_deadPMT = []
ListHDumpPMTid_pdf_deadPMT = []

for i in range(2):
    htemp = HDump_NLL.Clone()
    htemp.SetName("HDump_NLL_deadPMT"+str( dead_pmt_id[i] ))
    ListHDump_NLL_deadPMT.append(htemp)
    htemp2 = HDump_NLL.Clone()
    htemp2.SetName("HDump_Q_deadPMT"+str( dead_pmt_id[i] ))
    ListHDump_Q_deadPMT.append(htemp2)
    htemp3 = HDump_NLL.Clone()
    htemp3.SetName("HDump_pdf_deadPMT"+str( dead_pmt_id[i] ))
    ListHDump_pdf_deadPMT.append(htemp3)
    ## not by sorted pmt index 
    htemp1a = HDump_NLL.Clone()
    htemp1a.SetName("HDumpPMTid_NLL_deadPMT"+str( dead_pmt_id[i] ))
    ListHDumpPMTid_NLL_deadPMT.append(htemp1a)
    htemp2a = HDump_NLL.Clone()
    htemp2a.SetName("HDumpPMTid_Q_deadPMT"+str( dead_pmt_id[i] ))
    ListHDumpPMTid_Q_deadPMT.append(htemp2a)
    htemp3a = HDump_NLL.Clone()
    htemp3a.SetName("HDumpPMTid_pdf_deadPMT"+str( dead_pmt_id[i] ))
    ListHDumpPMTid_pdf_deadPMT.append(htemp3a)

ListHDump_NLL_gaps = []
ListHDump_Q_gaps = []
ListHDump_pdf_gaps = []
ListHDumpPMTid_NLL_gaps = []
ListHDumpPMTid_Q_gaps = []
ListHDumpPMTid_pdf_gaps = []

for i in range(11):
    htemp = HDump_NLL.Clone()
    htemp.SetName("HDump_NLL_gap"+str(i))
    ListHDump_NLL_gaps.append(htemp)
    htemp2 = HDump_NLL.Clone()
    htemp2.SetName("HDump_Q_gap"+str(i))
    ListHDump_Q_gaps.append(htemp2)
    htemp3 = HDump_NLL.Clone()
    htemp3.SetName("HDump_pdf_gap"+str(i))
    ListHDump_pdf_gaps.append(htemp3)
    ## not by sorted pmt index 
    htemp1a = HDump_NLL.Clone()
    htemp1a.SetName("HDumpPMTid_NLL_gap"+str(i))
    ListHDumpPMTid_NLL_gaps.append(htemp1a)
    htemp2a = HDump_NLL.Clone()
    htemp2a.SetName("HDumpPMTid_Q_gap"+str(i))
    ListHDumpPMTid_Q_gaps.append(htemp2a)
    htemp3a = HDump_NLL.Clone()
    htemp3a.SetName("HDumpPMTid_pdf_gap"+str(i))
    ListHDumpPMTid_pdf_gaps.append(htemp3a)

def pmts_gaps_ids():
    # construct arrays of the pentagonal gaps --> from SurfaceSearchProc.cc
    gaps_arrays = []
    gaps_arrays.append( pmtpos_offline[50] + pmtpos_offline[49] + pmtpos_offline[75] + pmtpos_offline[74] + pmtpos_offline[92] )
    gaps_arrays.append( pmtpos_offline[48] + pmtpos_offline[47] + pmtpos_offline[73] + pmtpos_offline[72] + pmtpos_offline[91] )
    gaps_arrays.append( pmtpos_offline[46] + pmtpos_offline[45] + pmtpos_offline[71] + pmtpos_offline[70] + pmtpos_offline[90] )
    gaps_arrays.append( pmtpos_offline[54] + pmtpos_offline[53] + pmtpos_offline[79] + pmtpos_offline[78] + pmtpos_offline[94] )
    gaps_arrays.append( pmtpos_offline[52] + pmtpos_offline[77] + pmtpos_offline[93] + pmtpos_offline[51] + pmtpos_offline[76] )

    gaps_arrays.append( pmtpos_offline[157] + pmtpos_offline[176] + pmtpos_offline[175] + pmtpos_offline[201] + pmtpos_offline[200] )
    gaps_arrays.append( pmtpos_offline[156] + pmtpos_offline[174] + pmtpos_offline[173] + pmtpos_offline[199] + pmtpos_offline[198] )
    gaps_arrays.append( pmtpos_offline[155] + pmtpos_offline[172] + pmtpos_offline[171] + pmtpos_offline[197] + pmtpos_offline[196] )
    gaps_arrays.append( pmtpos_offline[159] + pmtpos_offline[170] + pmtpos_offline[179] + pmtpos_offline[195] + pmtpos_offline[204] )
    gaps_arrays.append( pmtpos_offline[158] + pmtpos_offline[178] + pmtpos_offline[177] + pmtpos_offline[203] + pmtpos_offline[202] )

    # construct tvectors of the pentagonal gaps
    gaps_vectors=[TVector3() for i in range(11)]
    for i in range(10):
        gaps_vectors[i].SetXYZ(gaps_arrays[i][0], gaps_arrays[i][1], gaps_arrays[i][2])
        gaps_vectors[i].SetMag(851.0)
        if i==4:
            gaps_vectors[4].SetPhi(np.pi)
    gaps_vectors[10].SetXYZ(0.0,0.0,-851.0)

    # check if there are any dead/bad PMTs in this run!
    # if there are, construct their tvector
    dead_pmts = dead_pmt_id ## 149, 204
    dead_pmts_vectors = []
    if len(dead_pmts):
        dead_pmts_vectors = [TVector3() for i in range(len(dead_pmts))]
        for i in range(len(dead_pmts)):
            dead_pmts_vectors[i].SetXYZ(pmtpos_offline[dead_pmts[i]][0], pmtpos_offline[dead_pmts[i]][1], pmtpos_offline[dead_pmts[i]][2])
            dead_pmts_vectors[i].SetMag(851)
    pmts_angles = [{} for i in range(len(dead_pmts_vectors))]
    gaps_angles = [{} for i in range(len(gaps_vectors))]

    temp_pmt_pos = TVector3()
    pmts_angles_ids = []
    for iPMT in range(255):
        temp_pmt_pos.SetXYZ( pmtpos_offline[iPMT][0], pmtpos_offline[iPMT][1], pmtpos_offline[iPMT][2] )
        temp_pmt_pos.SetMag(851)

        for i in range(len(dead_pmts_vectors)):
            pmts_angles[i][iPMT] = np.cos(dead_pmts_vectors[i].Angle((temp_pmt_pos)))
        for i in range(len(gaps_vectors)):
            gaps_angles[i][iPMT] = np.cos(gaps_vectors[i].Angle((temp_pmt_pos)))

    for i in range(len(pmts_angles)):
        pmts_angles_dict = pmts_angles[i]
        pmts_angles_key = sorted(pmts_angles_dict.items(), key=operator.itemgetter(1), reverse=True)
        pmts_angles_ids.append( [item[0] for item in pmts_angles_key] )

    gaps_angles_ids = []
    for i in range(len(gaps_angles)):
        gaps_angles_dict = gaps_angles[i]
        gaps_angles_key  = sorted(gaps_angles_dict.items(), key=operator.itemgetter(1), reverse=True)
        gaps_angles_ids.append( [item[0] for item in gaps_angles_key] )

    del temp_pmt_pos
    del gaps_angles, gaps_angles_dict, gaps_angles_key, gaps_vectors
    if len(dead_pmts):
        del pmts_angles, pmts_angles_dict, pmts_angles_key, dead_pmts_vectors

    return pmts_angles_ids, gaps_angles_ids


####################################################################################################################################################
def get_poisson_nll(xbar, int_pmt_charge):
    if int_pmt_charge<170:
        #return -(int_pmt_charge*np.log(xbar) - xbar - np.log(np.double(math.factorial(int_pmt_charge))))
        return -( int_pmt_charge*math.log(xbar) - xbar - math.log(np.double(math.factorial(int_pmt_charge))) )
    ## Jie: there is a bug in python2 if int_pmt_charge is too large, get "for long object has no attribute 'log'"; should use math.log
    # see https://stackoverflow.com/questions/66404748/how-to-compute-log-factorial-of-an-array-of-numbers
    else:
        #return -(int_pmt_charge*np.log(xbar) - xbar - gammaln(int_pmt_charge))
        return -( int_pmt_charge*np.log(xbar) - xbar - (gammaln(int_pmt_charge) + log(int_pmt_charge)) )
    ## Jie: gammaln = log(gamma(x)) = log((x-1)!), so should be gammaln(x) + log(x) = log(x!)
    ## or why don't we just use gammaln(x), and x is not necessary to be integer
####################################################################################################################################################
def gof_calculator(xbar, pmt_charge):
    return (xbar - pmt_charge)**2/xbar

# Also sort the PMT indices since they are independent of the fitted event position
pmts_angles_ids, gaps_angles_ids = pmts_gaps_ids()
#-------------------------------------------------------------------------------------------------------------------------------------------------------
# to hold MBL and TF2 fitted positions
MBLEvent = TVector3()
TF2Event = TVector3()
MCEvent = TVector3()
########################################################################################################################################################
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
pmtChargeNSCB = array('f', Nmax*[0])
pmtChargePrompt = array('f', Nmax*[0])
listPmtID = array('I',Nmax*[0])

File = TFile(file0);
tree1 = File.Get("T1");#TTree
nentries = tree1.GetEntries();
tree1.SetBranchAddress("eventID",evtID) #ULong64_
tree1.SetBranchAddress("qpe", qpeVal)
tree1.SetBranchAddress("fprompt", fpromptVal)
tree1.SetBranchAddress("nSCBayes", nSCBayesVal)
tree1.SetBranchAddress("rprompt60Bayes", rprompt60BayesVal)
tree1.SetBranchAddress("mbx", mbx)
tree1.SetBranchAddress("mby", mby)
tree1.SetBranchAddress("mbz", mbz)

tree1.SetBranchAddress("mcx", mcx)
tree1.SetBranchAddress("mcy", mcy)
tree1.SetBranchAddress("mcz", mcz)

tree1.SetBranchAddress("tf2x", tf2x)
tree1.SetBranchAddress("tf2y", tf2y)
tree1.SetBranchAddress("tf2z", tf2z)
tree1.SetBranchAddress("nPMTs1",nPMTs1)
tree1.SetBranchAddress("pmtPhi", pmtPhi)
tree1.SetBranchAddress("pmtCosTheta", pmtCosTheta)
tree1.SetBranchAddress("pmtCharge", pmtCharge)
tree1.SetBranchAddress("pmtChargeNSCB", pmtChargeNSCB)
tree1.SetBranchAddress("pmtChargePrompt", pmtChargePrompt)
tree1.SetBranchAddress("listPmtID", listPmtID)

#ntupleSurfTag = TNtuple("ntupleSurfTag","surface tag variables","eventID:NLL_surfPDF_MBL:NLL_ar40PDF_MBL:NLL_surfPDF_TF2:NLL_ar40PDF_TF2:min(NLL_pmtsPDF):min(NLL_gapsPDF):gof_surfPDF_MBL:gof_ar40PDF_MBL:gof_surfPDF_TF2:gof_ar40PDF_TF2")

nentries = 5000
for event in range(nentries):
    if (event%20000 == 0): print "processed", event
    tree1.GetEntry(event)
    eventID = evtID[0]
    xmc = mcx[0]
    ymc = mcy[0]
    zmc = mcz[0]

    evtx = mbx[0]
    evty = mby[0]
    evtz = mbz[0]
    evtx_tf2 = tf2x[0]
    evty_tf2 = tf2y[0]
    evtz_tf2 = tf2z[0]

    if sqrt(evtx*evtx+evty*evty+evtz*evtz)>800: continue

    MBLEvent.SetXYZ(evtx, evty, evtz)
    TF2Event.SetXYZ(evtx_tf2, evty_tf2, evtz_tf2)
    
    MCEvent.SetXYZ(xmc, ymc, zmc)

    TF2_fit_found = 1 if TF2Event.Mag()<855 else 0
    MBL_fit_found = 1 if MBLEvent.Mag()<855 else 0
    if (not TF2_fit_found) and (not MBL_fit_found): continue
    #------------------------------------------------------------------------
    # sort PMT IDs!
    # NOTE 1: The sorting that is done within the *pmts_gaps_ids* function, sort the PMTs based on the pentagonal gap positions
    # and dead PMTS (if they exist)! This block here sorts the PMTs based on the MBL and TF2 positions.
    # NOTE 2: The cos theta is the angle between two vectors v1 and v2_i where v1 is from the centre of the detector to the event position and v2_i is from the centre to PMT_i.
    # Pentagonal gaps and dead PMTs are similar, where v1 is the vector from the centre to the position of the gap/dead PMT

    MBL_angles, TF2_angles = {}, {}
    MC_angles = {} ## only for truth position MC!!
    # first we need all the 255 PMTs to sort them based on cos_theta!
    temp_pmt_pos = TVector3()
    for iPMT in range(255):
        temp_pmt_pos.SetXYZ( pmtpos_offline[iPMT][0], pmtpos_offline[iPMT][1], pmtpos_offline[iPMT][2] )
        temp_pmt_pos.SetMag(851)

        MBL_angles[iPMT] = np.cos(MBLEvent.Angle((temp_pmt_pos)))
        TF2_angles[iPMT] = np.cos(TF2Event.Angle((temp_pmt_pos)))

        ### NOTE: ONLY for MC!!!!
        MC_angles[iPMT] = np.cos(MCEvent.Angle((temp_pmt_pos)))

    MBL_angles_dict = sorted(MBL_angles.items(), key=operator.itemgetter(1), reverse=True)
    MBL_angles_ids  = [item[0] for item in MBL_angles_dict]

    TF2_angles_dict = sorted(TF2_angles.items(), key=operator.itemgetter(1), reverse=True)
    TF2_angles_ids  = [item[0] for item in TF2_angles_dict]

    ###NOTE: ONLY for MC!!!!
    MC_angles_dict = sorted(MC_angles.items(), key=operator.itemgetter(1), reverse=True)
    MC_angles_ids  = [item[0] for item in MC_angles_dict]

    #------------------------------------------------------------------------
    # initialize NLL and gof variables
    NLL_surfPDF_MBL, NLL_ar40PDF_MBL = 0, 0 #use MBL position and the surface and ar40 PDF
    NLL_surfPDF_TF2, NLL_ar40PDF_TF2 = 0, 0 #use TF2 position and the surface and ar40 PDF

    gof_surfPDF_MBL , gof_ar40PDF_MBL  = 0, 0 #use MBL position and the surface and ar40 PDF
    gof_surfPDF_TF2 , gof_ar40PDF_TF2  = 0, 0 #use TF2 position and the surface and ar40 PDF

    NLL_pmtsPDF = [0 for i in range(len(pmts_angles_ids))] #use the dead PMT positions and the PMT PDFs
    NLL_gapsPDF = [0 for i in range(len(gaps_angles_ids))] #use the pentagonal gaps positions and the Gap PDFs

    gof_pmtsPDF_PMTs = [0 for i in range(len(pmts_angles_ids))] #use the dead PMT positions and the PMT PDFs
    gof_gapsPDF_gaps = [0 for i in range(len(gaps_angles_ids))] #use the pentagonal gaps positions and the Gap PDFs

    pmt_ids = []
    ### Warning !!! Here we consider the total charge of event, event qPE or nSCBayes
    qPE = qpeVal[0]
    nPMT = nPMTs1[0]
    list_255charges = []
    for iPMT in range(255):
        if iPMT==nPMT:
            # Complete the list of PMT IDs with the PMTs that did not see any charge in this event.
            # Since we are using Poisson distribution, PMTs not seeing charge can contain information for the likelihood.
            pmt_ids += [i for i in range(255) if i not in pmt_ids]

        if iPMT<nPMT:
            # NOTE: NLL calculation needs integer expected charge; but not the gof calculation.
            #pmt = cal.GetPMT(iPMT)
            #pmt_id = pmt.GetID()
            #pmt_qPE = pmt.qPE
            pmt_id = listPmtID[iPMT]
            ##### change charges !!!
            pmt_qPE = pmtCharge[iPMT]
            #pmt_qPE = pmtChargeNSCB[iPMT]
            list_255charges.append(pmt_qPE)
            pmt_ids.append(pmt_id)
        else:
            # PMTs that have not observed any charge.
            pmt_id = pmt_ids[iPMT]
            pmt_qPE = 0
            list_255charges.append(pmt_qPE)

        ##NOTE: only for MC!!!!
        sortedAngle = MC_angles_dict[MC_angles_ids.index(pmt_id)][1] ### costheta sorted by pmt.qPE, using pmt_id
        Hcostheta_MC.Fill(sortedAngle, pmt_qPE)
        HpmtID_pmtCharges.Fill(pmt_id, pmt_qPE)

        if iPMT == 149:
            Hcostheta_MC_badPMT149.Fill(sortedAngle, pmt_qPE)
        if iPMT == 204:
            Hcostheta_MC_badPMT204.Fill(sortedAngle, pmt_qPE)
        #Hsort_costheta_MC.Fill()
        #+1 is needed for GetBinContent because index 0 is underflow.
        #------------------------------------------------------------------------
        # calculate NLL and gof based on MBL position
        if MBL_fit_found:
            pdf_index = MBL_angles_ids.index(pmt_id)+1
            xbar_surfPDF = Hsurf_truth_angle.GetBinContent(pdf_index)*qPE
            xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
            nll_surf_MBL = get_poisson_nll(xbar_surfPDF, int(round(pmt_qPE)))
            NLL_surfPDF_MBL += nll_surf_MBL
            nll_ar40_MBL = get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))
            NLL_ar40PDF_MBL += nll_ar40_MBL
            
            gof_surfPDF_MBL += gof_calculator (xbar_surfPDF, pmt_qPE)
            gof_ar40PDF_MBL += gof_calculator (xbar_ar40PDF, pmt_qPE)

            Hcostheta_MBL.Fill(MBL_angles_dict[MBL_angles_ids.index(pmt_id)][1], pmt_qPE)

            Hsort_costheta_ar40MBL.Fill(pdf_index-1, nll_ar40_MBL)
            HpmtID_costheta_ar40MBL.Fill(pmt_id, nll_ar40_MBL)
            HpmtID_pdf_ar40MBL.Fill(pmt_id, xbar_ar40PDF)
            Hsort_costheta_surfMBL.Fill(pdf_index-1, nll_surf_MBL)
            HpmtID_costheta_surfMBL.Fill(pmt_id, nll_surf_MBL)
            HpmtID_pdf_surfMBL.Fill(pmt_id, xbar_surfPDF)
 
        #------------------------------------------------------------------------
        # calculate NLL and gof based on TF2 position
        if TF2_fit_found:
            pdf_index = TF2_angles_ids.index(pmt_id)+1
            xbar_surfPDF = Hsurf_truth_angle.GetBinContent(pdf_index)*qPE
            xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
            nll_surf_TF2 = get_poisson_nll(xbar_surfPDF, int(round(pmt_qPE)))
            NLL_surfPDF_TF2 += nll_surf_TF2
            nll_ar40_TF2 = get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))
            NLL_ar40PDF_TF2 += nll_ar40_TF2

            gof_surfPDF_TF2 += gof_calculator(xbar_surfPDF, pmt_qPE)
            gof_ar40PDF_TF2 += gof_calculator(xbar_ar40PDF, pmt_qPE)

            Hcostheta_TF2.Fill(TF2_angles_dict[TF2_angles_ids.index(pmt_id)][1], pmt_qPE)

            Hsort_costheta_ar40TF2.Fill(pdf_index-1, nll_ar40_TF2)
            HpmtID_costheta_ar40TF2.Fill(pmt_id, nll_ar40_TF2)
            HpmtID_pdf_ar40TF2.Fill(pmt_id, xbar_ar40PDF)
            Hsort_costheta_surfTF2.Fill(pdf_index-1, nll_surf_TF2)
            HpmtID_costheta_surfTF2.Fill(pmt_id, nll_surf_TF2)
            HpmtID_pdf_surfTF2.Fill(pmt_id, xbar_surfPDF)
        #------------------------------------------------------------------------
        # calculate NLL and gof based on dead PMT position(s)
        for i in range(len(pmts_angles_ids)): ##2 dead PMTs
            pdf_index = pmts_angles_ids[i].index(pmt_id)+1
            overwrite = True if pmt_id in dead_pmt_id else False # avoid X_deapPMT*X_deadPMT
            if not overwrite: ## pmt_id not in dead pmts
                xbar_deadPDF = Hsurf_PMTs_angle.GetBinContent(pdf_index)*qPE
            else:
                # If a (second) dead PMT, expection would be also zero charge.
                xbar_deadPDF = Hsurf_PMTs_angle.GetBinContent(1)*qPE
            ## here pmt_qPE is looped!!
            nllVal = get_poisson_nll(xbar_deadPDF, int(round(pmt_qPE)))
            NLL_pmtsPDF[i] += nllVal 
            gof_pmtsPDF_PMTs[i] += gof_calculator(xbar_deadPDF, pmt_qPE)
            #if eventID == checkEventID: continue
            ## sorted PMT index
            ListHDump_NLL_deadPMT[i].Fill(pdf_index-1, nllVal)
            ListHDump_Q_deadPMT[i].Fill(pdf_index-1, pmt_qPE)
            ListHDump_pdf_deadPMT[i].Fill(pdf_index-1, xbar_deadPDF)
            ## PMT ids
            ListHDumpPMTid_NLL_deadPMT[i].Fill(pmt_id, nllVal)
            ListHDumpPMTid_Q_deadPMT[i].Fill(pmt_id, pmt_qPE)
            ListHDumpPMTid_pdf_deadPMT[i].Fill(pmt_id, xbar_deadPDF)
        #------------------------------------------------------------------------
        # calculate NLL and gof based on pentagonal gap positions
        for i in range(len(gaps_angles_ids)):
            pdf_index = gaps_angles_ids[i].index(pmt_id)+1 ## iPMT is sorted to be pdf_index
            xbar_gapsPDF = Hsurf_gaps_angle.GetBinContent(pdf_index)*qPE
            nllVal = get_poisson_nll(xbar_gapsPDF, int(round(pmt_qPE)))
            NLL_gapsPDF[i] += nllVal
            gof_gapsPDF_gaps[i] += gof_calculator(xbar_gapsPDF, pmt_qPE)
            #if eventID != checkEventID: continue
            ## sorted PMT index 
            ListHDump_NLL_gaps[i].Fill(pdf_index-1, nllVal)
            ListHDump_Q_gaps[i].Fill(pdf_index-1, pmt_qPE)
            ListHDump_pdf_gaps[i].Fill(pdf_index-1, xbar_gapsPDF)
            ## PMT ids
            ListHDumpPMTid_NLL_gaps[i].Fill(pmt_id, nllVal)
            ListHDumpPMTid_Q_gaps[i].Fill(pmt_id, pmt_qPE)
            ListHDumpPMTid_pdf_gaps[i].Fill(pmt_id, xbar_gapsPDF)

    #------------------------------------------------------------------------

    # safety checks
    if not len(pmts_angles_ids):
        # to avoid min arg exception if we don't have any dead PMTs
        NLL_pmtsPDF.append(1e6)
        gof_pmtsPDF_PMTs.append(1e6)

    if not TF2_fit_found:
        NLL_surfPDF_TF2 = 1e6
        NLL_ar40PDF_TF2 = 1e6
    if not MBL_fit_found:
        NLL_surfPDF_MBL = 1e6
        NLL_ar40PDF_MBL = 1e6

    #print "nll dead pmtsPDF", NLL_pmtsPDF, gof_pmtsPDF_PMTs

    #------------------------------------------------------------------------
    # Set main NLL variables
    NLL_tag_MBL_list = [round(NLL_ar40PDF_MBL,3), round(NLL_surfPDF_MBL,3), round(min(NLL_pmtsPDF),3), round(min(NLL_gapsPDF),3)]
    NLL_tag_TF2_list = [round(NLL_ar40PDF_TF2,3), round(NLL_surfPDF_TF2,3), round(min(NLL_pmtsPDF),3), round(min(NLL_gapsPDF),3)]

    is_surface_like_MBL = 0 if min(NLL_tag_MBL_list)==NLL_tag_MBL_list[0] else 1
    is_surface_like_TF2 = 0 if min(NLL_tag_TF2_list)==NLL_tag_TF2_list[0] else 1

    # Mainly for events tag as ar40 --> How far/close was the next minimum NLL?!

    NLL_ar40_min_MBL = NLL_ar40PDF_MBL - min(NLL_tag_MBL_list[1:])
    NLL_ar40_min_TF2 = NLL_ar40PDF_TF2 - min(NLL_tag_TF2_list[1:])

    Hsave_pmtsPDF.Fill(min(NLL_pmtsPDF))
    Hsave_gapsPDF.Fill(min(NLL_gapsPDF))
    
    Hsave_ar40PDF_MBL.Fill(NLL_ar40PDF_MBL)
    Hsave_surfPDF_MBL.Fill(NLL_surfPDF_MBL)
    Hsave_ar40_min_MBL.Fill(NLL_ar40_min_MBL)
    
    Hsave_ar40PDF_TF2.Fill(NLL_ar40PDF_TF2)
    Hsave_surfPDF_TF2.Fill(NLL_surfPDF_TF2)
    Hsave_ar40_min_TF2.Fill(NLL_ar40_min_TF2)

    #------------------------------------------------------------------------
    # Set gof variables and the the two is_surface tags (one for MBL, one for TF2)
    #ev.Setgof_pmtsPDF    (min(gof_pmtsPDF_PMTs))
    #ev.Setgof_gapsPDF    (min(gof_gapsPDF_gaps))

    #ev.SetIs_surface_like_MBL(is_surface_like_MBL)
    #ev.SetIs_surface_like_TF2(is_surface_like_TF2)
    #print NLL_ar40_min_MBL, NLL_ar40_min_TF2

    #print "NLL_ar40_min MBL", round(NLL_ar40PDF_MBL, 3), "-", NLL_tag_MBL_list[1:], "=", round(NLL_ar40_min_MBL,3)
    #print "NLL_ar40_min_TF2", round(NLL_ar40PDF_TF2, 3), "-", NLL_tag_TF2_list[1:], "=", round(NLL_ar40_min_TF2,3)

    #ntupleSurfTag.Fill(eventID:NLL_surfPDF_MBL:NLL_ar40PDF_MBL:NLL_surfPDF_TF2:NLL_ar40PDF_TF2:min(NLL_pmtsPDF):min(NLL_gapsPDF):gof_surfPDF_MBL:gof_ar40PDF_MBL:gof_surfPDF_TF2:gof_ar40PDF_TF2)
    for i in range(255):
         Hsort_costheta_MC.SetBinContent(i+1, Hcostheta_MC.GetBinContent(255-i))
         Hsort_costheta_MBL.SetBinContent(i+1, Hcostheta_MBL.GetBinContent(255-i))
         Hsort_costheta_TF2.SetBinContent(i+1, Hcostheta_TF2.GetBinContent(255-i))

fout.cd()
for hist in ListHDump_NLL_gaps:
    hist.Write()

for hist in ListHDump_NLL_deadPMT:
    hist.Write()

for hist in ListHDump_Q_deadPMT:
    hist.Write()

for hist in ListHDump_pdf_deadPMT:
    hist.Write()


for hist in ListHDumpPMTid_NLL_gaps:
    hist.Write()

for hist in ListHDumpPMTid_Q_gaps:
    hist.Write()

for hist in ListHDumpPMTid_pdf_gaps:
    hist.Write()

for hist in ListHDumpPMTid_NLL_deadPMT:
    hist.Write()

for hist in ListHDumpPMTid_pdf_deadPMT:
    hist.Write()


Hsort_costheta_MC.Scale( 1/Hsort_costheta_MC.Integral() )
Hsort_costheta_MBL.Scale( 1/Hsort_costheta_MBL.Integral() )
Hsort_costheta_TF2.Scale( 1/Hsort_costheta_TF2.Integral() )

Hcostheta_TF2.Write()
Hcostheta_MBL.Write()
Hcostheta_MC.Write()
Hsort_costheta_MC.Write()
Hsort_costheta_MBL.Write()
Hsort_costheta_TF2.Write()
Hcostheta_MC_badPMT149.Write()
Hcostheta_MC_badPMT204.Write()

Hsort_costheta_ar40MBL.Scale( 1/Hsort_costheta_ar40MBL.Integral() )
Hsort_costheta_surfMBL.Scale( 1/Hsort_costheta_surfMBL.Integral() )
HpmtID_costheta_ar40MBL.Write()
Hsort_costheta_surfMBL.Write()
HpmtID_costheta_surfMBL.Write()

Hsort_costheta_ar40TF2.Scale( 1/Hsort_costheta_ar40TF2.Integral() )
Hsort_costheta_surfTF2.Scale( 1/Hsort_costheta_surfTF2.Integral() )
Hsort_costheta_ar40TF2.Write()
HpmtID_costheta_ar40TF2.Write()
Hsort_costheta_surfTF2.Write()
HpmtID_costheta_surfTF2.Write()
HpmtID_pdf_surfTF2.Write()
HpmtID_pdf_surfMBL.Write()

HpmtID_pmtCharges.Write()
## check NLL values
Hsave_pmtsPDF.Write()
Hsave_gapsPDF.Write()

Hsave_ar40PDF_MBL.Write()
Hsave_surfPDF_MBL.Write()
Hsave_ar40_min_MBL.Write()

Hsave_ar40PDF_TF2.Write()
Hsave_surfPDF_TF2.Write()
Hsave_ar40_min_TF2.Write()

HpmtID_pdf_ar40TF2.Write()
HpmtID_pdf_ar40MBL.Write()
