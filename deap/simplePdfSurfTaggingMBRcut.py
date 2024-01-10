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
### Hamamatsu5912 d = 204 mm; 204*2/850; 190 effective area

PI = np.pi
dist_angle = 16.5 # degree, angle between PMTs.
dist_infront = 101.6 # mm, how far infront of? = 4 inches
### NOTE: turn on MBRcut
useMBRcut = True 
MBRcut = 800

pmtpos = './Aksel_PMTpos_260.ratdb'
if not os.path.isfile(pmtpos):
    print 'ERROR: need Aksel_PMTpos_260.ratdb to work'
    print 'INFO: inside TH2DPMTfancy.py, edit the variable pmtpos to the right place'
    sys.exit(1)

deadPMT149pos = TVector3(87.805, 835.04, -132.345)
deadPMT204pos = TVector3(119.935,680.255,-495.295)
pmtIDsurround149 = [119,129,139,159,169,179]
pmtIDsurround204 = [170,179,194,195,219,220] #170 is close to gap 8

filePDF = TFile("surface_tagging_pdfs.root")

### filePDFnscb = TFile("")

print "now processing %s\n"%(sys.argv[1])
file0 = str(sys.argv[1])
fileName = os.path.basename(file0)

# load all the required PDFs
## MC position for alpha surface and ar40 PDF
Hsurf_truth_angle = filePDF.Get("surf_truth_angle_qPE")
Har40_truth_angle = filePDF.Get("ar40_truth_angle_qPE")

### !!!! Loading updated PDF
#fileNew = TFile("saveQPE_checkHistosOFF_MBR800.root")
#Har40_truth_angleReduced = fileNew.Get("Hsort_costhetaMC_scale")

#Har40_truth_angleReduced_reloaded = fileNew.Get("Hsort_costhetaMC_scale")
#Har40_truth_angleReduced_reloaded.SetName("HsortNew_costheta_MC")

### dead pmts
Hsurf_PMTs_angle  = filePDF.Get("surf_PMTs_angle_qPE")
Hsurf_gap10_angle = filePDF.Get("surf_gap10_angle_qPE")
Hsurf_gaps_angle  = filePDF.Get("surf_gaps_angle_qPE")

### build new PDFs
Hcostheta_TF2 = TH1F("Hcostheta_TF2","costheta TF2",255,-1,1)
Hcostheta_MBL = TH1F("Hcostheta_MBL","costheta MBL", 255,-1,1)

### before sorting; building truth_pdfs
HcosthetaMC  = TH1F("HcosthetaMC","costheta MC truth, pmt.qPE weighted by default", 255,-1,1) 
HcosthetaMC_qnscb = TH1F("HcosthetaMC_qnscb","pmt.qNSCB weighted costheta MC truth", 255,-1,1) 
HcosthetaMC_qprompt = TH1F("HcosthetaMC_qprompt","pmt.qPrompt weighted costheta MC truth", 255,-1,1) 

### sorted; building actual PDFs, ar40_truth_angle
Hsort_costhetaMC = TH1F("Hsort_costhetaMC","sorted costheta MC truth, default PDF", 255,0,255)
Hsort_costhetaMC_qnscb = TH1F("Hsort_costhetaMC_qnscb","sorted as PDF, pmt.qNSCB weighted costheta MC truth", 255, 0, 255)
Hsort_costhetaMC_qprompt = TH1F("Hsort_costhetaMC_qprompt","sorted as PDF, pmt.qPrompt weighted costheta MC truth", 255, 0, 255)

Hsort_costheta_MBL = TH1F("Hsort_costheta_MBL","sorted costheta MBL truth", 255,0,255)
Hsort_costheta_TF2 = TH1F("Hsort_costheta_TF2","sorted costheta TF2 truth", 255,0,255)
Hsort_costhetaMC_NLL = TH1F("Hsort_costhetaMC_NLL", "sorted costheta MC using NLL", 255, 0, 255)
Hsort_costhetaMC_NLL_scale = TH1F("Hsort_costhetaMC_NLL_scale","sorted costheta MC truth, scaled by dead PMTs", 255,0,255)

HpmtID_ar40TF2_NLL = TH1F("HpmtID_ar40TF2_NLL","old ar40 NLL",255,0,255)
HpmtID_surfTF2_NLL = TH1F("HpmtID_surfTF2_NLL","old surf NLL",255,0,255)
HpmtID_ar40TF2_NLL_scale = TH1F("HpmtID_ar40TF2_NLL_scale","scaled ar40 NLL",255,0,255)
HpmtID_surfTF2_NLL_scale = TH1F("HpmtID_surfTF2_NLL_scale","scaled surf NLL",255,0,255)

HcosthetaMC_badPMT149 = TH1F("HcosthetaMC_badPMT149","costheta MC truth, iPMT == 149", 255, -1, 1)
HcosthetaMC_badPMT204 = TH1F("HcosthetaMC_badPMT204","costheta MC truth, iPMT == 204", 255, -1, 1)
Hcostheta_deadPMT = TH1F("Hcostheta_deadPMT","costheta bad PMT", 255, -1, 1)
Hsort_costheta_deadPMT = TH1F("Hsort_costheta_deadPMT","sorted costheta bad PMT", 255, 0, 255)

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
HpmtIDnew_pdf_ar40TF2 = TH1F("HpmtIDnew_pdf_ar40TF2", "HpmtIDnew_pdf_ar40TF2", 255, 0, 255)
HpmtIDnew_pdf_ar40MBL = TH1F("HpmtIDnew_pdf_ar40MBL", "HpmtIDnew_pdf_ar40MBL", 255, 0, 255) 

### updated PDF
HsortNew_costheta_ar40TF2  = TH1F("HsortNew_costheta_ar40TF2", "new sorted nll ar40 TF2", 255, 0, 255)
HpmtIDnew_costheta_ar40TF2 = TH1F("HpmtIDnew_costheta_ar40TF2", "new pmtID nll ar40 TF2", 255, 0, 255)
HsortNew_costheta_ar40MBL  = TH1F("HsortNew_costheta_ar40MBL", "new sorted nll ar40 MBL", 255, 0, 255)
HpmtIDnew_costheta_ar40MBL = TH1F("HpmtIDnew_costheta_ar40MBL", "new pmtID nll ar40 MBL", 255, 0, 255)

Hsave_pmtsPDF = TH1F("Hsave_pmtsPDF","min(NLL_pmtsPDF)", 2000,-1000,1000)
Hsave_gapsPDF = TH1F("Hsave_gapsPDF","min(NLL_gapsPDF)", 2000,-1000,1000)

Hsave_ar40PDF_MBL = TH1F("Hsave_ar40PDF_MBL","ar40PDF",2000,-1000,1000)
Hsave_surfPDF_MBL = TH1F("Hsave_surfPDF_MBL","surfPDF",2000,-1000,1000)
Hsave_ar40_min_MBL = TH1F("Hsave_ar40_min_MBL", "ar40_min_MBL", 2000,-1000,1000)

Hsave_ar40PDF_TF2 = TH1F("Hsave_ar40PDF_TF2","ar40PDF",2000,-1000,1000)
Hsave_surfPDF_TF2 = TH1F("Hsave_surfPDF_TF2","surfPDF",2000,-1000,1000)
Hsave_ar40_min_TF2 = TH1F("Hsave_ar40min_TF2","ar40_min_TF2",2000,-1000,1000) 

#HsaveNew_ar40PDF_TF2 = TH1F("HsaveNew_ar40PDF_TF2","new ar40PDF",2000,-1000,1000)
#HsaveNew_surfPDF_TF2 = TH1F("HsaveNew_surfPDF_TF2","surfPDF",2000,-1000,1000)

HsaveNew_ar40_min_TF2 = TH1F("HsaveNew_ar40_min_TF2","scaled by DEAD pmts, ar40_min_TF2",2000,-1000,1000)
HsaveNew_surf_min_TF2 = TH1F("HsaveNew_surf_min_TF2","scaled by DEAD pmts, surf_min_TF2",2000,-1000,1000)

HdeadPMTtoAr40    = TH1F("HdeadPMTtoAr40", "deadPMT - Ar40", 2000,-1000, 1000)
HdeadPMTtoAr40new = TH1F("HdeadPMTtoAr40new","deadPMT - Ar40, new", 2000,-1000, 1000) 

Hpdf_deadPMTminusAr40    = TH1F("Hpdf_deadPMTminusAr40", "deadPMT - Ar40",255,0,255)
Hpdf_deadPMTminusAr40new = TH1F("Hpdf_deadPMTminusAr40new","deadPMT - Ar40, new",255,0,255)

### same for nscb
### Save to a new tree!!!
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
nPMTs = array('I',[0])

checkInfront149 = array('I',[0])
checkInfront204 = array('I',[0])

## sum all PMTs for ar40 and surf pdfs
sumAr40MBL_qBar = array('f',[0])
sumAr40MBL_NLL  = array('f',[0])
sumSurfMBL_qBar = array('f',[0])
sumSurfMBL_NLL = array('f',[0])

sumAr40TF2_qBar = array('f',[0])
sumAr40TF2_NLL  = array('f',[0])
sumSurfTF2_qBar = array('f',[0])
sumSurfTF2_NLL = array('f',[0])

### separate 149 and 204 for ar40 and surf pdfs
sumAr40MBL_149_qBar = array('f',[0])
sumAr40MBL_149_NLL  = array('f',[0])
sumAr40MBL_204_qBar = array('f',[0])
sumAr40MBL_204_NLL  = array('f',[0])
sumSurfMBL_149_qBar = array('f',[0])
sumSurfMBL_149_NLL = array('f',[0])
sumSurfMBL_204_qBar = array('f',[0])
sumSurfMBL_204_NLL = array('f',[0])

sumAr40TF2_149_qBar = array('f',[0])
sumAr40TF2_149_NLL  = array('f',[0])
sumAr40TF2_204_qBar = array('f',[0])
sumAr40TF2_204_NLL  = array('f',[0])
sumSurfTF2_149_qBar = array('f',[0])
sumSurfTF2_149_NLL = array('f',[0])
sumSurfTF2_204_qBar = array('f',[0])
sumSurfTF2_204_NLL = array('f',[0])

### surrounding 149 or 204 PMTs for ar40 and surf pdfs
sumSurroundAr40MBL_149_qBar = array('f',[0])
sumSurroundAr40MBL_149_NLL  = array('f',[0])
sumSurroundAr40MBL_204_qBar = array('f',[0])
sumSurroundAr40MBL_204_NLL  = array('f',[0])
sumSurroundAr40TF2_149_qBar = array('f',[0])
sumSurroundAr40TF2_149_NLL  = array('f',[0])
sumSurroundAr40TF2_204_qBar = array('f',[0])
sumSurroundAr40TF2_204_NLL  = array('f',[0])
## -- surf --
sumSurroundSurfMBL_149_qBar = array('f',[0])
sumSurroundSurfMBL_149_NLL = array('f',[0])
sumSurroundSurfMBL_204_qBar = array('f',[0])
sumSurroundSurfMBL_204_NLL = array('f',[0])
sumSurroundSurfTF2_149_qBar = array('f',[0])
sumSurroundSurfTF2_149_NLL = array('f',[0])
sumSurroundSurfTF2_204_qBar = array('f',[0])
sumSurroundSurfTF2_204_NLL = array('f',[0])

## sum two dead PMTs
sumDead149_qBar = array('f',[0])
sumDead149_NLL = array('f',[0])
sumDead204_qBar = array('f',[0])
sumDead204_NLL = array('f',[0])
## surround 149 or 204 PMTs
sumSurroundDead149_qBar = array('f',[0])
sumSurroundDead149_NLL  = array('f',[0])
sumSurroundDead204_qBar = array('f',[0])
sumSurroundDead204_NLL = array('f',[0])

##### NOTE: checking for qnscb!!!
## sum all PMTs for ar40 and surf pdfs
sumAr40MBL_qnscb = array('f',[0])
sumAr40MBL_NLL_qnscb  = array('f',[0])
sumSurfMBL_qnscb = array('f',[0])
sumSurfMBL_NLL_qnscb = array('f',[0])

sumAr40TF2_qnscb = array('f',[0])
sumAr40TF2_NLL_qnscb  = array('f',[0])
sumSurfTF2_qnscb = array('f',[0])
sumSurfTF2_NLL_qnscb = array('f',[0])

### separate 149 and 204 for ar40 and surf pdfs
sumAr40MBL_149_qnscb = array('f',[0])
sumAr40MBL_149_NLL_qnscb  = array('f',[0])
sumAr40MBL_204_qnscb = array('f',[0])
sumAr40MBL_204_NLL_qnscb  = array('f',[0])
sumSurfMBL_149_qnscb = array('f',[0])
sumSurfMBL_149_NLL_qnscb = array('f',[0])
sumSurfMBL_204_qnscb = array('f',[0])
sumSurfMBL_204_NLL_qnscb = array('f',[0])

sumAr40TF2_149_qnscb = array('f',[0])
sumAr40TF2_149_NLL_qnscb  = array('f',[0])
sumAr40TF2_204_qnscb = array('f',[0])
sumAr40TF2_204_NLL_qnscb  = array('f',[0])
sumSurfTF2_149_qnscb = array('f',[0])
sumSurfTF2_149_NLL_qnscb = array('f',[0])
sumSurfTF2_204_qnscb = array('f',[0])
sumSurfTF2_204_NLL_qnscb = array('f',[0])

### surrounding 149 or 204 PMTs for ar40 and surf pdfs
sumSurroundAr40MBL_149_qnscb = array('f',[0])
sumSurroundAr40MBL_149_NLL_qnscb  = array('f',[0])
sumSurroundAr40MBL_204_qnscb = array('f',[0])
sumSurroundAr40MBL_204_NLL_qnscb  = array('f',[0])
### -- surf --
sumSurroundSurfMBL_149_qnscb = array('f',[0])
sumSurroundSurfMBL_149_NLL_qnscb = array('f',[0])
sumSurroundSurfMBL_204_qnscb = array('f',[0])
sumSurroundSurfMBL_204_NLL_qnscb = array('f',[0])

sumSurroundAr40TF2_149_qnscb = array('f',[0])
sumSurroundAr40TF2_149_NLL_qnscb  = array('f',[0])
sumSurroundAr40TF2_204_qnscb = array('f',[0])
sumSurroundAr40TF2_204_NLL_qnscb  = array('f',[0])
### -- surf --
sumSurroundSurfTF2_149_qnscb = array('f',[0])
sumSurroundSurfTF2_149_NLL_qnscb = array('f',[0])
sumSurroundSurfTF2_204_qnscb = array('f',[0])
sumSurroundSurfTF2_204_NLL_qnscb = array('f',[0])

## sum two dead PMTs
sumDead149_qnscb = array('f',[0])
sumDead149_NLL_qnscb = array('f',[0])
sumDead204_qnscb = array('f',[0])
sumDead204_NLL_qnscb = array('f',[0])
## surround 149 or 204 PMTs
sumSurroundDead149_qnscb = array('f',[0])
sumSurroundDead149_NLL_qnscb  = array('f',[0])
sumSurroundDead204_qnscb = array('f',[0])
sumSurroundDead204_NLL_qnscb = array('f',[0])
#####################################################

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

## NOTE: fout must put before tree2
fout = TFile("testSurround_MBRcut800_"+fileName,"recreate")
## new tree2
tree2 = TTree("T2","save nll")
tree2.Branch("runID",runIDval,"runID/L")
tree2.Branch("subrunID", subrunIDval, "subrunID/L")
tree2.Branch("eventID",evtID,"eventID/L") #ULong64_
tree2.Branch("qpe", qpeVal, "qpe/F")
tree2.Branch("fprompt", fpromptVal, "fprompt/F")
tree2.Branch("mcx", mcx, "mcx/F")
tree2.Branch("mcy", mcy, "mcy/F")
tree2.Branch("mcz", mcz, "mcz/F")
tree2.Branch("mbx", mbx, "mbx/F")
tree2.Branch("mby", mby, "mby/F")
tree2.Branch("mbz", mbz, "mbz/F")
tree2.Branch("nPMTs",nPMTs, "nPMTs/I")
tree2.Branch("pmtPhi", pmtPhi, "pmtPhi[nPMTs]/F")
tree2.Branch("pmtCosTheta", pmtCosTheta, "pmtCosTheta[nPMTs]/F")
tree2.Branch("checkInfront149", checkInfront149, "checkInfront149/I")
tree2.Branch("checkInfront204", checkInfront204, "checkInfront204/I")

#### NOTE: save total and surrounding, using pmt.qpe
####################################################################################################
## ar40 and surf pdf for MBL, all PMTs
tree2.Branch("sumAr40MBL_qBar", sumAr40MBL_qBar, "sumAr40MBL_qBar/F")
tree2.Branch("sumAr40MBL_NLL", sumAr40MBL_NLL, "sumAr40MBL_NLL/F")

tree2.Branch("sumSurfMBL_qBar", sumSurfMBL_qBar,"sumSurfMBL_qBar/F")
tree2.Branch("sumSurfMBL_NLL", sumSurfMBL_204_NLL, "sumSurfMBL_NLL/F")
## ar40 pdf for MBL, 149 or 204 PMTs (for scaling)
tree2.Branch("sumAr40MBL_149_qBar", sumAr40MBL_149_qBar, "sumAr40MBL_149_qBar/F") 
tree2.Branch("sumAr40MBL_149_NLL", sumAr40MBL_149_NLL, "sumAr40MBL_149_NLL/F")
tree2.Branch("sumAr40MBL_204_qBar", sumAr40MBL_204_qBar, "sumAr40MBL_204_qBar/F")
tree2.Branch("sumAr40MBL_204_NLL", sumAr40MBL_204_NLL, "sumAr40MBL_204_NLL/F") 

tree2.Branch("sumSurfMBL_149_qBar", sumSurfMBL_149_qBar, "sumSurfMBL_149_qBar/F")
tree2.Branch("sumSurfMBL_149_NLL", sumSurfMBL_149_NLL, "sumSurfMBL_149_NLL/F")
tree2.Branch("sumSurfMBL_204_qBar", sumSurfMBL_204_qBar, "sumSurfMBL_204_qBar/F")
tree2.Branch("sumSurfMBL_204_NLL", sumSurfMBL_204_NLL, "sumSurfMBL_204_NLL/F")
###### surrounding 149 or 204 for MBL
tree2.Branch("sumSurroundAr40MBL_149_qBar", sumSurroundAr40MBL_149_qBar, "sumSurroundAr40MBL_149_qBar/F")
tree2.Branch("sumSurroundAr40MBL_149_NLL", sumSurroundAr40MBL_149_NLL, "sumSurroundAr40MBL_149_NLL/F")
tree2.Branch("sumSurroundAr40MBL_204_qBar", sumSurroundAr40MBL_204_qBar,"sumSurroundAr40MBL_204_qBar/F")
tree2.Branch("sumSurroundAr40MBL_204_NLL", sumSurroundAr40MBL_204_NLL, "sumSurroundAr40MBL_204_NLL/F")
tree2.Branch("sumSurroundSurfMBL_149_qBar", sumSurroundSurfMBL_149_qBar, "sumSurroundSurfMBL_149_qBar/F")
tree2.Branch("sumSurroundSurfMBL_149_NLL", sumSurroundSurfMBL_149_NLL, "sumSurroundSurfMBL_149_NLL/F")
tree2.Branch("sumSurroundSurfMBL_204_qBar", sumSurroundSurfMBL_204_qBar,"sumSurroundSurfMBL_204_qBar/F")
tree2.Branch("sumSurroundSurfMBL_204_NLL", sumSurroundSurfMBL_204_NLL, "sumSurroundSurfMBL_204_NLL/F")

## ar40 and surf pdf for TF2, all PMTs
tree2.Branch("sumAr40TF2_qBar", sumAr40TF2_qBar, "sumAr40TF2_qBar/F")
tree2.Branch("sumAr40TF2_NLL", sumAr40TF2_NLL, "sumAr40TF2_NLL/F")

tree2.Branch("sumSurfTF2_qBar", sumSurfTF2_qBar,"sumSurfTF2_qBar/F")
tree2.Branch("sumSurfTF2_NLL", sumSurfTF2_204_NLL, "sumSurfTF2_NLL/F")
## ar40 pdf for TF2, 149 or 204 PMTs (for scaling)
tree2.Branch("sumAr40TF2_149_qBar", sumAr40TF2_149_qBar, "sumAr40TF2_149_qBar/F") 
tree2.Branch("sumAr40TF2_149_NLL", sumAr40TF2_149_NLL, "sumAr40TF2_149_NLL/F")
tree2.Branch("sumAr40TF2_204_qBar", sumAr40TF2_204_qBar, "sumAr40TF2_204_qBar/F")
tree2.Branch("sumAr40TF2_204_NLL", sumAr40TF2_204_NLL, "sumAr40TF2_204_NLL/F") 

tree2.Branch("sumSurfTF2_149_qBar", sumSurfTF2_149_qBar, "sumSurfTF2_149_qBar/F")
tree2.Branch("sumSurfTF2_149_NLL", sumSurfTF2_149_NLL, "sumSurfTF2_149_NLL/F")
tree2.Branch("sumSurfTF2_204_qBar", sumSurfTF2_204_qBar, "sumSurfTF2_204_qBar/F")
tree2.Branch("sumSurfTF2_204_NLL", sumSurfTF2_204_NLL, "sumSurfTF2_204_NLL/F")
###### surrounding 149 or 204 for TF2
tree2.Branch("sumSurroundAr40TF2_149_qBar", sumSurroundAr40TF2_149_qBar, "sumSurroundAr40TF2_149_qBar/F")
tree2.Branch("sumSurroundAr40TF2_149_NLL", sumSurroundAr40TF2_149_NLL, "sumSurroundAr40TF2_149_NLL/F")
tree2.Branch("sumSurroundAr40TF2_204_qBar", sumSurroundAr40TF2_204_qBar,"sumSurroundAr40TF2_204_qBar/F")
tree2.Branch("sumSurroundAr40TF2_204_NLL", sumSurroundAr40TF2_204_NLL, "sumSurroundAr40TF2_204_NLL/F")
tree2.Branch("sumSurroundSurfTF2_149_qBar", sumSurroundSurfTF2_149_qBar, "sumSurroundSurfTF2_149_qBar/F")
tree2.Branch("sumSurroundSurfTF2_149_NLL", sumSurroundSurfTF2_149_NLL, "sumSurroundSurfTF2_149_NLL/F")
tree2.Branch("sumSurroundSurfTF2_204_qBar", sumSurroundSurfTF2_204_qBar,"sumSurroundSurfTF2_204_qBar/F")
tree2.Branch("sumSurroundSurfTF2_204_NLL", sumSurroundSurfTF2_204_NLL, "sumSurroundSurfTF2_204_NLL/F")

#################### for 2 dead PMTs: sum and surrounding
tree2.Branch("sumDead149_qBar",sumDead149_qBar, "sumDead149_qBar/F")
tree2.Branch("sumDead149_NLL",sumDead149_NLL, "sumDead149_NLL/F")
tree2.Branch("sumDead204_qBar", sumDead204_qBar, "sumDead204_qBar/F")
tree2.Branch("sumDead204_NLL",sumDead204_NLL, "sumDead204_NLL/F")

tree2.Branch("sumSurroundDead149_qBar",sumSurroundDead149_qBar,"sumSurroundDead149_qBar/F")
tree2.Branch("sumSurroundDead149_NLL",sumSurroundDead149_NLL,"sumSurroundDead149_NLL/F")
tree2.Branch("sumSurroundDead204_qBar",sumSurroundDead204_qBar,"sumSurroundDead204_qBar/F")
tree2.Branch("sumSurroundDead204_NLL",sumSurroundDead204_NLL,"sumSurroundDead204_NLL/F")
####################################################################################################
## !!! using pmt.qnscb
## ar40 and surf pdf for MBL, all PMTs
tree2.Branch("sumAr40MBL_qnscb", sumAr40MBL_qnscb, "sumAr40MBL_qnscb/F")
tree2.Branch("sumAr40MBL_NLL_qnscb", sumAr40MBL_NLL_qnscb, "sumAr40MBL_NLL_qnscb/F")

tree2.Branch("sumSurfMBL_qnscb", sumSurfMBL_qnscb,"sumSurfMBL_qnscb/F")
tree2.Branch("sumSurfMBL_NLL_qnscb", sumSurfMBL_204_NLL_qnscb, "sumSurfMBL_NLL_qnscb/F")
## ar40 pdf for MBL, 149 or 204 PMTs (for scaling)
tree2.Branch("sumAr40MBL_149_qnscb", sumAr40MBL_149_qnscb, "sumAr40MBL_149_qnscb/F") 
tree2.Branch("sumAr40MBL_149_NLL_qnscb", sumAr40MBL_149_NLL_qnscb, "sumAr40MBL_149_NLL_qnscb/F")
tree2.Branch("sumAr40MBL_204_qnscb", sumAr40MBL_204_qnscb, "sumAr40MBL_204_qnscb/F")
tree2.Branch("sumAr40MBL_204_NLL_qnscb", sumAr40MBL_204_NLL_qnscb, "sumAr40MBL_204_NLL_qnscb/F") 
###### surrounding 149 or 204 for MBL
tree2.Branch("sumSurroundAr40MBL_149_qnscb", sumSurroundAr40MBL_149_qnscb, "sumSurroundAr40MBL_149_qnscb/F")
tree2.Branch("sumSurroundAr40MBL_149_NLL_qnscb", sumSurroundAr40MBL_149_NLL_qnscb, "sumSurroundAr40MBL_149_NLL_qnscb/F")
tree2.Branch("sumSurroundAr40MBL_204_qnscb", sumSurroundAr40MBL_204_qnscb,"sumSurroundAr40MBL_204_qnscb/F")
tree2.Branch("sumSurroundAr40MBL_204_NLL_qnscb", sumSurroundAr40MBL_204_NLL_qnscb, "sumSurroundAr40MBL_204_NLL_qnscb/F")

tree2.Branch("sumSurfMBL_149_qnscb", sumSurfMBL_149_qnscb, "sumSurfMBL_149_qnscb/F")
tree2.Branch("sumSurfMBL_149_NLL_qnscb", sumSurfMBL_149_NLL_qnscb, "sumSurfMBL_149_NLL_qnscb/F")
tree2.Branch("sumSurfMBL_204_qnscb", sumSurfMBL_204_qnscb, "sumSurfMBL_204_qnscb/F")
tree2.Branch("sumSurfMBL_204_NLL_qnscb", sumSurfMBL_204_NLL_qnscb, "sumSurfMBL_204_NLL_qnscb/F")
###### surrounding 149 or 204 for surf TF2
tree2.Branch("sumSurroundSurfMBL_149_qnscb", sumSurroundSurfMBL_149_qnscb, "sumSurroundSurfMBL_149_qnscb/F")
tree2.Branch("sumSurroundSurfMBL_149_NLL_qnscb", sumSurroundSurfMBL_149_NLL_qnscb, "sumSurroundSurfMBL_149_NLL_qnscb/F")
tree2.Branch("sumSurroundSurfMBL_204_qnscb", sumSurroundSurfMBL_204_qnscb,"sumSurroundSurfMBL_204_qnscb/F")
tree2.Branch("sumSurroundSurfMBL_204_NLL_qnscb", sumSurroundSurfMBL_204_NLL_qnscb, "sumSurroundSurfMBL_204_NLL_qnscb/F")

## ar40 and surf pdf for TF2, all PMTs
tree2.Branch("sumAr40TF2_qnscb", sumAr40TF2_qnscb, "sumAr40TF2_qnscb/F")
tree2.Branch("sumAr40TF2_NLL_qnscb", sumAr40TF2_NLL_qnscb, "sumAr40TF2_NLL_qnscb/F")

tree2.Branch("sumSurfTF2_qnscb", sumSurfTF2_qnscb,"sumSurfTF2_qnscb/F")
tree2.Branch("sumSurfTF2_NLL_qnscb", sumSurfTF2_204_NLL_qnscb, "sumSurfTF2_NLL_qnscb/F")
## ar40 pdf for TF2, 149 or 204 PMTs (for scaling)
tree2.Branch("sumAr40TF2_149_qnscb", sumAr40TF2_149_qnscb, "sumAr40TF2_149_qnscb/F") 
tree2.Branch("sumAr40TF2_149_NLL_qnscb", sumAr40TF2_149_NLL_qnscb, "sumAr40TF2_149_NLL_qnscb/F")
tree2.Branch("sumAr40TF2_204_qnscb", sumAr40TF2_204_qnscb, "sumAr40TF2_204_qnscb/F")
tree2.Branch("sumAr40TF2_204_NLL_qnscb", sumAr40TF2_204_NLL_qnscb, "sumAr40TF2_204_NLL_qnscb/F") 
###### surrounding 149 or 204 for ar40 TF2
tree2.Branch("sumSurroundAr40TF2_149_qnscb", sumSurroundAr40TF2_149_qnscb, "sumSurroundAr40TF2_149_qnscb/F")
tree2.Branch("sumSurroundAr40TF2_149_NLL_qnscb", sumSurroundAr40TF2_149_NLL_qnscb, "sumSurroundAr40TF2_149_NLL_qnscb/F")
tree2.Branch("sumSurroundAr40TF2_204_qnscb", sumSurroundAr40TF2_204_qnscb,"sumSurroundAr40TF2_204_qnscb/F")
tree2.Branch("sumSurroundAr40TF2_204_NLL_qnscb", sumSurroundAr40TF2_204_NLL_qnscb, "sumSurroundAr40TF2_204_NLL_qnscb/F")

tree2.Branch("sumSurfTF2_149_qnscb", sumSurfTF2_149_qnscb, "sumSurfTF2_149_qnscb/F")
tree2.Branch("sumSurfTF2_149_NLL_qnscb", sumSurfTF2_149_NLL_qnscb, "sumSurfTF2_149_NLL_qnscb/F")
tree2.Branch("sumSurfTF2_204_qnscb", sumSurfTF2_204_qnscb, "sumSurfTF2_204_qnscb/F")
tree2.Branch("sumSurfTF2_204_NLL_qnscb", sumSurfTF2_204_NLL_qnscb, "sumSurfTF2_204_NLL_qnscb/F")
###### surrounding 149 or 204 for surf TF2
tree2.Branch("sumSurroundSurfTF2_149_qnscb", sumSurroundSurfTF2_149_qnscb, "sumSurroundSurfTF2_149_qnscb/F")
tree2.Branch("sumSurroundSurfTF2_149_NLL_qnscb", sumSurroundSurfTF2_149_NLL_qnscb, "sumSurroundSurfTF2_149_NLL_qnscb/F")
tree2.Branch("sumSurroundSurfTF2_204_qnscb", sumSurroundSurfTF2_204_qnscb,"sumSurroundSurfTF2_204_qnscb/F")
tree2.Branch("sumSurroundSurfTF2_204_NLL_qnscb", sumSurroundSurfTF2_204_NLL_qnscb, "sumSurroundSurfTF2_204_NLL_qnscb/F")

#################### for 2 dead PMTs: sum and surrounding
tree2.Branch("sumDead149_qnscb",sumDead149_qnscb, "sumDead149_qnscb/F")
tree2.Branch("sumDead149_NLL_qnscb",sumDead149_NLL_qnscb, "sumDead149_NLL_qnscb/F")
tree2.Branch("sumDead204_qnscb", sumDead204_qnscb, "sumDead204_qnscb/F")
tree2.Branch("sumDead204_NLL_qnscb",sumDead204_NLL_qnscb, "sumDead204_NLL_qnscb/F")

tree2.Branch("sumSurroundDead149_qnscb",sumSurroundDead149_qnscb,"sumSurroundDead149_qnscb/F")
tree2.Branch("sumSurroundDead149_NLL_qnscb",sumSurroundDead149_NLL_qnscb,"sumSurroundDead149_NLL_qnscb/F")
tree2.Branch("sumSurroundDead204_qnscb",sumSurroundDead204_qnscb,"sumSurroundDead204_qnscb/F")
tree2.Branch("sumSurroundDead204_NLL_qnscb",sumSurroundDead204_NLL_qnscb,"sumSurroundDead204_NLL_qnscb/F")
##################################################################################################

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

## two dead PMTs
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

gaps_vectors=[TVector3() for i in range(11)]
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
    gaps_vectors=[TVector3() for i in range(11)] ##Jie: move this to global
    for i in range(10):
        gaps_vectors[i].SetXYZ(gaps_arrays[i][0], gaps_arrays[i][1], gaps_arrays[i][2])
        gaps_vectors[i].SetMag(851.0)
        if i==4:
            gaps_vectors[4].SetPhi(PI)
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
    del gaps_angles, gaps_angles_dict, gaps_angles_key#, gaps_vectors # Jie: check gaps_vectors
    if len(dead_pmts):
        del pmts_angles, pmts_angles_dict, pmts_angles_key, dead_pmts_vectors

    return pmts_angles_ids, gaps_angles_ids, gaps_vectors # Jie: add gaps_vectors just for simulations


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
pmts_angles_ids, gaps_angles_ids, gaps_vectors = pmts_gaps_ids()
#-------------------------------------------------------------------------------------------------------------------------------------------------------
# to hold MBL and TF2 fitted positions
MBLEvent = TVector3()
TF2Event = TVector3()
MCEvent = TVector3()
########################################################################################################################################################
tree1 = File.Get("T1");#TTree
nentries = tree1.GetEntries();

#ntupleSurfTag = TNtuple("ntupleSurfTag","surface tag variables","eventID:NLL_surfPDF_MBL:NLL_ar40PDF_MBL:NLL_surfPDF_TF2:NLL_ar40PDF_TF2:min(NLL_pmtsPDF):min(NLL_gapsPDF):gof_surfPDF_MBL:gof_ar40PDF_MBL:gof_surfPDF_TF2:gof_ar40PDF_TF2")

#if nentries>5000:
#    nentries = 15000
#nentries = 15000
ListProb = []
for event in range(nentries):
    if (event%20000 == 0): print "processed", event
    tree1.GetEntry(event)
    eventID = tree1.eventID 
    #if eventID != 37: continue
    #print "---------- eventID=", eventID
    xmc = tree1.mcx
    ymc = tree1.mcy
    zmc = tree1.mcz

    evtx = tree1.mbx
    evty = tree1.mby
    evtz = tree1.mbz
    evtx_tf2 = tree1.tf2x
    evty_tf2 = tree1.tf2y
    evtz_tf2 = tree1.tf2z

    ### !! turn on MBR cut here
    if useMBRcut and sqrt(evtx*evtx+evty*evty+evtz*evtz)>MBRcut: continue

    MBLEvent.SetXYZ(evtx, evty, evtz)
    TF2Event.SetXYZ(evtx_tf2, evty_tf2, evtz_tf2)

    MCEvent.SetXYZ(xmc, ymc, zmc)

    pmtPhi          = tree1.pmtPhi
    pmtCosTheta     = tree1.pmtCosTheta
    pmtCharge       = tree1.pmtCharge
    pmtChargeNSCB   = tree1.pmtChargeNSCB
    pmtChargePrompt = tree1.pmtChargePrompt
    listPmtID       = tree1.listPmtID

    ###NOTE: for an MC event, check whether the event is infront of dead pmts
    is_infront_dead149 = ( acos( deadPMT149pos.Unit()*MCEvent.Unit() )< dist_angle*PI/180 and (deadPMT149pos - MCEvent).Mag()<dist_infront )
    is_infront_dead204 = ( acos( deadPMT204pos.Unit()*MCEvent.Unit() )< dist_angle*PI/180 and (deadPMT204pos - MCEvent).Mag()<dist_infront )
   
    checkInfront149[0] = int(is_infront_dead149)
    checkInfront204[0] = int(is_infront_dead204)
    
    #print "in front of what?", is_infront_dead149, "=149", is_infront_dead204, "=204"
    ### MC truth position is close to gap position
    igap = 0
    ### check events close to gaps
    #for gapPos in gaps_vectors:
    #   radius = gapPos.Mag() 
    #   gapPos.SetMag(radius - 4*25.4)
    #   #print igap, round(gapPos[0],3), round(gapPos[1],3), round(gapPos[2],3)
    #   igap += 1
    #   distance = (MCEvent - gapPos).Mag()
    #   if distance<=4*25.2 and distance>1: # 1 mm<distance<101.6 mm (4 inches)
    #      print "gap events!!"

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

    ## checking updated PDF
    NLL00_surfPDF_MBL, NLL00_ar40PDF_MBL = 0, 0 #use MBL position and the surface and ar40 PDF
    NLL00_surfPDF_TF2, NLL00_ar40PDF_TF2 = 0, 0 #use TF2 position and the surface and ar40 PDF

    gof_surfPDF_MBL , gof_ar40PDF_MBL  = 0, 0 #use MBL position and the surface and ar40 PDF
    gof_surfPDF_TF2 , gof_ar40PDF_TF2  = 0, 0 #use TF2 position and the surface and ar40 PDF

    NLL_pmtsPDF = [0 for i in range(len(pmts_angles_ids))] #use the dead PMT positions and the PMT PDFs
    NLL_gapsPDF = [0 for i in range(len(gaps_angles_ids))] #use the pentagonal gaps positions and the Gap PDFs

    gof_pmtsPDF_PMTs = [0 for i in range(len(pmts_angles_ids))] #use the dead PMT positions and the PMT PDFs
    gof_gapsPDF_gaps = [0 for i in range(len(gaps_angles_ids))] #use the pentagonal gaps positions and the Gap PDFs

    pmt_ids = []
    ### Warning !!! Here we consider the total charge of event, event qPE or nSCBayes
    qPE = tree1.qpe  # event qPE
    nSCB = tree1.nSCBayes # event nscb
    nPMT = tree1.nPMTs1
    
    ### NOTE: 90<nSCB<200 cut?
    #if nSCB>200: continue
    list_255charges = []

    probSurround1 = 0; probSurround2 = 0; ## calculate surrounding PMTs NLL
    sumProbDead = 0; ##NOTE: two dead PMTs
    listNLLmc = []
    listPDFid = []
    ### Looping 255 PMTs
    ### NOTE: To modify NLL by dead PMTs. ONLY for TF2 !!! 
    saveNLL_ar40_TF2 = {}
    saveNLL_surf_TF2 = {}
    pmt_qPE = 0
    pmt_qncb = 0
    pmt_qprompt = 0
    #### using pmt.qpe
    sumAr40MBL_qBar[0] = 0;sumAr40MBL_NLL[0] = 0
    sumSurfMBL_qBar[0] = 0;sumSurfMBL_NLL[0] = 0
    sumAr40MBL_149_qBar[0] = 0;sumAr40MBL_149_NLL[0] = 0
    sumAr40MBL_204_qBar[0] = 0;sumAr40MBL_204_NLL[0] = 0
    sumSurfMBL_149_qBar[0] = 0;sumSurfMBL_149_NLL[0] = 0
    sumSurfMBL_204_qBar[0] = 0;sumSurfMBL_204_NLL[0] = 0
    sumSurroundAr40MBL_149_qBar[0] = 0
    sumSurroundAr40MBL_149_NLL[0] = 0
    sumSurroundAr40MBL_204_qBar[0] = 0
    sumSurroundAr40MBL_204_NLL[0] = 0
    sumSurroundSurfMBL_149_qBar[0] = 0
    sumSurroundSurfMBL_149_NLL[0] = 0
    sumSurroundSurfMBL_204_qBar[0] = 0
    sumSurroundSurfMBL_204_NLL[0] = 0

    sumAr40TF2_qBar[0] = 0;sumAr40TF2_NLL[0] = 0
    sumSurfTF2_qBar[0] = 0;sumSurfTF2_NLL[0] = 0
    sumAr40TF2_149_qBar[0] = 0;sumAr40TF2_149_NLL[0] = 0
    sumAr40TF2_204_qBar[0] = 0;sumAr40TF2_204_NLL[0] = 0
    sumSurfTF2_149_qBar[0] = 0;sumSurfTF2_149_NLL[0] = 0
    sumSurfTF2_204_qBar[0] = 0;sumSurfTF2_204_NLL[0] = 0
    sumSurroundAr40TF2_149_qBar[0] = 0
    sumSurroundAr40TF2_149_NLL[0] = 0
    sumSurroundAr40TF2_204_qBar[0] = 0
    sumSurroundAr40TF2_204_NLL[0] = 0
    sumSurroundSurfTF2_149_qBar[0] = 0
    sumSurroundSurfTF2_149_NLL[0] = 0
    sumSurroundSurfTF2_204_qBar[0] = 0
    sumSurroundSurfTF2_204_NLL[0] = 0

    ### dead pmt pdfs
    sumDead149_qBar[0] = 0; sumDead149_NLL[0] = 0
    sumDead204_qBar[0] = 0; sumDead204_NLL[0] = 0
    sumSurroundDead149_qBar[0] = 0; sumSurroundDead149_NLL[0] = 0 
    sumSurroundDead204_qBar[0] = 0; sumSurroundDead204_NLL[0] = 0

    ##############################################################
    #### using pmt.qnscb
    sumAr40MBL_qnscb[0] = 0;sumAr40MBL_NLL_qnscb[0] = 0
    sumSurfMBL_qnscb[0] = 0;sumSurfMBL_NLL_qnscb[0] = 0
    sumAr40MBL_149_qnscb[0] = 0;sumAr40MBL_149_NLL_qnscb[0] = 0
    sumAr40MBL_204_qnscb[0] = 0;sumAr40MBL_204_NLL_qnscb[0] = 0
    sumSurfMBL_149_qnscb[0] = 0;sumSurfMBL_149_NLL_qnscb[0] = 0
    sumSurfMBL_204_qnscb[0] = 0;sumSurfMBL_204_NLL_qnscb[0] = 0
    sumSurroundAr40MBL_149_qnscb[0] = 0
    sumSurroundAr40MBL_149_NLL_qnscb[0] = 0
    sumSurroundAr40MBL_204_qnscb[0] = 0
    sumSurroundAr40MBL_204_NLL_qnscb[0] = 0
    sumSurroundSurfMBL_149_qnscb[0] = 0
    sumSurroundSurfMBL_149_NLL_qnscb[0] = 0
    sumSurroundSurfMBL_204_qnscb[0] = 0
    sumSurroundSurfMBL_204_NLL_qnscb[0] = 0

    sumAr40TF2_qnscb[0] = 0;sumAr40TF2_NLL_qnscb[0] = 0
    sumSurfTF2_qnscb[0] = 0;sumSurfTF2_NLL_qnscb[0] = 0
    sumAr40TF2_149_qnscb[0] = 0;sumAr40TF2_149_NLL_qnscb[0] = 0
    sumAr40TF2_204_qnscb[0] = 0;sumAr40TF2_204_NLL_qnscb[0] = 0
    sumSurfTF2_149_qnscb[0] = 0;sumSurfTF2_149_NLL_qnscb[0] = 0
    sumSurfTF2_204_qnscb[0] = 0;sumSurfTF2_204_NLL_qnscb[0] = 0
    sumSurroundAr40TF2_149_qnscb[0] = 0
    sumSurroundAr40TF2_149_NLL_qnscb[0] = 0
    sumSurroundAr40TF2_204_qnscb[0] = 0
    sumSurroundAr40TF2_204_NLL_qnscb[0] = 0
    sumSurroundSurfTF2_149_qnscb[0] = 0
    sumSurroundSurfTF2_149_NLL_qnscb[0] = 0
    sumSurroundSurfTF2_204_qnscb[0] = 0
    sumSurroundSurfTF2_204_NLL_qnscb[0] = 0
 
    ### dead pmt pdfs
    sumDead149_qnscb[0] = 0; sumDead149_NLL_qnscb[0] = 0
    sumDead204_qnscb[0] = 0; sumDead204_NLL_qnscb[0] = 0
    sumSurroundDead149_qnscb[0] = 0; sumSurroundDead149_NLL_qnscb[0] = 0 
    sumSurroundDead204_qnscb[0] = 0; sumSurroundDead204_NLL_qnscb[0] = 0

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
            pmt_qnscb = pmtChargeNSCB[iPMT]
            pmt_qprompt = pmtChargePrompt[iPMT]
            #print "what ??",pmt_id, pmt_qPE, pmt_qpropmt
            list_255charges.append(pmt_qPE)
            pmt_ids.append(pmt_id)
        else:
            # PMTs that have not observed any charge.
            pmt_id = pmt_ids[iPMT]
            pmt_qPE = 0
            pmt_qnscb = 0
            pmt_qprompt = 0
            list_255charges.append(pmt_qPE)

        ## NOTE: only for MC!!!!
        sortedAngle = MC_angles_dict[MC_angles_ids.index(pmt_id)][1] ### costheta sorted by pmt.qPE, using pmt_id
        HcosthetaMC.Fill(sortedAngle, pmt_qPE)
        HcosthetaMC_qnscb.Fill(sortedAngle, pmt_qnscb)
        HcosthetaMC_qprompt.Fill(sortedAngle, pmt_qprompt)
        HpmtID_pmtCharges.Fill(pmt_id, pmt_qPE)
         
        if MC_angles_ids.index(pmt_id) == 0: ## assume the first is dead !!
            Hcostheta_deadPMT.Fill(sortedAngle, 0)
        else: Hcostheta_deadPMT.Fill(sortedAngle, pmt_qPE)    

        #####################################################################
        ## NOTE: calculate NLL for MC to modify PDF
        pdf_index = MC_angles_ids.index(pmt_id)+1
        xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
        nll_ar40_MC = get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))

        #if pmt_id == 149:
        #    prob1 = nll_ar40_MC
        #if pmt_id == 204:
        #    prob2 = nll_ar40_MC
        #sumProbDead = prob1 + prob2
        #if sumProbDead != 0: print "sum prob", sumProb
        #Hsort_costhetaMC_scale.Fill(pdf_index-1, nll_ar40_MC)
        listPDFid.append( pdf_index-1 )
        listNLLmc.append( nll_ar40_MC )
        #if pmt_id != 149 or pmt_id != 204:
        #    Hsort_costhetaMC_scale.Fill( pdf_index-1, pmt_qPE/(1-sumProbDead) )
        #Hsort_costhetaMC_NLL.Fill( pdf_index-1, nll_ar40_MC)
        #if pmt_id == 149:
        #    print "bad!!!!", pmt_id, pmt_qPE
        #    HcosthetaMC_badPMT149.Fill(sortedAngle, pmt_qPE)
        #if pmt_id == 204:
        #    HcosthetaMC_badPMT204.Fill(sortedAngle, pmt_qPE)

        #Hsort_costhetaMC.Fill()
        #+1 is needed for GetBinContent because index 0 is underflow.
        #------------------------------------------------------------------------
        # calculate NLL and gof based on MBL position
        if MBL_fit_found:
            pdf_index = MBL_angles_ids.index(pmt_id)+1

            xbar_surfPDF = Hsurf_truth_angle.GetBinContent(pdf_index)*qPE
            xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
            nll_surf_MBL = get_poisson_nll(xbar_surfPDF, int(round(pmt_qPE)))
            nll_ar40_MBL = get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))
             
            ##NOTE: plot surrounding PMTs!!!
            sumAr40MBL_qBar[0] += xbar_ar40PDF
            sumAr40MBL_NLL[0] += nll_ar40_MBL
            sumSurfMBL_qBar[0] += xbar_surfPDF
            sumSurfMBL_NLL[0] += nll_surf_MBL

            if pmt_id == 149:
                sumAr40MBL_149_qBar[0] = xbar_ar40PDF
                sumAr40MBL_149_NLL[0] = nll_ar40_MBL
                sumSurfMBL_149_qBar[0] = xbar_surfPDF
                sumSurfMBL_149_NLL[0] = nll_surf_MBL

            if is_infront_dead149:
                if pmt_id in pmtIDsurround149:
                    sumSurroundAr40MBL_149_qBar[0] += xbar_ar40PDF
                    sumSurroundAr40MBL_149_NLL[0] += nll_ar40_MBL 
                    sumSurroundSurfMBL_149_qBar[0] += xbar_surfPDF
                    sumSurroundSurfMBL_149_NLL[0] += nll_surf_MBL

            if pmt_id == 204:
                sumAr40MBL_204_qBar[0] = xbar_ar40PDF
                sumAr40MBL_204_NLL[0] = nll_ar40_MBL
                sumSurfMBL_204_qBar[0] = xbar_surfPDF
                sumSurfMBL_204_NLL[0] = nll_surf_MBL

            if is_infront_dead204:
                if pmt_id in pmtIDsurround204:
                    sumSurroundAr40MBL_204_qBar[0] += xbar_ar40PDF
                    sumSurroundAr40MBL_204_NLL[0] += nll_ar40_MBL
                    sumSurroundSurfMBL_204_qBar[0] += xbar_surfPDF
                    sumSurroundSurfMBL_204_NLL[0] += nll_surf_MBL

            ## --------------- using pmt.qnscb
            xbar_surfPDFqnscb = Hsurf_truth_angle.GetBinContent(pdf_index)*nSCB ## change pdf??
            xbar_ar40PDFqnscb = Har40_truth_angle.GetBinContent(pdf_index)*nSCB
            nll_surf_MBLqnscb = get_poisson_nll(xbar_surfPDFqnscb, int(round(pmt_qnscb)))
            nll_ar40_MBLqnscb = get_poisson_nll(xbar_ar40PDFqnscb, int(round(pmt_qnscb)))

            sumAr40MBL_qnscb[0] += xbar_ar40PDFqnscb
            sumAr40MBL_NLL_qnscb[0] += nll_ar40_MBLqnscb
            sumSurfMBL_qnscb[0] += xbar_surfPDFqnscb
            sumSurfMBL_NLL_qnscb[0] += nll_surf_MBLqnscb

            if pmt_id == 149:
                sumAr40MBL_149_qnscb[0] = xbar_ar40PDFqnscb
                sumAr40MBL_149_NLL_qnscb[0] = nll_ar40_MBLqnscb
                sumSurfMBL_149_qnscb[0] = xbar_surfPDFqnscb
                sumSurfMBL_149_NLL_qnscb[0] = nll_surf_MBLqnscb

            if is_infront_dead149:
                if pmt_id in pmtIDsurround149:
                    sumSurroundAr40MBL_149_qnscb[0] += xbar_ar40PDFqnscb
                    sumSurroundAr40MBL_149_NLL_qnscb[0] += nll_ar40_MBLqnscb
                    sumSurroundSurfMBL_149_qnscb[0] += xbar_surfPDFqnscb
                    sumSurroundSurfMBL_149_NLL_qnscb[0] += nll_surf_MBLqnscb

            if pmt_id == 204:
                sumAr40MBL_204_qnscb[0] = xbar_ar40PDFqnscb
                sumAr40MBL_204_NLL_qnscb[0] = nll_ar40_MBLqnscb
                sumSurfMBL_204_qnscb[0] = xbar_surfPDFqnscb
                sumSurfMBL_204_NLL_qnscb[0] = nll_surf_MBLqnscb

            if is_infront_dead204:
                if pmt_id in pmtIDsurround204:
                    sumSurroundAr40MBL_204_qnscb[0] += xbar_ar40PDFqnscb
                    sumSurroundAr40MBL_204_NLL_qnscb[0] += nll_ar40_MBLqnscb
                    sumSurroundSurfMBL_204_qnscb[0] += xbar_surfPDFqnscb
                    sumSurroundSurfMBL_204_NLL_qnscb[0] += nll_surf_MBLqnscb

            #if pmt_id != 149 and pmt_id != 204:
            #   nll_ar40_MBL = nll_ar40_MBL/(1 - sumProbDead) 
            #   nll_surf_MBL = nll_surf_MBL/(1 - sumProbDead)

            NLL_ar40PDF_MBL += nll_ar40_MBL
            NLL_surfPDF_MBL += nll_surf_MBL

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
        # calculate NLL and gof based on just TF2 position
        saveNLL = 0; saveNLLnew = 0;
        if TF2_fit_found:
            pdf_index = TF2_angles_ids.index(pmt_id)+1
            xbar_surfPDF = Hsurf_truth_angle.GetBinContent(pdf_index)*qPE
            xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
            
            nll_ar40_TF2 = get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))
            nll_surf_TF2 = get_poisson_nll(xbar_surfPDF, int(round(pmt_qPE)))

            ##NOTE: plot surrounding PMTs!!!
            sumAr40TF2_qBar[0] += xbar_ar40PDF
            sumAr40TF2_NLL[0] += nll_ar40_TF2
            sumSurfTF2_qBar[0] += xbar_surfPDF
            sumSurfTF2_NLL[0] += nll_surf_TF2
            #print "each pmtID", pmt_id, pmt_qPE, xbar_ar40PDF, nll_ar40_TF2
            if pmt_id == 149:
                sumAr40TF2_149_qBar[0] += xbar_ar40PDF
                sumAr40TF2_149_NLL[0] += nll_ar40_TF2
                sumSurfTF2_149_qBar[0] += xbar_surfPDF
                sumSurfTF2_149_NLL[0] += nll_surf_TF2

            if is_infront_dead149:
                if pmt_id in pmtIDsurround149:
                    sumSurroundAr40TF2_149_qBar[0] += xbar_ar40PDF
                    sumSurroundAr40TF2_149_NLL[0] += nll_ar40_TF2
                    sumSurroundSurfTF2_149_qBar[0] += xbar_surfPDF
                    sumSurroundSurfTF2_149_NLL[0] += nll_surf_TF2

            if pmt_id == 204:
                sumAr40TF2_204_qBar[0] += xbar_ar40PDF
                sumAr40TF2_204_NLL[0] += nll_ar40_TF2
                sumSurfTF2_204_qBar[0] += xbar_surfPDF
                sumSurfTF2_204_NLL[0] += nll_surf_TF2

            if is_infront_dead204:
                if pmt_id in pmtIDsurround204:
                    sumSurroundAr40TF2_204_qBar[0] += xbar_ar40PDF
                    sumSurroundAr40TF2_204_NLL[0] += nll_ar40_TF2
                    sumSurroundSurfTF2_204_qBar[0] += xbar_surfPDF
                    sumSurroundSurfTF2_204_NLL[0] += nll_surf_TF2

            ## --------------- using pmt.qnscb
            xbar_surfPDFqnscb = Hsurf_truth_angle.GetBinContent(pdf_index)*nSCB ## change pdf??
            xbar_ar40PDFqnscb = Har40_truth_angle.GetBinContent(pdf_index)*nSCB

            nll_surf_TF2qnscb = get_poisson_nll(xbar_surfPDFqnscb, int(round(pmt_qnscb)))
            nll_ar40_TF2qnscb = get_poisson_nll(xbar_ar40PDFqnscb, int(round(pmt_qnscb)))

            sumAr40TF2_qnscb[0] += xbar_ar40PDFqnscb
            sumAr40TF2_NLL_qnscb[0] += nll_ar40_TF2qnscb
            sumSurfTF2_qnscb[0] += xbar_surfPDFqnscb
            sumSurfTF2_NLL_qnscb[0] += nll_surf_TF2qnscb

            if pmt_id == 149: ## single dead NLL
                sumAr40TF2_149_qnscb[0] = xbar_ar40PDFqnscb
                sumAr40TF2_149_NLL_qnscb[0] = nll_ar40_TF2qnscb
                sumSurfTF2_149_qnscb[0] = xbar_surfPDFqnscb
                sumSurfTF2_149_NLL_qnscb[0] = nll_surf_TF2qnscb

            if is_infront_dead149:
                if pmt_id in pmtIDsurround149:
                    sumSurroundAr40TF2_149_qnscb[0] += xbar_ar40PDFqnscb
                    sumSurroundAr40TF2_149_NLL_qnscb[0] += nll_ar40_TF2qnscb
                    sumSurroundSurfTF2_149_qnscb[0] += xbar_surfPDFqnscb
                    sumSurroundSurfTF2_149_NLL_qnscb[0] += nll_surf_TF2qnscb

            if pmt_id == 204: ## single dead NLL
                sumAr40TF2_204_qnscb[0] = xbar_ar40PDFqnscb
                sumAr40TF2_204_NLL_qnscb[0] = nll_ar40_TF2qnscb
                sumSurfTF2_204_qnscb[0] = xbar_surfPDFqnscb
                sumSurfTF2_204_NLL_qnscb[0] = nll_surf_TF2qnscb

            if is_infront_dead204:
                if pmt_id in pmtIDsurround204:
                    sumSurroundAr40TF2_204_qnscb[0] += xbar_ar40PDFqnscb
                    sumSurroundAr40TF2_204_NLL_qnscb[0] += nll_ar40_TF2qnscb
                    sumSurroundSurfTF2_204_qnscb[0] += xbar_surfPDFqnscb
                    sumSurroundSurfTF2_204_NLL_qnscb[0] += nll_surf_TF2qnscb

            #print nll_surf_TF2, nll_ar40_TF2, 1-sumProbDead
            #if pmt_id != 149 and pmt_id != 204:
            #   nll_ar40_MBL = nll_ar40_MBL/(1 - sumProbDead)
            #   nll_surf_MBL = nll_surf_MBL/(1 - sumProbDead)

            saveNLL_ar40_TF2[pmt_id] = nll_ar40_TF2
            saveNLL_surf_TF2[pmt_id] = nll_surf_TF2

            NLL_ar40PDF_TF2 += nll_ar40_TF2
            NLL_surfPDF_TF2 += nll_surf_TF2

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
        saveListDeadPMTs = []

        for i in range(len(pmts_angles_ids)): ##2 dead PMTs
            pdf_index = pmts_angles_ids[i].index(pmt_id)+1
            overwrite = True if pmt_id in dead_pmt_id else False # avoid X_deapPMT*X_deadPMT
            xbar_deadPDF = 0; xbar_deadPDFqnscb = 0;
            if not overwrite: ## pmt_id not in dead pmts
                xbar_deadPDF = Hsurf_PMTs_angle.GetBinContent(pdf_index)*qPE
                xbar_deadPDFqnscb = Hsurf_PMTs_angle.GetBinContent(pdf_index)*nSCB
            else:
                # If a (second) dead PMT, expection would be also zero charge.
                xbar_deadPDF = Hsurf_PMTs_angle.GetBinContent(1)*qPE
                xbar_deadPDFqnscb = Hsurf_PMTs_angle.GetBinContent(1)*nSCB
            
            ## here pmt_qPE is looped!!
            nllVal = get_poisson_nll(xbar_deadPDF, int(round(pmt_qPE)))
            NLL_pmtsPDF[i] += nllVal 
            gof_pmtsPDF_PMTs[i] += gof_calculator(xbar_deadPDF, pmt_qPE)
            #if eventID != checkEventID: continue
            ## sorted PMT index
            ListHDump_NLL_deadPMT[i].Fill(pdf_index-1, nllVal)
            ListHDump_Q_deadPMT[i].Fill(pdf_index-1, pmt_qPE)
            ListHDump_pdf_deadPMT[i].Fill(pdf_index-1, xbar_deadPDF)
            ## PMT ids
            ListHDumpPMTid_NLL_deadPMT[i].Fill(pmt_id, nllVal)
            saveListDeadPMTs.append(nllVal)
            ListHDumpPMTid_Q_deadPMT[i].Fill(pmt_id, pmt_qPE)
            ListHDumpPMTid_pdf_deadPMT[i].Fill(pmt_id, xbar_deadPDF)

            if i == 0: # dead 149
                sumDead149_qBar[0] += xbar_deadPDF
                sumDead149_NLL[0] += nllVal
            
            if i == 1: # dead 204 
                sumDead204_qBar[0] += xbar_deadPDF
                sumDead204_NLL[0] += nllVal

            ##NOTE: plot surrounding PMTs!!!
            if i == 0 and is_infront_dead149:
                if pmt_id in pmtIDsurround149: 
                    sumSurroundDead149_qBar[0] += xbar_deadPDF
                    sumSurroundDead149_NLL[0] += nllVal

            if i==1 and is_infront_dead204:
                if pmt_id in pmtIDsurround204:
                    sumSurroundDead204_qBar[0] += xbar_deadPDF
                    sumSurroundDead204_NLL[0] += nllVal

            nllVal_qnscb = get_poisson_nll(xbar_deadPDFqnscb, int(round(pmt_qnscb)))
            if i == 0: # dead 149
                sumDead149_qnscb[0] += xbar_deadPDFqnscb
                sumDead149_NLL_qnscb[0] += nllVal_qnscb
            if i == 1:
                sumDead204_qnscb[0] += xbar_deadPDFqnscb
                sumDead204_NLL_qnscb[0] += nllVal_qnscb

            if i == 0 and is_infront_dead149:
                if pmt_id in pmtIDsurround149:
                    sumSurroundDead149_qnscb[0] += xbar_deadPDFqnscb
                    sumSurroundDead149_NLL_qnscb[0] += nllVal_qnscb

            if i == 1 and is_infront_dead204:
                if pmt_id in pmtIDsurround204:
                    sumSurroundDead204_qnscb[0] += xbar_deadPDFqnscb
                    sumSurroundDead204_NLL_qnscb[0] += nllVal_qnscb

        ### check the difference!!!!
        if TF2_fit_found:
            pdf_index = TF2_angles_ids.index(pmt_id)+1   
            saveNllMinDeadPMT = min(saveListDeadPMTs)
            Hpdf_deadPMTminusAr40.Fill(pmt_id, saveNllMinDeadPMT - saveNLL)
            Hpdf_deadPMTminusAr40new.Fill(pmt_id, saveNllMinDeadPMT - saveNLLnew)
    
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
    NLL_tag_MBL_list = [NLL_ar40PDF_MBL, NLL_surfPDF_MBL, min(NLL_pmtsPDF), min(NLL_gapsPDF)]
    NLL_tag_TF2_list = [NLL_ar40PDF_TF2, NLL_surfPDF_TF2, min(NLL_pmtsPDF), min(NLL_gapsPDF)]

    is_surface_like_MBL = 0 if min(NLL_tag_MBL_list)==NLL_tag_MBL_list[0] else 1
    is_surface_like_TF2 = 0 if min(NLL_tag_TF2_list)==NLL_tag_TF2_list[0] else 1

    # Mainly for events tag as ar40 --> How far/close was the next minimum NLL?!

    NLL_ar40_min_MBL = NLL_ar40PDF_MBL - min(NLL_tag_MBL_list[1:])
    NLL_ar40_min_TF2 = NLL_ar40PDF_TF2 - min(NLL_tag_TF2_list[1:])
    
    #######################################################################
    ### exist PMT loop and scale PDF
    #----------------------------------------------------------------------
    sumNLLmc = sum(listNLLmc)
    #print "check scale", listPDFid
    scale1 = 1 - sumProbDead/(sumNLLmc - sumProbDead)
    scale2 = 1 - sumProbDead/sumNLLmc
    #print prob1, prob2, sumNLLmc, scale1, scale2
    #for i in range( len(listNLLmc) ):
    #    Hsort_costhetaMC_NLL.Fill(listPDFid[i], listNLLmc[i])
    #    Hsort_costhetaMC_NLL_scale.Fill(listPDFid[i], listNLLmc[i]/scale2)
    #######################################################################
    ##NOTE: recalculate NLL!!!!
    new_NLL_ar40 = 0
    new_NLL_surf = 0
    for pmt_id in saveNLL_ar40_TF2: ## same to saveNLL_surf_TF2, 0 to 255
         nll_ar40 = saveNLL_ar40_TF2[pmt_id]
         nll_surf = saveNLL_surf_TF2[pmt_id]
         #HpmtID_ar40TF2_NLL.Fill(pmt_id, nll_ar40)
         #HpmtID_surfTF2_NLL.Fill(pmt_id, nll_surf)

         if pmt_id != 149 and pmt_id != 204:
             nll_ar40 = nll_ar40/scale2
             nll_surf = nll_surf/scale2
         new_NLL_ar40 += nll_ar40
         new_NLL_surf += nll_surf
         #HpmtID_ar40TF2_NLL_scale.Fill(pmt_id, nll_ar40)
         #HpmtID_surfTF2_NLL_scale.Fill(pmt_id, nll_surf)

    new_NLL_ar40_min_TF2 = new_NLL_ar40 - min([new_NLL_surf, min(NLL_pmtsPDF), min(NLL_gapsPDF)])
    #print "sum ar40/surf NLL", qPE, sumAr40TF2_qBar[0], sumAr40TF2_NLL[0], NLL_ar40PDF_TF2, sumSurfTF2_qBar[0], sumSurfTF2_NLL[0], NLL_surfPDF_TF2
    
    #print NLL_ar40_min_TF2, "new scaled what happened", new_NLL_ar40_min_TF2

    HdeadPMTtoAr40.Fill(NLL_ar40PDF_TF2 - min(NLL_pmtsPDF))

    Hsave_pmtsPDF.Fill(min(NLL_pmtsPDF))
    Hsave_gapsPDF.Fill(min(NLL_gapsPDF))

    Hsave_ar40PDF_MBL.Fill(NLL_ar40PDF_MBL)
    Hsave_surfPDF_MBL.Fill(NLL_surfPDF_MBL)
    Hsave_ar40_min_MBL.Fill(NLL_ar40_min_MBL)

    Hsave_ar40PDF_TF2.Fill(NLL_ar40PDF_TF2)
    Hsave_surfPDF_TF2.Fill(NLL_surfPDF_TF2)
    Hsave_ar40_min_TF2.Fill(NLL_ar40_min_TF2)

    HsaveNew_ar40_min_TF2.Fill(new_NLL_ar40_min_TF2)
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
         Hsort_costhetaMC.SetBinContent(i+1, HcosthetaMC.GetBinContent(255-i))
         Hsort_costhetaMC_qnscb.SetBinContent(i+1, HcosthetaMC_qnscb.GetBinContent(255-i))
         #Hsort_costhetaMC_qprompt.SetBinContent(i+1, HcosthetaMC_qprompt.GetBinContent(255-i))

         Hsort_costheta_MBL.SetBinContent(i+1, Hcostheta_MBL.GetBinContent(255-i))
         Hsort_costheta_TF2.SetBinContent(i+1, Hcostheta_TF2.GetBinContent(255-i))
         Hsort_costheta_deadPMT.SetBinContent(i+1, Hcostheta_deadPMT.GetBinContent(255-i))
    
    ### save to tree2
    runIDval[0] = tree1.runID
    subrunIDval[0] = tree1.subrunID
    evtID[0] = tree1.eventID
    mbx[0] = tree1.mbx
    mby[0] = tree1.mby
    mbz[0] = tree1.mbz
    mcx[0] = xmc
    mcy[0] = ymc
    mcz[0] = zmc
    qpeVal[0] = qPE
    nPMTs[0] = nPMT
    fpromptVal[0] = tree1.fprompt
    nSCBayesVal[0] = nSCB 
    rprompt60BayesVal[0] = tree1.rprompt60Bayes

    #print sumDead149_qBar[0], sumDead149_NLL[0]
    #print sumDead204_qBar[0], sumDead204_NLL[0]
    #print sumSurroundDead149_qBar[0], sumSurroundDead149_NLL[0]
    #print sumSurroundDead204_qBar[0], sumSurroundDead204_NLL[0]
    #print mbx[0], mby[0]
    tree2.Fill()

fout.cd()
tree2.Write()
#for hist in ListHDump_NLL_deadPMT:
#    hist.Write()

#for hist in ListHDump_Q_deadPMT:
#    hist.Write()

#for hist in ListHDump_pdf_deadPMT:
#    hist.Write()

#for hist in ListHDumpPMTid_NLL_deadPMT:
#    hist.Write()

#for hist in ListHDumpPMTid_pdf_deadPMT:
#    hist.Write()

Hsort_costhetaMC.Scale(1./Hsort_costhetaMC.Integral())
Hsort_costhetaMC_qnscb.Scale(1./Hsort_costhetaMC_qnscb.Integral())
#Hsort_costhetaMC_qprompt.Scale(1./Hsort_costhetaMC_qprompt.Integral())

Hsort_costhetaMC.Write()
Hsort_costhetaMC_qnscb.Write()

#Hsort_costhetaMC_qprompt.Write()
## save surround PMTs
#HqBar_surround_dead149.Write()
#HqBar_surround_dead204.Write()
#HNLLsurround_dead149.Write()
#HNLLsurround_dead204.Write()
#
#HqBar_surround_ar40TF2.Write()
#HNLLsurround_ar40TF2.Write()
#HNLLsurround_surfTF2.Write()
#
#HqBar_surround_ar40MBL.Write()
#HNLLsurround_ar40MBL.Write()
#HNLLsurround_surfMBL.Write()

fout.Close()
