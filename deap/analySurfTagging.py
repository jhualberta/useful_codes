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


# load all the required PDFs
Hsurf_truth_angle = filePDF.Get("surf_truth_angle_qPE")
Har40_truth_angle = filePDF.Get("ar40_truth_angle_qPE")
Hsurf_PMTs_angle  = filePDF.Get("surf_PMTs_angle_qPE")
Hsurf_gap10_angle = filePDF.Get("surf_gap10_angle_qPE")
Hsurf_gaps_angle  = filePDF.Get("surf_gaps_angle_qPE")

# load offline pmt positions
ff = file(pmtpos,'r')
pmtpos_offline = []
xpos = []
ypos = []
zpos = []

dead_pmt_id = [149, 204]
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
    ## or why don't we just use gammaln(x)
####################################################################################################################################################
def gof_calculator(xbar, pmt_charge):
    return (xbar - pmt_charge)**2/xbar

# Also sort the PMT indices since they are independent of the fitted event position
pmts_angles_ids, gaps_angles_ids = pmts_gaps_ids()
#-------------------------------------------------------------------------------------------------------------------------------------------------------
# to hold MBL and TF2 fitted positions
MBLEvent = TVector3()
TF2Event = TVector3()

File = TFile(file0);
data = File.Get("T");
nentries = data.GetEntries();


######## 
## skim and save surface tagging variables
## ntupleSurfTag = TNtuple("ntupleSurfTag","surface tag variables","eventID:NLL_surfPDF_MBL:NLL_ar40PDF_MBL:NLL_surfPDF_TF2:NLL_ar40PDF_TF2:min(NLL_pmtsPDF):min(NLL_gapsPDF):gof_surfPDF_MBL:gof_ar40PDF_MBL:gof_surfPDF_TF2:gof_ar40PDF_TF2")

#runIDval = array('l',[0])
#subrunIDval = array('l',[0])
#evtID = array('l',[0])
#mbx = array('f',[0])
#mby = array('f',[0])
#mbz = array('f',[0])
#
#tf2x = array('f',[0])
#tf2y = array('f',[0])
#tf2z = array('f',[0])
#
#mcx = array('f', [0])
#mcy = array('f', [0])
#mcz = array('f', [0])
#
#qpeVal = array('f',[0]) #unsigned double 
#fpromptVal = array('f',[0])
#nSCBayesVal = array('f',[0])
#rprompt60BayesVal = array('f',[0])
#fmaxpeVal = array('f',[0])
#nhitVal = array('i',[0])
#pulseGar = array('f',[0])
#cft2r = array('f',[0]) ## (chargetopring + chargesecondring)/qpe < 0.04
#cfb3r = array('f',[0]) ## (chargebottomring + chargesecondbottomring + chargethirdbottomring)/qpe<0.1
#neckVetoVal = array('f',[0])
#nPMTs1 = array('I',[0])

tree1 = TTree("T1","saveEvent")
update = int(nentries/10000);
#pmtInfo = RAT.PMTInfoUtil.GetPMTInfoUtil()

#tree1.Branch("eventID",evtID) #ULong64_
#tree1.Branch("qpe", qpeVal)
#tree1.Branch("fprompt", fpromptVal)
#tree1.Branch("mbx", mbx)
#tree1.Branch("mby", mby)
#tree1.Branch("mbz", mbz)
#
#tree1.Branch("mcx", mcx)
#tree1.Branch("mcy", mcy)
#tree1.Branch("mcz", mcz)
#
#tree1.Branch("tf2x", tf2x)
#tree1.Branch("tf2y", tf2y)
#tree1.Branch("tf2z", tf2z)
#tree1.Branch("nPMTs1",nPMTs1)
#tree1.Branch("pmtPhi", pmtPhi)
#tree1.Branch("pmtCosTheta", pmtCosTheta)
#tree1.Branch("pmtCharge", pmtCharge)
#tree1.Branch("listPmtID", listPmtID)






for event in range(nentries):
    if (event%20000 == 0): print "processed", event
    data.GetEntry(event)
    DS = data.ds;#(required)
    if DS.GetTSCount() > 0 and DS.GetEVCount() > 0 and DS.GetCALCount() > 0:
    TS = DS.GetTS(0); EV = DS.GetEV(0); cal = DS.GetCAL(0);
    try: CalTrigTime = CAL.GetEventTime()
    eventID = TS.GetEventID();

    ## pre-cuts
    if (cal.GetQPE() <= 0 or cal.GetSubeventCount() > 1 or cal.GetEventTime() <= 2250 or cal.GetEventTime() >= 2700 or cal.GetNumEarlyPulses() > 3):
       continue ## no tagging process

    #runID = DS.GetRunID();
    #subrunID = DS.GetSubrunID();
    evtx = 0; evty = 0; evtz = 0;
    tf2x = 0; tf2y = 0; tf2z = 0;
    ### MB likelihood valid
    try:
        evtx = EV.mblikelihood.GetPosition().X(); evty = EV.mblikelihood.GetPosition().Y(); evtz = EV.mblikelihood.GetPosition().Z();
    except:
        continue

     ### timefit2 valid (fmaxpe<0.75 automatically applied)
     try:
         #time_tf2_ch = EV.timefit2.ch_t
         tf2x = EV.GetTimeFit2().GetX(); tf2y = EV.GetTimeFit2().GetY(); tf2z = EV.GetTimeFit2().GetZ();
     except:
         ## print "timefit2 not valid"
         continue

    #evtx = EV.mblikelihood.GetPosition().X(); evty = EV.mblikelihood.GetPosition().Y(); evtz = EV.mblikelihood.GetPosition().Z();
    #posPhi = EV.mblikelihood.GetPosition().Phi(); posCosTheta = EV.mblikelihood.GetPosition().CosTheta();
    
    MBLEvent.SetXYZ(evtx, evty, evtz)
    TF2Event.SetXYZ(tf2x, tf2y, tf2z)

    TF2_fit_found = 1 if TF2Event.Mag()<855 else 0
    MBL_fit_found = 1 if MBLEvent.Mag()<855 else 0
    if (not TF2_fit_found) and (not MBL_fit_found): continue


    if MBLEvent.Mag()>800: continue

    #------------------------------------------------------------------------
    # sort PMT IDs!
    # NOTE 1: The sorting that is done within the *pmts_gaps_ids* function, sort the PMTs based on the pentagonal gap positions
    # and dead PMTS (if they exist)! This block here sorts the PMTs based on the MBL and TF2 positions.
    # NOTE 2: The cos theta is the angle between two vectors v1 and v2_i where v1 is from the centre of the detector to the event position and v2_i is from the centre to PMT_i.
    # Pentagonal gaps and dead PMTs are similar, where v1 is the vector from the centre to the position of the gap/dead PMT

    MBL_angles, TF2_angles = {}, {}

    # first we need all the 255 PMTs to sort them based on cos_theta!
    temp_pmt_pos = TVector3()
    for iPMT in range(255):
        temp_pmt_pos.SetXYZ( pmtpos_offline[iPMT][0], pmtpos_offline[iPMT][1], pmtpos_offline[iPMT][2] )
        temp_pmt_pos.SetMag(851)

        MBL_angles[iPMT] = np.cos(MBLEvent.Angle((temp_pmt_pos)))
        TF2_angles[iPMT] = np.cos(TF2Event.Angle((temp_pmt_pos)))

    MBL_angles_dict = sorted(MBL_angles.items(), key=operator.itemgetter(1), reverse=True)
    MBL_angles_ids  = [item[0] for item in MBL_angles_dict]

    TF2_angles_dict = sorted(TF2_angles.items(), key=operator.itemgetter(1), reverse=True)
    TF2_angles_ids  = [item[0] for item in TF2_angles_dict]
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
    qPE = qpeVal[0]  
    nPMT = nPMTs1[0]
    for iPMT in range(255):
        if iPMT==nPMT:
            # Complete the list of PMT IDs with the PMTs that did not see any charge in this event.
            # Since we are using Poisson distribution, PMTs not seeing charge can contain information for the likelihood.
            pmt_ids += [i for i in range(255) if i not in pmt_ids]

        if iPMT<nPMT:
            # NOTE: NLL calculation needs integer expected charge; but not the gof calculation.
            #pmt = cal.GetPMT(iPMT)
            #pmt = cal.GetPMT(iPMT)
            #pmt_id = pmt.GetID()
            #pmt_qPE = pmt.qPE

            pmt_id = listPmtID[iPMT] 
            pmt_qPE = pmtCharge[iPMT] 
            pmt_ids.append(pmt_id)
        else:
            # PMTs that have not observed any charge.
            pmt_id = pmt_ids[iPMT]
            pmt_qPE = 0

        #+1 is needed for GetBinContent because index 0 is underflow.
        #------------------------------------------------------------------------
        # calculate NLL and gof based on MBL position
        if MBL_fit_found:
            pdf_index = MBL_angles_ids.index(pmt_id)+1
            xbar_surfPDF = Hsurf_truth_angle.GetBinContent(pdf_index)*qPE
            xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
            NLL_surfPDF_MBL += get_poisson_nll(xbar_surfPDF, int(round(pmt_qPE)))
            NLL_ar40PDF_MBL += get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))

            gof_surfPDF_MBL += gof_calculator (xbar_surfPDF, pmt_qPE)
            gof_ar40PDF_MBL += gof_calculator (xbar_ar40PDF, pmt_qPE)
        #------------------------------------------------------------------------
        # calculate NLL and gof based on TF2 position
        if TF2_fit_found:
            pdf_index = TF2_angles_ids.index(pmt_id)+1
            xbar_surfPDF = Hsurf_truth_angle.GetBinContent(pdf_index)*qPE
            xbar_ar40PDF = Har40_truth_angle.GetBinContent(pdf_index)*qPE
            NLL_surfPDF_TF2 += get_poisson_nll(xbar_surfPDF, int(round(pmt_qPE)))
            NLL_ar40PDF_TF2 += get_poisson_nll(xbar_ar40PDF, int(round(pmt_qPE)))

            gof_surfPDF_TF2 += gof_calculator (xbar_surfPDF, pmt_qPE)
            gof_ar40PDF_TF2 += gof_calculator (xbar_ar40PDF, pmt_qPE)
        #------------------------------------------------------------------------
        # calculate NLL and gof based on dead PMT position(s)
        for i in range(len(pmts_angles_ids)):
            pdf_index = pmts_angles_ids[i].index(pmt_id)+1
            overwrite = True if pmt_id in dead_pmt_id else False # avoid X_deapPMT*X_deadPMT
            if not overwrite: ## pmt_id not in dead pmts
                xbar_deadPDF = Hsurf_PMTs_angle.GetBinContent(pdf_index)*qPE
            else:
                # If a (second) dead PMT, expection would be also zero charge.
                xbar_deadPDF = Hsurf_PMTs_angle.GetBinContent(1)*qPE
            NLL_pmtsPDF[i] += get_poisson_nll(xbar_deadPDF, int(round(pmt_qPE)))
            gof_pmtsPDF_PMTs[i] += gof_calculator(xbar_deadPDF, pmt_qPE)
        #------------------------------------------------------------------------
        # calculate NLL and gof based on pentagonal gap positions
        for i in range(len(gaps_angles_ids)):
            pdf_index = gaps_angles_ids[i].index(pmt_id)+1
            xbar_gapsPDF = Hsurf_gaps_angle.GetBinContent(pdf_index)*qPE
            NLL_gapsPDF[i] += get_poisson_nll(xbar_gapsPDF, int(round(pmt_qPE)))
            gof_gapsPDF_gaps[i] += gof_calculator(xbar_gapsPDF, pmt_qPE)
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
    #------------------------------------------------------------------------
    # Set main NLL variables
    NLL_tag_MBL_list = [NLL_ar40PDF_MBL, NLL_surfPDF_MBL, min(NLL_pmtsPDF), min(NLL_gapsPDF)]
    NLL_tag_TF2_list = [NLL_ar40PDF_TF2, NLL_surfPDF_TF2, min(NLL_pmtsPDF), min(NLL_gapsPDF)]

    is_surface_like_MBL = 0 if min(NLL_tag_MBL_list)==NLL_tag_MBL_list[0] else 1
    is_surface_like_TF2 = 0 if min(NLL_tag_TF2_list)==NLL_tag_TF2_list[0] else 1

    # Mainly for events tag as ar40 --> How far/close was the next minimum NLL?!

    NLL_ar40_min_MBL = NLL_ar40PDF_MBL - min(NLL_tag_MBL_list[1:])
    NLL_ar40_min_TF2 = NLL_ar40PDF_TF2 - min(NLL_tag_TF2_list[1:])

    #------------------------------------------------------------------------
    # Set gof variables and the the two is_surface tags (one for MBL, one for TF2)
    #ev.Setgof_pmtsPDF    (min(gof_pmtsPDF_PMTs))
    #ev.Setgof_gapsPDF    (min(gof_gapsPDF_gaps))

    #ev.SetIs_surface_like_MBL(is_surface_like_MBL)
    #ev.SetIs_surface_like_TF2(is_surface_like_TF2)
    #print NLL_ar40_min_MBL, NLL_ar40_min_TF2

    print "NLL_ar40_min MBL", round(NLL_ar40PDF_MBL, 3), "-", NLL_tag_MBL_list[1:], "=", round(NLL_ar40_min_MBL,3)
    print "NLL_ar40_min_TF2", round(NLL_ar40PDF_TF2, 3), "-", NLL_tag_TF2_list[1:], "=", round(NLL_ar40_min_TF2,3)

    #ntupleSurfTag.Fill(eventID:NLL_surfPDF_MBL:NLL_ar40PDF_MBL:NLL_surfPDF_TF2:NLL_ar40PDF_TF2:min(NLL_pmtsPDF):min(NLL_gapsPDF):gof_surfPDF_MBL:gof_ar40PDF_MBL:gof_surfPDF_TF2:gof_ar40PDF_TF2)

########################################################################################################################################################
