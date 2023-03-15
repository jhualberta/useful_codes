import sys
import glob
from rat import *
from ROOT import *
import numpy as np
from numpy import sqrt
import struct
#from run_numbers import runs
from datetime import datetime, time as datetime_time, timedelta
## Sina: expanding MBLR cut, use toproi05 instead of top30, and remove pulseindexfirstGar!
ffroi = TFile("saveSideBandROI.root")
#### saveSideBandROI.root")
####roi_802_days_11March2020_nsc_rp60.root")
## from loose cut to tight cut
roicut_top05 = ffroi.Get("top05")
roicut_top30 = ffroi.Get("top30")
roicut_top55 = ffroi.Get("top55")

sidebandcut = ffroi.Get("cutSideBand")
sidebandcut.SetName("cutSideBand")

#runs_name = ["U232Source_SLB005_2020_CalC_L0", "U232Source_SLB005_2020_CalC_RT465_L0", "U232Source_SLB005_2020_CalC_RT480_L0", "U232Source_SLB005_2020_CalF_L0", "U232Source_SLB005_2020_CalFpos6_RT465_L0", "U232Source_SLB005_2020_CalFpos6_RT480_L0", "U232Source_SLB005_2020_CalFpos7_RT465_L0", "U232Source_SLB005_2020_CalFpos7_RT480_L0", "U232Source_SLB005_2021_CalFpos6_RT480_L0"]
#
#RunList = runs

#Run            = sys.argv[1]
file0 = str(sys.argv[1])
fileName = os.path.basename(file0)
### 7-level cuts
###############################################################################################################################
H2_nSCBayes_rprompt60Bayes_0 = TH2D( "H2_nSCBayes_rprompt60Bayes_0", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_0             = TH2D( "H2_qpe_fprompt_0", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_0 = TH2D("H2_rhoZ_0", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_0 = TH2D("H2_nhitQPE_0", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_0 = TH2D("H2_nhitNSCB_0", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_0 = TH1D("H1_fmaxpe_0","fmaxpe",100,0,1)
H1_scintlike_0 = TH1D("H1_scintlike_0","scintlike",1000,0,1000)
H1_llneutron_0 = TH1D("H1_llneutron_0","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_1 = TH2D("H2_nSCBayes_rprompt60Bayes_1", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_1             = TH2D("H2_qpe_fprompt_1", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_1 = TH2D("H2_rhoZ_1", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_1 = TH2D("H2_nhitQPE_1", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_1 = TH2D("H2_nhitNSCB_1", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_1 = TH1D("H1_fmaxpe_1","fmaxpe",100,0,1)
H1_scintlike_1 = TH1D("H1_scintlike_1","scintlike",1000,0,1000)
H1_llneutron_1 = TH1D("H1_llneutron_1","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_2 = TH2D("H2_nSCBayes_rprompt60Bayes_2", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_2             = TH2D("H2_qpe_fprompt_2", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_2 = TH2D("H2_rhoZ_2", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_2 = TH2D("H2_nhitQPE_2", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_2 = TH2D("H2_nhitNSCB_2", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_2 = TH1D("H1_fmaxpe_2","fmaxpe",100,0,1)
H1_scintlike_2 = TH1D("H1_scintlike_2","scintlike",1000,0,1000)
H1_llneutron_2 = TH1D("H1_llneutron_2","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_3 = TH2D("H2_nSCBayes_rprompt60Bayes_3", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_3             = TH2D("H2_qpe_fprompt_3", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_3 = TH2D("H2_rhoZ_3", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_3 = TH2D("H2_nhitQPE_3", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_3 = TH2D("H2_nhitNSCB_3", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_3 = TH1D("H1_fmaxpe_3","fmaxpe",100,0,1)
H1_scintlike_3 = TH1D("H1_scintlike_3","scintlike",1000,0,1000)
H1_llneutron_3 = TH1D("H1_llneutron_3","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_4 = TH2D("H2_nSCBayes_rprompt60Bayes_4", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_4             = TH2D("H2_qpe_fprompt_4", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_4 = TH2D("H2_rhoZ_4", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_4 = TH2D("H2_nhitQPE_4", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_4 = TH2D("H2_nhitNSCB_4", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_4 = TH1D("H1_fmaxpe_4","fmaxpe",100,0,1)
H1_scintlike_4 = TH1D("H1_scintlike_4","scintlike",1000,0,1000)
H1_llneutron_4 = TH1D("H1_llneutron_4","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_5 = TH2D("H2_nSCBayes_rprompt60Bayes_5", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_5             = TH2D("H2_qpe_fprompt_5", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_5 = TH2D("H2_rhoZ_5", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_5 = TH2D("H2_nhitQPE_5", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_5 = TH2D("H2_nhitNSCB_5", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_5 = TH1D("H1_fmaxpe_5","fmaxpe",100,0,1)
H1_scintlike_5 = TH1D("H1_scintlike_5","scintlike",1000,0,1000)
H1_llneutron_5 = TH1D("H1_llneutron_5","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_6 = TH2D("H2_nSCBayes_rprompt60Bayes_6", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_6             = TH2D("H2_qpe_fprompt_6", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_6 = TH2D("H2_rhoZ_6", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_6 = TH2D("H2_nhitQPE_6", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_6 = TH2D("H2_nhitNSCB_6", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_6 = TH1D("H1_fmaxpe_6","fmaxpe",100,0,1)
H1_scintlike_6 = TH1D("H1_scintlike_6","scintlike",1000,0,1000)
H1_llneutron_6 = TH1D("H1_llneutron_6","llneutron",200,-100,100)

H2_nSCBayes_rprompt60Bayes_7 = TH2D("H2_nSCBayes_rprompt60Bayes_7", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7             = TH2D("H2_qpe_fprompt_7", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7 = TH2D("H2_rhoZ_7", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7 = TH2D("H2_nhitQPE_7", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7 = TH2D("H2_nhitNSCB_7", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7 = TH1D("H1_fmaxpe_7","fmaxpe",100,0,1)
H1_scintlike_7 = TH1D("H1_scintlike_7","scintlike",1000,0,1000)
H1_llneutron_7 = TH1D("H1_llneutron_7","llneutron",200,-100,100)

###########################################################################################################################
### Level 8, 7-level cuts + STR levels, start with fmaxpe first !!
H2_nSCBayes_rprompt60Bayes_7_fmaxpe = TH2D( "H2_nSCBayes_rprompt60Bayes_7_fmaxpe", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_fmaxpe = TH2D( "H2_qpe_fprompt_7_fmaxpe", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_fmaxpe = TH2D("H2_rhoZ_7_fmaxpe", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_fmaxpe = TH2D("H2_nhitQPE_7_fmaxpe", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_fmaxpe = TH2D("H2_nhitNSCB_7_fmaxpe", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_fmaxpe = TH1D("H1_fmaxpe_7_fmaxpe","fmaxpe",100,0,1)
H1_scintlike_7_fmaxpe = TH1D("H1_scintlike_7_fmaxpe","scintlike",1000,0,1000)
H1_llneutron_7_fmaxpe = TH1D("H1_llneutron_7_fmaxpe","llneutron",200,-100,100)
## Level 9, fmaxpe + neckVeto
H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck = TH2D( "H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_fmaxpe_neck = TH2D( "H2_qpe_fprompt_7_fmaxpe_neck", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_fmaxpe_neck = TH2D("H2_rhoZ_7_fmaxpe_neck", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_fmaxpe_neck = TH2D("H2_nhitQPE_7_fmaxpe_neck", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_fmaxpe_neck = TH2D("H2_nhitNSCB_7_fmaxpe_neck", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_fmaxpe_neck = TH1D("H1_fmaxpe_7_fmaxpe_neck","fmaxpe",100,0,1)
H1_scintlike_7_fmaxpe_neck = TH1D("H1_scintlike_7_fmaxpe_neck","scintlike",1000,0,1000)
H1_llneutron_7_fmaxpe_neck = TH1D("H1_llneutron_7_fmaxpe_neck","llneutron",200,-100,100)
## Level 10, fmaxpe + neckVeto + mbR, mbR is enlarged to 800 mm
H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck_mbR = TH2D( "H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck_mbR", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_fmaxpe_neck_mbR = TH2D( "H2_qpe_fprompt_7_fmaxpe_neck_mbR", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_fmaxpe_neck_mbR = TH2D("H2_rhoZ_7_fmaxpe_neck_mbR", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_fmaxpe_neck_mbR = TH2D("H2_nhitQPE_7_fmaxpe_neck_mbR", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_fmaxpe_neck_mbR = TH2D("H2_nhitNSCB_7_fmaxpe_neck_mbR", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_fmaxpe_neck_mbR = TH1D("H1_fmaxpe_7_fmaxpe_neck_mbR","fmaxpe", 100,0,1)
H1_scintlike_7_fmaxpe_neck_mbR = TH1D("H1_scintlike_7_fmaxpe_neck_mbR","scintlike",1000,0,1000)
H1_llneutron_7_fmaxpe_neck_mbR = TH1D("H1_llneutron_7_fmaxpe_neck_mbR","llneutron",200,-100,100)

###Level 11, 7-level cuts + STR(fmaxpe+neckVeto+mbR) + start LAr with pulseG
## H2_nSCBayes_rprompt60Bayes_7_str_pulseG = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_pulseG", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_pulseG = TH2D( "H2_qpe_fprompt_7_str_pulseG", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_pulseG = TH2D("H2_rhoZ_7_str_pulseG", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_pulseG = TH2D("H2_nhitQPE_7_str_pulseG", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_pulseG = TH2D("H2_nhitNSCB_7_str_pulseG", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_pulseG = TH1D("H1_fmaxpe_7_str_pulseG","fmaxpe",100,0,1)
## H1_scintlike_7_str_pulseG = TH1D("H1_scintlike_7_str_pulseG","scintlike",1000,0,1000)
## H1_llneutron_7_str_pulseG = TH1D("H1_llneutron_7_str_pulseG","llneutron",200,-100,100)

### Level 12, LAr: pulse G + cft2r_mbZ 
## H2_nSCBayes_rprompt60Bayes_7_str_pulseG_cft2r_mbZ = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_pulseG_cft2r_mbZ", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_pulseG_cft2r_mbZ = TH2D( "H2_qpe_fprompt_7_str_pulseG_cft2r_mbZ", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_pulseG_cft2r_mbZ = TH2D("H2_rhoZ_7_str_pulseG_cft2r_mbZ", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_pulseG_cft2r_mbZ = TH2D("H2_nhitQPE_7_str_pulseG_cft2r_mbZ", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_pulseG_cft2r_mbZ = TH2D("H2_nhitNSCB_7_str_pulseG_cft2r_mbZ", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_pulseG_cft2r_mbZ = TH1D("H1_fmaxpe_7_str_pulseG_cft2r_mbZ","fmaxpe",100,0,1)
## H1_scintlike_7_str_pulseG_cft2r_mbZ = TH1D("H1_scintlike_7_str_pulseG_cft2r_mbZ","scintlike",1000,0,1000)
## H1_llneutron_7_str_pulseG_cft2r_mbZ = TH1D("H1_llneutron_7_str_pulseG_cft2r_mbZ","llneutron",200,-100,100)

### Level 13, LAr: pulse G + cft2r_mbZ + cfb3r
## H2_nSCBayes_rprompt60Bayes_7_str_pulseG_cft2r_mbZ_cfb3r = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_pulseG_cft2r_mbZ_cfb3r", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_pulseG_cft2r_mbZ_cfb3r = TH2D( "H2_qpe_fprompt_7_str_pulseG_cft2r_mbZ_cfb3r", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_pulseG_cft2r_mbZ_cfb3r = TH2D("H2_rhoZ_7_str_pulseG_cft2r_mbZ_cfb3r", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_pulseG_cft2r_mbZ_cfb3r = TH2D("H2_nhitQPE_7_str_pulseG_cft2r_mbZ_cfb3r", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_pulseG_cft2r_mbZ_cfb3r = TH2D("H2_nhitNSCB_7_str_pulseG_cft2r_mbZ_cfb3r", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_pulseG_cft2r_mbZ_cfb3r = TH1D("H1_fmaxpe_7_str_pulseG_cft2r_mbZ_cfb3r","fmaxpe",100,0,1)
## H1_scintlike_7_str_pulseG_cft2r_mbZ_cfb3r = TH1D("H1_scintlike_7_str_pulseG_cft2r_mbZ_cfb3r","scintlike",1000,0,1000)
## H1_llneutron_7_str_pulseG_cft2r_mbZ_cfb3r = TH1D("H1_llneutron_7_str_pulseG_cft2r_mbZ_cfb3r","llneutron", 200,-100,100)
## 
## #### !!!! Be careful that TF2 must exist then it is fair to do this cut
## #### TF2 - MB consistency cuts, Level 14: z consistency
## H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_LAr_zCon = TH2D( "H2_qpe_fprompt_7_str_LAr_zCon", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_LAr_zCon = TH2D("H2_rhoZ_7_str_LAr_zCon", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_LAr_zCon = TH2D("H2_nhitQPE_7_str_LAr_zCon", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_LAr_zCon = TH2D("H2_nhitNSCB_7_str_LAr_zCon", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_LAr_zCon = TH1D("H1_fmaxpe_7_str_LAr_zCon","fmaxpe",100,0,1)
## H1_scintlike_7_str_LAr_zCon = TH1D("H1_scintlike_7_str_LAr_zCon","scintlike",1000,0,1000)
## H1_llneutron_7_str_LAr_zCon = TH1D("H1_llneutron_7_str_LAr_zCon","llneutron", 200,-100,100)
## 
## ### TF2 - MB consistency cuts, Level 15
## H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_LAr_zrCon = TH2D( "H2_qpe_fprompt_7_str_LAr_zrCon", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_LAr_zrCon = TH2D("H2_rhoZ_7_str_LAr_zrCon", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_LAr_zrCon = TH2D("H2_nhitQPE_7_str_LAr_zrCon", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_LAr_zrCon = TH2D("H2_nhitNSCB_7_str_LAr_zrCon", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_LAr_zrCon = TH1D("H1_fmaxpe_7_str_LAr_zrCon","fmaxpe",100,0,1)
## H1_scintlike_7_str_LAr_zrCon = TH1D("H1_scintlike_7_str_LAr_zrCon","scintlike",1000,0,1000)
## H1_llneutron_7_str_LAr_zrCon = TH1D("H1_llneutron_7_str_LAr_zrCon","llneutron", 200,-100,100)
## 
## #### TF2 - MB consistency cuts, Level 14 a or level 16: z consistency
## H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon_tf2Cerenkov = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon_tf2Cerenkov", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_LAr_zCon_tf2Cerenkov = TH2D( "H2_qpe_fprompt_7_str_LAr_zCon_tf2Cerenkov", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_LAr_zCon_tf2Cerenkov = TH2D("H2_rhoZ_7_str_LAr_zCon_tf2Cerenkov", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_LAr_zCon_tf2Cerenkov = TH2D("H2_nhitQPE_7_str_LAr_zCon_tf2Cerenkov", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_LAr_zCon_tf2Cerenkov = TH2D("H2_nhitNSCB_7_str_LAr_zCon_tf2Cerenkov", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_LAr_zCon_tf2Cerenkov = TH1D("H1_fmaxpe_7_str_LAr_zCon_tf2Cerenkov","fmaxpe",100,0,1)
## H1_scintlike_7_str_LAr_zCon_tf2Cerenkov = TH1D("H1_scintlike_7_str_LAr_zCon_tf2Cerenkov","scintlike",1000,0,1000)
## H1_llneutron_7_str_LAr_zCon_tf2Cerenkov = TH1D("H1_llneutron_7_str_LAr_zCon_tf2Cerenkov","llneutron", 200,-100,100)
## 
## ### TF2 - MB consistency cuts, Level 15 a or level 17
## H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon_tf2Cerenkov = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon_tf2Cerenkov", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_qpe_fprompt_7_str_LAr_zrCon_tf2Cerenkov = TH2D( "H2_qpe_fprompt_7_str_LAr_zrCon_tf2Cerenkov", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
## H2_rhoZ_7_str_LAr_zrCon_tf2Cerenkov = TH2D("H2_rhoZ_7_str_LAr_zrCon_tf2Cerenkov", "rho vs z", 850, 0, 850, 1700, -850, 850)
## H2_nhitQPE_7_str_LAr_zrCon_tf2Cerenkov = TH2D("H2_nhitQPE_7_str_LAr_zrCon_tf2Cerenkov", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
## H2_nhitNSCB_7_str_LAr_zrCon_tf2Cerenkov = TH2D("H2_nhitNSCB_7_str_LAr_zrCon_tf2Cerenkov", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
## H1_fmaxpe_7_str_LAr_zrCon_tf2Cerenkov = TH1D("H1_fmaxpe_7_str_LAr_zrCon_tf2Cerenkov","fmaxpe",100,0,1)
## H1_scintlike_7_str_LAr_zrCon_tf2Cerenkov = TH1D("H1_scintlike_7_str_LAr_zrCon_tf2Cerenkov","scintlike",1000,0,1000)
## H1_llneutron_7_str_LAr_zrCon_tf2Cerenkov = TH1D("H1_llneutron_7_str_LAr_zrCon_tf2Cerenkov","llneutron", 200,-100,100)

#Files = AllFiles[subrun];
File = TFile(file0)
print ("\n********************************** File " , File , " loaded... **********************************\n")
data = File.Get("data_satCorr");#TTree
nentries = data.GetEntries();
update = int(nentries/10);
fout = TFile("Sina_"+fileName, "RECREATE")
###########################################################################################################################
#foutfile="/project/6004969/dpapi/cherenkov/skim/{}/{}".format(runs_name[int(Run)-1], value.split("/")[-1])
#fout = TFile(foutfile,"RECREATE")
dstreeclone = data.CloneTree(0);
dstreeclone.SetDirectory(fout);

level = 15 
countTotal = [0 for i in range(level)]
countROI05 = [0 for i in range(level)]
countROI30 = [0 for i in range(level)]
countROI55 = [0 for i in range(level)]
countSideband = [0 for i in range(level)]

H_countTotal = TH1F("H_countTotal", "total", level, 0, level)
H_countROI05 = TH1F("H_countROI05", "top05", level, 0, level)
H_countROI30 = TH1F("H_countROI30", "top30", level, 0, level)
H_countROI55 = TH1F("H_countROI55", "top55", level, 0, level)
H_countSideBand = TH1F("H_countSideBand", "sideband cuts", level, 0, level)

### this ntuple saves the event info for tf2 invalid
### type == 0 ->timefit2,  type ==1 ->timefit2_cerenkov
ntupleTF2 = TNtuple("ntupleTF2","save TF2 info","runID:subrunID:eventID:nscb:rprompt60Bayes:qpe:fprompt:mbX:mbY:mbZ:tf2X:tf2Y:tf2Z:pulseGar")

###########################################################################################################################

for event in range(nentries):
    if (event+1)%update == 0:
        print (event+1), "Analyzed..."
    data.GetEntry(event)
    flag=0
    runID = data.runID;
    subrunID = data.subrunID;
    eventID = data.eventID;
    nscb = data.nSCBayes; 
    qpe = data.qPE;
    ## print qpe
    fprompt = data.fprompt; rprompt60Bayes = data.rprompt60Bayes;
    nhit00 = data.nhit
    nhit = ord(nhit00) #### Nhit issue!!!!
    ### print "nhitsss", nhit 
    fmaxpe = data.fmaxpe
    llneutronflash = data.llneutronflash
    scintlike = data.scintlike
    mbX = data.mblikelihoodX
    mbY = data.mblikelihoodY
    mbZ = data.mblikelihoodZ
    mbR = data.mblikelihoodR
    H2_nSCBayes_rprompt60Bayes_0.Fill(nscb, rprompt60Bayes)
    H2_qpe_fprompt_0.Fill(qpe, fprompt)
    H2_rhoZ_0.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
    if qpe>0: H2_nhitQPE_0.Fill(qpe, nhit/qpe)
    if nscb>0: H2_nhitNSCB_0.Fill(nscb, nhit/nscb)
    H1_fmaxpe_0.Fill(fmaxpe)
    H1_scintlike_0.Fill(scintlike)
    H1_llneutron_0.Fill(llneutronflash)
    countTotal[0] += 1
    if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[0] += 1
    if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[0] += 1
    if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[0] += 1
    if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[0] += 1

    if data.dtmTrigSrc&0x82 == 0: ## level 1
        H2_nSCBayes_rprompt60Bayes_1.Fill(nscb, rprompt60Bayes)
        H2_qpe_fprompt_1.Fill(qpe, fprompt)
        H2_rhoZ_1.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
        if qpe>0: H2_nhitQPE_1.Fill(qpe, nhit/qpe)
        if nscb>0: H2_nhitNSCB_1.Fill(nscb, nhit/nscb)
        H1_fmaxpe_1.Fill(fmaxpe)
        H1_scintlike_1.Fill(scintlike)
        H1_llneutron_1.Fill(llneutronflash)
        countTotal[1] += 1
        if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[1] += 1
        if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[1] += 1
        if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[1] += 1
        if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[1] += 1

        if data.calcut&0x31f8 == 0: ## level 2
           H2_nSCBayes_rprompt60Bayes_2.Fill(nscb, rprompt60Bayes)
           H2_qpe_fprompt_2.Fill(qpe, fprompt)
           H2_rhoZ_2.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
           if qpe>0: H2_nhitQPE_2.Fill(qpe, nhit/qpe)
           if nscb>0: H2_nhitNSCB_2.Fill(nscb, nhit/nscb)
           H1_fmaxpe_2.Fill(fmaxpe)
           H1_scintlike_2.Fill(scintlike)
           H1_llneutron_2.Fill(llneutronflash)
           countTotal[2] += 1
           if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[2] += 1
           if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[2] += 1
           if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[2] += 1
           if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[2] += 1

           if data.deltat > 20000: ## level 3
              H2_nSCBayes_rprompt60Bayes_3.Fill(nscb, rprompt60Bayes)
              H2_qpe_fprompt_3.Fill(qpe, fprompt)
              H2_rhoZ_3.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
              if qpe>0: H2_nhitQPE_3.Fill(qpe, nhit/qpe)
              if nscb>0: H2_nhitNSCB_3.Fill(nscb, nhit/nscb)
              H1_fmaxpe_3.Fill(fmaxpe)
              H1_scintlike_3.Fill(scintlike)
              H1_llneutron_3.Fill(llneutronflash)
              countTotal[3] += 1
              if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[3] += 1
              if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[3] += 1
              if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[3] += 1
              if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[3] += 1

              if data.numEarlyPulses <= 3: ## level 4
                 H2_nSCBayes_rprompt60Bayes_4.Fill(nscb, rprompt60Bayes)
                 H2_qpe_fprompt_4.Fill(qpe, fprompt)
                 H2_rhoZ_4.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                 if qpe>0: H2_nhitQPE_4.Fill(qpe, nhit/qpe)
                 if nscb>0: H2_nhitNSCB_4.Fill(nscb, nhit/nscb)
                 H1_fmaxpe_4.Fill(fmaxpe)
                 H1_scintlike_4.Fill(scintlike)
                 H1_llneutron_4.Fill(llneutronflash)
                 countTotal[4] += 1
                 if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[4] += 1
                 if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[4] += 1
                 if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[4] += 1
                 if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[4] += 1

                 if data.subeventN==1: ## level 5
                    H2_nSCBayes_rprompt60Bayes_5.Fill(nscb, rprompt60Bayes)
                    H2_qpe_fprompt_5.Fill(qpe, fprompt)
                    H2_rhoZ_5.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                    if qpe>0: H2_nhitQPE_5.Fill(qpe, nhit/qpe)
                    if nscb>0: H2_nhitNSCB_5.Fill(nscb, nhit/nscb)
                    H1_fmaxpe_5.Fill(fmaxpe)
                    H1_scintlike_5.Fill(scintlike)
                    H1_llneutron_5.Fill(llneutronflash)
                    countTotal[5] += 1
                    if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[5] += 1
                    if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[5] += 1
                    if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[5] += 1
                    if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[5] += 1

                    if data.eventTime > 2250  and data.eventTime < 2700: ## level 6
                        H2_nSCBayes_rprompt60Bayes_6.Fill(nscb, rprompt60Bayes)
                        H2_qpe_fprompt_6.Fill(qpe, fprompt)
                        H2_rhoZ_6.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                        if qpe>0: H2_nhitQPE_6.Fill(qpe, nhit/qpe)
                        if nscb>0: H2_nhitNSCB_6.Fill(nscb, nhit/nscb)
                        H1_fmaxpe_6.Fill(fmaxpe)
                        H1_scintlike_6.Fill(scintlike)
                        H1_llneutron_6.Fill(llneutronflash)
                        countTotal[6] += 1
                        if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[6] += 1
                        if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[6] += 1
                        if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[6] += 1
                        if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[6] += 1

                        if qpe>60: ## level 7
                          H2_nSCBayes_rprompt60Bayes_7.Fill(nscb, rprompt60Bayes)
                          H2_qpe_fprompt_7.Fill(qpe, fprompt)
                          H2_rhoZ_7.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                          if qpe>0: H2_nhitQPE_7.Fill(qpe, nhit/qpe)
                          if nscb>0: H2_nhitNSCB_7.Fill(nscb, nhit/nscb)
                          H1_fmaxpe_7.Fill(fmaxpe)
                          H1_scintlike_7.Fill(scintlike)
                          H1_llneutron_7.Fill(llneutronflash)
                          countTotal[7] += 1

                          if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[7] += 1
                          if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[7] += 1
                          if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[7] += 1
                          if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[7] += 1
                          ### Note: for some events, qpe could be 0
                          cft2r = (data.chargetopring + data.chargesecondring)/qpe
                          cfb3r = (data.chargebottomring + data.chargesecondbottomring + data.chargethirdbottomring)/qpe

                          ### Level 8, 7-level cuts + STR levels, start with fmaxpe
                          if data.fmaxpe<0.4:
                            H2_nSCBayes_rprompt60Bayes_7_fmaxpe.Fill(nscb, rprompt60Bayes)
                            H2_qpe_fprompt_7_fmaxpe.Fill(qpe, fprompt)
                            H2_rhoZ_7_fmaxpe.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                            if qpe>0: H2_nhitQPE_7_fmaxpe.Fill(qpe, nhit/qpe)
                            if nscb>0: H2_nhitNSCB_7_fmaxpe.Fill(nscb, nhit/nscb)
                            H1_fmaxpe_7_fmaxpe.Fill(fmaxpe)
                            H1_scintlike_7_fmaxpe.Fill(scintlike)
                            H1_llneutron_7_fmaxpe.Fill(llneutronflash)
                            countTotal[8] += 1
                            if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[8] += 1
                            if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[8] += 1
                            if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[8] += 1
                            if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[8] += 1

                            if data.neckVetoN  == 0: ## level 9
                              H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck.Fill(nscb, rprompt60Bayes)
                              H2_qpe_fprompt_7_fmaxpe_neck.Fill(qpe, fprompt)
                              H2_rhoZ_7_fmaxpe_neck.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                              if qpe>0: H2_nhitQPE_7_fmaxpe_neck.Fill(qpe, nhit/qpe)
                              if nscb>0: H2_nhitNSCB_7_fmaxpe_neck.Fill(nscb, nhit/nscb)
                              H1_fmaxpe_7_fmaxpe_neck.Fill(fmaxpe)
                              H1_scintlike_7_fmaxpe_neck.Fill(scintlike)
                              H1_llneutron_7_fmaxpe_neck.Fill(llneutronflash)
                              countTotal[9] += 1
                              if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[9] += 1
                              if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[9] += 1
                              if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[9] += 1
                              if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[9] += 1

                              if mbR<800: ## level 10 
                                H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck_mbR.Fill(nscb, rprompt60Bayes)
                                H2_qpe_fprompt_7_fmaxpe_neck_mbR.Fill(qpe, fprompt)
                                H2_rhoZ_7_fmaxpe_neck_mbR.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                if qpe>0: H2_nhitQPE_7_fmaxpe_neck_mbR.Fill(qpe, nhit/qpe)
                                if nscb>0: H2_nhitNSCB_7_fmaxpe_neck_mbR.Fill(nscb, nhit/nscb)
                                H1_fmaxpe_7_fmaxpe_neck_mbR.Fill(fmaxpe)
                                H1_scintlike_7_fmaxpe_neck_mbR.Fill(scintlike)
                                H1_llneutron_7_fmaxpe_neck_mbR.Fill(llneutronflash)
                                countTotal[10] += 1
                                if roicut_top05.IsInside(nscb, rprompt60Bayes):countROI05[10] += 1
                                if roicut_top30.IsInside(nscb, rprompt60Bayes):countROI30[10] += 1
                                if roicut_top55.IsInside(nscb, rprompt60Bayes):countROI55[10] += 1
                                if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[10] += 1
                        
                                ######## just go to Level 10 !!! 
                                ## promptPE = nscb*rprompt60Bayes;
                                mbPos = TVector3(mbX,mbY,mbZ);
                                tf2X = data.timefit2X; tf2Y = data.timefit2Y; tf2Z = data.timefit2Z;
                                tf2cerenkovX = data.timefit2cerenkovX; tf2cerenkovY = data.timefit2cerenkovY; tf2cerenkovZ = data.timefit2cerenkovZ;
                                ## just save tf2 info for further checking
                                pulseGar = data.pulseindexfirstgar 
                                ##             runID:subrunID:eventID:nscb:rprompt60Bayes:qpe:fprompt:mbX:mbY:mbZ:tf2X:tf2Y:tf2Z:tf2cerenkovX:tf2cerenkovY:tf2cerenkovZ:pulseGar
                                ntupleTF2.Fill(runID, subrunID, eventID, nscb, rprompt60Bayes, qpe, fprompt, mbX, mbY, mbZ, tf2X, tf2Y, tf2Z, pulseGar)#tf2cerenkovX, tf2cerenkovY, tf2cerenkovZ, pulseGar)


## fill the counts in ROI for each level cuts
for i in range(0, level):
   H_countTotal.Fill(i,countTotal[i])
   H_countROI05.Fill(i, countROI05[i])
   H_countROI30.Fill(i, countROI30[i])
   H_countROI55.Fill(i, countROI55[i])
   H_countSideBand.Fill(i, countSideband[i])

#    if flag:
#        dstreeclone.Fill()
#################################################################################################################################

#dstreeclone.Write()
#fout.Close()
#File.Close()

fout.cd()
#dstreeclone.Write()

#H2_nSCBayes_rprompt60Bayes_0.Write("H2_nSCBayes_rprompt60Bayes_0")
#H2_qpe_fprompt_0.Write("H2_qpe_fprompt_0")
#
#H2_nSCBayes_rprompt60Bayes_1.Write("H2_nSCBayes_rprompt60Bayes_1")
#H2_qpe_fprompt_1.Write("H2_qpe_fprompt_1")
#
#H2_nSCBayes_rprompt60Bayes_2.Write("H2_nSCBayes_rprompt60Bayes_2")
#H2_qpe_fprompt_2.Write("H2_qpe_fprompt_2")
#
#H2_nSCBayes_rprompt60Bayes_3.Write("H2_nSCBayes_rprompt60Bayes_3")
#H2_qpe_fprompt_3.Write("H2_qpe_fprompt_3")
#
#H2_nSCBayes_rprompt60Bayes_4.Write("H2_nSCBayes_rprompt60Bayes_4")
#H2_qpe_fprompt_4.Write("H2_qpe_fprompt_4")
#
#H2_nSCBayes_rprompt60Bayes_5.Write("H2_nSCBayes_rprompt60Bayes_5")
#H2_qpe_fprompt_5.Write("H2_qpe_fprompt_5")
#
#H2_nSCBayes_rprompt60Bayes_6.Write("H2_nSCBayes_rprompt60Bayes_6")
#H2_qpe_fprompt_6.Write("H2_qpe_fprompt_6")

H2_nSCBayes_rprompt60Bayes_7.Write("H2_nSCBayes_rprompt60Bayes_7")
H2_qpe_fprompt_7.Write("H2_qpe_fprompt_7")

## Level 8
H2_nSCBayes_rprompt60Bayes_7_fmaxpe.Write()
H2_qpe_fprompt_7_fmaxpe.Write()
## Level 9
H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck.Write()
H2_qpe_fprompt_7_fmaxpe_neck.Write()
## Level 10
H2_nSCBayes_rprompt60Bayes_7_fmaxpe_neck_mbR.Write()
H2_qpe_fprompt_7_fmaxpe_neck_mbR.Write()

#H2_rhoZ_0.Write()
#H2_nhitQPE_0.Write()
#H2_nhitNSCB_0.Write()
#H1_fmaxpe_0.Write()
#H1_scintlike_0.Write()
#H1_llneutron_0.Write()
#
#H2_rhoZ_1.Write()
#H2_nhitQPE_1.Write()
#H2_nhitNSCB_1.Write()
#H1_fmaxpe_1.Write()
#H1_scintlike_1.Write()
#H1_llneutron_1.Write()
#
#H2_rhoZ_2.Write()
#H2_nhitQPE_2.Write()
#H2_nhitNSCB_2.Write()
#H1_fmaxpe_2.Write()
#H1_scintlike_2.Write()
#H1_llneutron_2.Write()
#
#H2_rhoZ_3.Write()
#H2_nhitQPE_3.Write()
#H2_nhitNSCB_3.Write()
#H1_fmaxpe_3.Write()
#H1_scintlike_3.Write()
#H1_llneutron_3.Write()
#
#H2_rhoZ_4.Write()
#H2_nhitQPE_4.Write()
#H2_nhitNSCB_4.Write()
#H1_fmaxpe_4.Write()
#H1_scintlike_4.Write()
#H1_llneutron_4.Write()
#
#H2_rhoZ_5.Write()
#H2_nhitQPE_5.Write()
#H2_nhitNSCB_5.Write()
#H1_fmaxpe_5.Write()
#H1_scintlike_5.Write()
#H1_llneutron_5.Write()
#
#H2_rhoZ_6.Write()
#H2_nhitQPE_6.Write()
#H2_nhitNSCB_6.Write()
#H1_fmaxpe_6.Write()
#H1_scintlike_6.Write()
#H1_llneutron_6.Write()
#
#H2_rhoZ_7.Write()
#H2_nhitQPE_7.Write()
#H2_nhitNSCB_7.Write()
#H1_fmaxpe_7.Write()
#H1_scintlike_7.Write()
#H1_llneutron_7.Write()

H2_rhoZ_7_fmaxpe.Write()
H2_nhitQPE_7_fmaxpe.Write()
H2_nhitNSCB_7_fmaxpe.Write()
H1_fmaxpe_7_fmaxpe.Write()
H1_scintlike_7_fmaxpe.Write()
H1_llneutron_7_fmaxpe.Write()

H2_rhoZ_7_fmaxpe_neck.Write()
H2_nhitQPE_7_fmaxpe_neck.Write()
H2_nhitNSCB_7_fmaxpe_neck.Write()
H1_fmaxpe_7_fmaxpe_neck.Write()
H1_scintlike_7_fmaxpe_neck.Write()
H1_llneutron_7_fmaxpe_neck.Write()

H2_rhoZ_7_fmaxpe_neck_mbR.Write()
H2_nhitQPE_7_fmaxpe_neck_mbR.Write()
H2_nhitNSCB_7_fmaxpe_neck_mbR.Write()
H1_fmaxpe_7_fmaxpe_neck_mbR.Write()
H1_scintlike_7_fmaxpe_neck_mbR.Write()
H1_llneutron_7_fmaxpe_neck_mbR.Write()

H_countTotal.Write()
H_countROI05.Write()
H_countROI30.Write()
H_countROI55.Write()
H_countSideBand.Write()
ntupleTF2.Write()
fout.Close()
