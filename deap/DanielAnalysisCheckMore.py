import sys
import glob
from rat import *
from ROOT import *
import numpy as np
from numpy import sqrt
import struct
#from run_numbers import runs
from datetime import datetime, time as datetime_time, timedelta
ffroi = TFile("saveSideBandROI.root")
####roi_802_days_11March2020_nsc_rp60.root")
roicut = ffroi.Get("roi")
roicut.SetName("roicut")
sidebandcut = ffroi.Get("cutSideBand")
sidebandcut.SetName("cutSideBand")

#runs_name = ["U232Source_SLB005_2020_CalC_L0", "U232Source_SLB005_2020_CalC_RT465_L0", "U232Source_SLB005_2020_CalC_RT480_L0", "U232Source_SLB005_2020_CalF_L0", "U232Source_SLB005_2020_CalFpos6_RT465_L0", "U232Source_SLB005_2020_CalFpos6_RT480_L0", "U232Source_SLB005_2020_CalFpos7_RT465_L0", "U232Source_SLB005_2020_CalFpos7_RT480_L0", "U232Source_SLB005_2021_CalFpos6_RT480_L0"]
#
#RunList = runs
#
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
### 7-level cuts + STR levels, Level 8
H2_nSCBayes_rprompt60Bayes_7_neck = TH2D( "H2_nSCBayes_rprompt60Bayes_7_neck", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_neck = TH2D( "H2_qpe_fprompt_7_neck", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_neck = TH2D("H2_rhoZ_7_neck", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_neck = TH2D("H2_nhitQPE_7_neck", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_neck = TH2D("H2_nhitNSCB_7_neck", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_neck = TH1D("H1_fmaxpe_7_neck","fmaxpe",100,0,1)
H1_scintlike_7_neck = TH1D("H1_scintlike_7_neck","scintlike",1000,0,1000)
H1_llneutron_7_neck = TH1D("H1_llneutron_7_neck","llneutron",200,-100,100)
## Level 9
H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe = TH2D( "H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_neck_fmaxpe = TH2D( "H2_qpe_fprompt_7_neck_fmaxpe", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_neck_fmaxpe = TH2D("H2_rhoZ_7_neck_fmaxpe", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_neck_fmaxpe = TH2D("H2_nhitQPE_7_neck_fmaxpe", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_neck_fmaxpe = TH2D("H2_nhitNSCB_7_neck_fmaxpe", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_neck_fmaxpe = TH1D("H1_fmaxpe_7_neck_fmaxpe","fmaxpe",100,0,1)
H1_scintlike_7_neck_fmaxpe = TH1D("H1_scintlike_7_neck_fmaxpe","scintlike",1000,0,1000)
H1_llneutron_7_neck_fmaxpe = TH1D("H1_llneutron_7_neck_fmaxpe","llneutron",200,-100,100)
## Level 10
H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe_mbR = TH2D( "H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe_mbR", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_neck_fmaxpe_mbR = TH2D( "H2_qpe_fprompt_7_neck_fmaxpe_mbR", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_neck_fmaxpe_mbR = TH2D("H2_rhoZ_7_neck_fmaxpe_mbR", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_neck_fmaxpe_mbR = TH2D("H2_nhitQPE_7_neck_fmaxpe_mbR", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_neck_fmaxpe_mbR = TH2D("H2_nhitNSCB_7_neck_fmaxpe_mbR", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_neck_fmaxpe_mbR = TH1D("H1_fmaxpe_7_neck_fmaxpe_mbR","fmaxpe", 100,0,1)
H1_scintlike_7_neck_fmaxpe_mbR = TH1D("H1_scintlike_7_neck_fmaxpe_mbR","scintlike",1000,0,1000)
H1_llneutron_7_neck_fmaxpe_mbR = TH1D("H1_llneutron_7_neck_fmaxpe_mbR","llneutron",200,-100,100)

### 7-level cuts + STR + LAr level cuts, Level 11
H2_nSCBayes_rprompt60Bayes_7_str_qRatio = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_qRatio", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_qRatio = TH2D( "H2_qpe_fprompt_7_str_qRatio", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_qRatio = TH2D("H2_rhoZ_7_str_qRatio", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_qRatio = TH2D("H2_nhitQPE_7_str_qRatio", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_qRatio = TH2D("H2_nhitNSCB_7_str_qRatio", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_qRatio = TH1D("H1_fmaxpe_7_str_qRatio","fmaxpe",100,0,1)
H1_scintlike_7_str_qRatio = TH1D("H1_scintlike_7_str_qRatio","scintlike",1000,0,1000)
H1_llneutron_7_str_qRatio = TH1D("H1_llneutron_7_str_qRatio","llneutron",200,-100,100)
### Level 12
H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_qRatio_mbZ = TH2D( "H2_qpe_fprompt_7_str_qRatio_mbZ", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_qRatio_mbZ = TH2D("H2_rhoZ_7_str_qRatio_mbZ", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_qRatio_mbZ = TH2D("H2_nhitQPE_7_str_qRatio_mbZ", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_qRatio_mbZ = TH2D("H2_nhitNSCB_7_str_qRatio_mbZ", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_qRatio_mbZ = TH1D("H1_fmaxpe_7_str_qRatio_mbZ","fmaxpe",100,0,1)
H1_scintlike_7_str_qRatio_mbZ = TH1D("H1_scintlike_7_str_qRatio_mbZ","scintlike",1000,0,1000)
H1_llneutron_7_str_qRatio_mbZ = TH1D("H1_llneutron_7_str_qRatio_mbZ","llneutron",200,-100,100)
### Level 13
H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ_pulseG = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ_pulseG", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_qRatio_mbZ_pulseG = TH2D( "H2_qpe_fprompt_7_str_qRatio_mbZ_pulseG", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_qRatio_mbZ_pulseG = TH2D("H2_rhoZ_7_str_qRatio_mbZ_pulseG", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_qRatio_mbZ_pulseG = TH2D("H2_nhitQPE_7_str_qRatio_mbZ_pulseG", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_qRatio_mbZ_pulseG = TH2D("H2_nhitNSCB_7_str_qRatio_mbZ_pulseG", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_qRatio_mbZ_pulseG = TH1D("H1_fmaxpe_7_str_qRatio_mbZ_pulseG","fmaxpe",100,0,1)
H1_scintlike_7_str_qRatio_mbZ_pulseG = TH1D("H1_scintlike_7_str_qRatio_mbZ_pulseG","scintlike",1000,0,1000)
H1_llneutron_7_str_qRatio_mbZ_pulseG = TH1D("H1_llneutron_7_str_qRatio_mbZ_pulseG","llneutron", 200,-100,100)

#### TF2 - MB consistency cuts, Level 14: z consistency
H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_LAr_zCon = TH2D( "H2_qpe_fprompt_7_str_LAr_zCon", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_LAr_zCon = TH2D("H2_rhoZ_7_str_LAr_zCon", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_LAr_zCon = TH2D("H2_nhitQPE_7_str_LAr_zCon", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_LAr_zCon = TH2D("H2_nhitNSCB_7_str_LAr_zCon", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_LAr_zCon = TH1D("H1_fmaxpe_7_str_LAr_zCon","fmaxpe",100,0,1)
H1_scintlike_7_str_LAr_zCon = TH1D("H1_scintlike_7_str_LAr_zCon","scintlike",1000,0,1000)
H1_llneutron_7_str_LAr_zCon = TH1D("H1_llneutron_7_str_LAr_zCon","llneutron", 200,-100,100)

### TF2 - MB consistency cuts, Level 15
H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_LAr_zrCon = TH2D( "H2_qpe_fprompt_7_str_LAr_zrCon", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_LAr_zrCon = TH2D("H2_rhoZ_7_str_LAr_zrCon", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_LAr_zrCon = TH2D("H2_nhitQPE_7_str_LAr_zrCon", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_LAr_zrCon = TH2D("H2_nhitNSCB_7_str_LAr_zrCon", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_LAr_zrCon = TH1D("H1_fmaxpe_7_str_LAr_zrCon","fmaxpe",100,0,1)
H1_scintlike_7_str_LAr_zrCon = TH1D("H1_scintlike_7_str_LAr_zrCon","scintlike",1000,0,1000)
H1_llneutron_7_str_LAr_zrCon = TH1D("H1_llneutron_7_str_LAr_zrCon","llneutron", 200,-100,100)

#### TF2 - MB consistency cuts, Level 14 a: z consistency
H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon_tf2Cerenkov = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon_tf2Cerenkov", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_LAr_zCon_tf2Cerenkov = TH2D( "H2_qpe_fprompt_7_str_LAr_zCon_tf2Cerenkov", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_LAr_zCon_tf2Cerenkov = TH2D("H2_rhoZ_7_str_LAr_zCon_tf2Cerenkov", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_LAr_zCon_tf2Cerenkov = TH2D("H2_nhitQPE_7_str_LAr_zCon_tf2Cerenkov", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_LAr_zCon_tf2Cerenkov = TH2D("H2_nhitNSCB_7_str_LAr_zCon_tf2Cerenkov", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_LAr_zCon_tf2Cerenkov = TH1D("H1_fmaxpe_7_str_LAr_zCon_tf2Cerenkov","fmaxpe",100,0,1)
H1_scintlike_7_str_LAr_zCon_tf2Cerenkov = TH1D("H1_scintlike_7_str_LAr_zCon_tf2Cerenkov","scintlike",1000,0,1000)
H1_llneutron_7_str_LAr_zCon_tf2Cerenkov = TH1D("H1_llneutron_7_str_LAr_zCon_tf2Cerenkov","llneutron", 200,-100,100)

### TF2 - MB consistency cuts, Level 15 a
H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon_tf2Cerenkov = TH2D( "H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon_tf2Cerenkov", "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7_str_LAr_zrCon_tf2Cerenkov = TH2D( "H2_qpe_fprompt_7_str_LAr_zrCon_tf2Cerenkov", "; qpe [1/bin]; fprompt [0.01/bin]", 2000, 0, 2000, 100, 0.0, 1.0 )
H2_rhoZ_7_str_LAr_zrCon_tf2Cerenkov = TH2D("H2_rhoZ_7_str_LAr_zrCon_tf2Cerenkov", "rho vs z", 850, 0, 850, 1700, -850, 850)
H2_nhitQPE_7_str_LAr_zrCon_tf2Cerenkov = TH2D("H2_nhitQPE_7_str_LAr_zrCon_tf2Cerenkov", "qpe vs nhit/qpe", 2000, 0, 2000, 100, 0, 1)
H2_nhitNSCB_7_str_LAr_zrCon_tf2Cerenkov = TH2D("H2_nhitNSCB_7_str_LAr_zrCon_tf2Cerenkov", "qpe vs nhit/nSCB", 2000, 0, 2000, 100, 0, 1)
H1_fmaxpe_7_str_LAr_zrCon_tf2Cerenkov = TH1D("H1_fmaxpe_7_str_LAr_zrCon_tf2Cerenkov","fmaxpe",100,0,1)
H1_scintlike_7_str_LAr_zrCon_tf2Cerenkov = TH1D("H1_scintlike_7_str_LAr_zrCon_tf2Cerenkov","scintlike",1000,0,1000)
H1_llneutron_7_str_LAr_zrCon_tf2Cerenkov = TH1D("H1_llneutron_7_str_LAr_zrCon_tf2Cerenkov","llneutron", 200,-100,100)

#Files = AllFiles[subrun];
File = TFile(file0)
print ("\n********************************** File " , File , " loaded... **********************************\n")
data = File.Get("data_satCorr");#TTree
nentries = data.GetEntries();
update = int(nentries/10);
fout = TFile("Daniel_"+fileName, "RECREATE")
###########################################################################################################################
#foutfile="/project/6004969/dpapi/cherenkov/skim/{}/{}".format(runs_name[int(Run)-1], value.split("/")[-1])
#fout = TFile(foutfile,"RECREATE")
dstreeclone = data.CloneTree(0);
dstreeclone.SetDirectory(fout);

level = 20
countTotal = [0 for i in range(level)]
countROI = [0 for i in range(level)]
countSideband = [0 for i in range(level)]

H_countTotal = TH1F("H_countTotal", "total", level, 0, level)
H_countROI = TH1F("H_countROI", "roi cuts", level, 0, level)
H_countSideBand = TH1F("H_countSideBand", "sideband cuts", level, 0, level)


###########################################################################################################################

for event in range(nentries):
    if (event+1)%update == 0:
        print (event+1), "Analyzed..."
    data.GetEntry(event)
    flag=0
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
    if roicut.IsInside(nscb, rprompt60Bayes): countROI[0] += 1
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
        if roicut.IsInside(nscb, rprompt60Bayes): countROI[1] += 1
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
           if roicut.IsInside(nscb, rprompt60Bayes): countROI[2] += 1
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
              if roicut.IsInside(nscb, rprompt60Bayes): countROI[3] += 1
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
                 if roicut.IsInside(nscb, rprompt60Bayes): countROI[4] += 1
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
                    if roicut.IsInside(nscb, rprompt60Bayes): countROI[5] += 1
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
                        if roicut.IsInside(nscb, rprompt60Bayes): countROI[6] += 1
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
                          if roicut.IsInside(nscb, rprompt60Bayes): countROI[7] += 1
                          if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[7] += 1
                          flag=1
                          ### Note: for some events, qpe could be 0
                          cft2r = (data.chargetopring + data.chargesecondring)/qpe
                          cfb3r = (data.chargebottomring + data.chargesecondbottomring + data.chargethirdbottomring)/qpe

                          ### 7-level cuts + STR levels
                          if data.neckVetoN  == 0: ## level 8 
                            H2_nSCBayes_rprompt60Bayes_7_neck.Fill(nscb, rprompt60Bayes)
                            H2_qpe_fprompt_7_neck.Fill(qpe, fprompt)
                            H2_rhoZ_7_neck.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                            if qpe>0: H2_nhitQPE_7_neck.Fill(qpe, nhit/qpe)
                            if nscb>0: H2_nhitNSCB_7_neck.Fill(nscb, nhit/nscb)
                            H1_fmaxpe_7_neck.Fill(fmaxpe)
                            H1_scintlike_7_neck.Fill(scintlike)
                            H1_llneutron_7_neck.Fill(llneutronflash)
                            countTotal[8] += 1
                            if roicut.IsInside(nscb, rprompt60Bayes): countROI[8] += 1
                            if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[8] += 1

                            if data.fmaxpe<0.4: ## level 9
                              H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe.Fill(nscb, rprompt60Bayes)
                              H2_qpe_fprompt_7_neck_fmaxpe.Fill(qpe, fprompt)
                              H2_rhoZ_7_neck_fmaxpe.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                              if qpe>0: H2_nhitQPE_7_neck_fmaxpe.Fill(qpe, nhit/qpe)
                              if nscb>0: H2_nhitNSCB_7_neck_fmaxpe.Fill(nscb, nhit/nscb)
                              H1_fmaxpe_7_neck_fmaxpe.Fill(fmaxpe)
                              H1_scintlike_7_neck_fmaxpe.Fill(scintlike)
                              H1_llneutron_7_neck_fmaxpe.Fill(llneutronflash)
                              countTotal[9] += 1
                              if roicut.IsInside(nscb, rprompt60Bayes): countROI[9] += 1
                              if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[9] += 1

                              if mbR<630: ## level 10 
                                H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe_mbR.Fill(nscb, rprompt60Bayes)
                                H2_qpe_fprompt_7_neck_fmaxpe_mbR.Fill(qpe, fprompt)
                                H2_rhoZ_7_neck_fmaxpe_mbR.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                if qpe>0: H2_nhitQPE_7_neck_fmaxpe_mbR.Fill(qpe, nhit/qpe)
                                if nscb>0: H2_nhitNSCB_7_neck_fmaxpe_mbR.Fill(nscb, nhit/nscb)
                                H1_fmaxpe_7_neck_fmaxpe_mbR.Fill(fmaxpe)
                                H1_scintlike_7_neck_fmaxpe_mbR.Fill(scintlike)
                                H1_llneutron_7_neck_fmaxpe_mbR.Fill(llneutronflash)
                                countTotal[10] += 1
                                if roicut.IsInside(nscb, rprompt60Bayes): countROI[10] += 1
                                if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[10] += 1
                        
                                ### LAr cuts
                                if cft2r<0.04 and cfb3r<0.1: ## level 11
                                     H2_nSCBayes_rprompt60Bayes_7_str_qRatio.Fill(nscb, rprompt60Bayes)
                                     H2_qpe_fprompt_7_str_qRatio.Fill(qpe, fprompt)
                                     H2_rhoZ_7_str_qRatio.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                     if qpe>0: H2_nhitQPE_7_str_qRatio.Fill(qpe, nhit/qpe)
                                     if nscb>0: H2_nhitNSCB_7_str_qRatio.Fill(nscb, nhit/nscb)
                                     H1_fmaxpe_7_str_qRatio.Fill(fmaxpe)
                                     H1_scintlike_7_str_qRatio.Fill(scintlike)
                                     H1_llneutron_7_str_qRatio.Fill(llneutronflash)
                                     countTotal[11] += 1
                                     if roicut.IsInside(nscb, rprompt60Bayes): countROI[11] += 1
                                     if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[11] += 1

                                     if mbZ<550:## level 12
                                       H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ.Fill(nscb, rprompt60Bayes)
                                       H2_qpe_fprompt_7_str_qRatio_mbZ.Fill(qpe, fprompt)
                                       H2_rhoZ_7_str_qRatio_mbZ.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                       if qpe>0: H2_nhitQPE_7_str_qRatio_mbZ.Fill(qpe, nhit/qpe)
                                       if nscb>0: H2_nhitNSCB_7_str_qRatio_mbZ.Fill(nscb, nhit/nscb)
                                       H1_fmaxpe_7_str_qRatio_mbZ.Fill(fmaxpe)
                                       H1_scintlike_7_str_qRatio_mbZ.Fill(scintlike)
                                       H1_llneutron_7_str_qRatio_mbZ.Fill(llneutronflash)
                                       countTotal[12] += 1
                                       if roicut.IsInside(nscb, rprompt60Bayes): countROI[12] += 1
                                       if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[12] += 1

                                       if data.pulseindexfirstgar>2: ## level 13
                                          H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ_pulseG.Fill(nscb, rprompt60Bayes)
                                          H2_qpe_fprompt_7_str_qRatio_mbZ_pulseG.Fill(qpe, fprompt)
                                          H2_rhoZ_7_str_qRatio_mbZ_pulseG.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                          if qpe>0: H2_nhitQPE_7_str_qRatio_mbZ_pulseG.Fill(qpe, nhit/qpe)
                                          if nscb>0: H2_nhitNSCB_7_str_qRatio_mbZ_pulseG.Fill(nscb, nhit/nscb)
                                          H1_fmaxpe_7_str_qRatio_mbZ_pulseG.Fill(fmaxpe)
                                          H1_scintlike_7_str_qRatio_mbZ_pulseG.Fill(scintlike)
                                          H1_llneutron_7_str_qRatio_mbZ_pulseG.Fill(llneutronflash)
                                          countTotal[13] += 1
                                          if roicut.IsInside(nscb, rprompt60Bayes): countROI[13] += 1
                                          if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[13] += 1

                                          ### level 14, Position reconstruction consistency contours require agreement between MBLikelihood and TimeFit2
                                          promptPE = nscb*rprompt60Bayes;
                                          mbPos = TVector3(mbX,mbY,mbZ);
                                          tf2X = data.timefit2X; tf2Y = data.timefit2Y; tf2Z = data.timefit2Z;
                                          tf2cerenkovX = data.timefit2cerenkovX; tf2cerenkovY = data.timefit2cerenkovY; tf2cerenkovZ = data.timefit2cerenkovZ; 

                                          tf2Pos = TVector3(tf2X,tf2Y,tf2Z);
                                          tf2CerenPos = TVector3(tf2cerenkovX,tf2cerenkovY,tf2cerenkovZ);

                                          ### Contency based cut, reconstructed Z (90% Ar39 acceptance)
                                          zContF = TFile("tf2mb_nSCBayes_deltaz_contours.root","READ");
                                          zCut = zContF.Get("cont90_cut;1");

                                          rContF = TFile("tf2mb_nSCBayes_dist_after_deltaz90_contours.root","READ");
                                          rCut = rContF.Get("cont85_cut;1");

                                          if (zCut.IsInside(promptPE,(tf2Pos.Z()-mbPos.Z()))):
                                             H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon.Fill(nscb, rprompt60Bayes)
                                             H2_qpe_fprompt_7_str_LAr_zCon.Fill(qpe, fprompt)
                                             H2_rhoZ_7_str_LAr_zCon.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                             if qpe>0: H2_nhitQPE_7_str_LAr_zCon.Fill(qpe, nhit/qpe)
                                             if nscb>0: H2_nhitNSCB_7_str_LAr_zCon.Fill(nscb, nhit/nscb)
                                             H1_fmaxpe_7_str_LAr_zCon.Fill(fmaxpe)
                                             H1_scintlike_7_str_LAr_zCon.Fill(scintlike)
                                             H1_llneutron_7_str_LAr_zCon.Fill(llneutronflash)
                                             countTotal[14] += 1
                                             if roicut.IsInside(nscb, rprompt60Bayes): countROI[14] += 1
                                             if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[14] += 1

                                             ### level 15, Contency based cut, 3D distance between MB and TF2 vertices (85% Ar39 acceptance)
                                             if (rCut.IsInside(promptPE,(tf2Pos-mbPos).Mag())):
                                                  H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon.Fill(nscb, rprompt60Bayes)
                                                  H2_qpe_fprompt_7_str_LAr_zrCon.Fill(qpe, fprompt)
                                                  H2_rhoZ_7_str_LAr_zrCon.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                                  if qpe>0: H2_nhitQPE_7_str_LAr_zrCon.Fill(qpe, nhit/qpe)
                                                  if nscb>0: H2_nhitNSCB_7_str_LAr_zrCon.Fill(nscb, nhit/nscb)
                                                  H1_fmaxpe_7_str_LAr_zrCon.Fill(fmaxpe)
                                                  H1_scintlike_7_str_LAr_zrCon.Fill(scintlike)
                                                  H1_llneutron_7_str_LAr_zrCon.Fill(llneutronflash)
                                                  countTotal[15] += 1
                                                  if roicut.IsInside(nscb, rprompt60Bayes): countROI[15] += 1
                                                  if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[15] += 1

                                          ### using tf2 cerenkov as position consistency, level 16, acutally level 14-a
                                          if (zCut.IsInside(promptPE,(tf2CerenPos.Z()-mbPos.Z()))):
                                             H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon_tf2Cerenkov.Fill(nscb, rprompt60Bayes)
                                             H2_qpe_fprompt_7_str_LAr_zCon_tf2Cerenkov.Fill(qpe, fprompt)
                                             H2_rhoZ_7_str_LAr_zCon_tf2Cerenkov.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                             if qpe>0: H2_nhitQPE_7_str_LAr_zCon_tf2Cerenkov.Fill(qpe, nhit/qpe)
                                             if nscb>0: H2_nhitNSCB_7_str_LAr_zCon_tf2Cerenkov.Fill(nscb, nhit/nscb)
                                             H1_fmaxpe_7_str_LAr_zCon_tf2Cerenkov.Fill(fmaxpe)
                                             H1_scintlike_7_str_LAr_zCon_tf2Cerenkov.Fill(scintlike)
                                             H1_llneutron_7_str_LAr_zCon_tf2Cerenkov.Fill(llneutronflash)
                                             countTotal[16] += 1
                                             if roicut.IsInside(nscb, rprompt60Bayes): countROI[16] += 1
                                             if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[16] += 1

                                             ### level 17 (or level 15-a), Contency based cut, 3D distance between MB and TF2 vertices (85% Ar39 acceptance)
                                             if (rCut.IsInside(promptPE,(tf2CerenPos-mbPos).Mag())):
                                                  H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon_tf2Cerenkov.Fill(nscb, rprompt60Bayes)
                                                  H2_qpe_fprompt_7_str_LAr_zrCon_tf2Cerenkov.Fill(qpe, fprompt)
                                                  H2_rhoZ_7_str_LAr_zrCon_tf2Cerenkov.Fill(sqrt(mbX*mbX+mbY*mbY),mbZ)
                                                  if qpe>0: H2_nhitQPE_7_str_LAr_zrCon_tf2Cerenkov.Fill(qpe, nhit/qpe)
                                                  if nscb>0: H2_nhitNSCB_7_str_LAr_zrCon_tf2Cerenkov.Fill(nscb, nhit/nscb)
                                                  H1_fmaxpe_7_str_LAr_zrCon_tf2Cerenkov.Fill(fmaxpe)
                                                  H1_scintlike_7_str_LAr_zrCon_tf2Cerenkov.Fill(scintlike)
                                                  H1_llneutron_7_str_LAr_zrCon_tf2Cerenkov.Fill(llneutronflash)
                                                  countTotal[17] += 1
                                                  if roicut.IsInside(nscb, rprompt60Bayes): countROI[17] += 1
                                                  if sidebandcut.IsInside(nscb, rprompt60Bayes): countSideband[17] += 1


## fill the counts in ROI for each level cuts
for i in range(0,Level):
   H_countTotal.Fill(i+1,countTotal[i])
   H_countROI.Fill(i+1, countROI[i])
   H_countSideBand.Fill(i+1, countSideband[i])

#    if flag:
#        dstreeclone.Fill()
#################################################################################################################################

#dstreeclone.Write()
#fout.Close()
#File.Close()

fout.cd()
#dstreeclone.Write()

H2_nSCBayes_rprompt60Bayes_0.Write("H2_nSCBayes_rprompt60Bayes_0")
H2_qpe_fprompt_0.Write("H2_qpe_fprompt_0")

H2_nSCBayes_rprompt60Bayes_1.Write( "H2_nSCBayes_rprompt60Bayes_1")
H2_qpe_fprompt_1.Write("H2_qpe_fprompt_1")

H2_nSCBayes_rprompt60Bayes_2.Write( "H2_nSCBayes_rprompt60Bayes_2")
H2_qpe_fprompt_2.Write("H2_qpe_fprompt_2")

H2_nSCBayes_rprompt60Bayes_3.Write( "H2_nSCBayes_rprompt60Bayes_3")
H2_qpe_fprompt_3.Write("H2_qpe_fprompt_3")

H2_nSCBayes_rprompt60Bayes_4.Write( "H2_nSCBayes_rprompt60Bayes_4")
H2_qpe_fprompt_4.Write("H2_qpe_fprompt_4")

H2_nSCBayes_rprompt60Bayes_5.Write( "H2_nSCBayes_rprompt60Bayes_5")
H2_qpe_fprompt_5.Write("H2_qpe_fprompt_5")

H2_nSCBayes_rprompt60Bayes_6.Write( "H2_nSCBayes_rprompt60Bayes_6")
H2_qpe_fprompt_6.Write("H2_qpe_fprompt_6")

H2_nSCBayes_rprompt60Bayes_7.Write("H2_nSCBayes_rprompt60Bayes_7")
H2_qpe_fprompt_7.Write("H2_qpe_fprompt_7")

H2_nSCBayes_rprompt60Bayes_7_neck.Write("H2_nSCBayes_rprompt60Bayes_7_neck")
H2_qpe_fprompt_7_neck.Write("H2_qpe_fprompt_7_neck")

H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe.Write("H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe")
H2_qpe_fprompt_7_neck_fmaxpe.Write("H2_qpe_fprompt_7_neck_fmaxpe")

H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe_mbR.Write("H2_nSCBayes_rprompt60Bayes_7_neck_fmaxpe_mbR")
H2_qpe_fprompt_7_neck_fmaxpe_mbR.Write("H2_qpe_fprompt_7_neck_fmaxpe_mbR")

H2_nSCBayes_rprompt60Bayes_7_str_qRatio.Write("H2_nSCBayes_rprompt60Bayes_7_str_qRatio")
H2_qpe_fprompt_7_str_qRatio.Write("H2_qpe_fprompt_7_str_qRatio")

H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ.Write("H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ")
H2_qpe_fprompt_7_str_qRatio_mbZ.Write("H2_qpe_fprompt_7_str_qRatio_mbZ")

H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ_pulseG.Write("H2_nSCBayes_rprompt60Bayes_7_str_qRatio_mbZ_pulseG")
H2_qpe_fprompt_7_str_qRatio_mbZ_pulseG.Write("H2_qpe_fprompt_7_str_qRatio_mbZ_pulseG")

H2_rhoZ_0.Write()
H2_nhitQPE_0.Write()
H2_nhitNSCB_0.Write()
H1_fmaxpe_0.Write()
H1_scintlike_0.Write()
H1_llneutron_0.Write()

H2_rhoZ_1.Write()
H2_nhitQPE_1.Write()
H2_nhitNSCB_1.Write()
H1_fmaxpe_1.Write()
H1_scintlike_1.Write()
H1_llneutron_1.Write()

H2_rhoZ_2.Write()
H2_nhitQPE_2.Write()
H2_nhitNSCB_2.Write()
H1_fmaxpe_2.Write()
H1_scintlike_2.Write()
H1_llneutron_2.Write()

H2_rhoZ_3.Write()
H2_nhitQPE_3.Write()
H2_nhitNSCB_3.Write()
H1_fmaxpe_3.Write()
H1_scintlike_3.Write()
H1_llneutron_3.Write()

H2_rhoZ_4.Write()
H2_nhitQPE_4.Write()
H2_nhitNSCB_4.Write()
H1_fmaxpe_4.Write()
H1_scintlike_4.Write()
H1_llneutron_4.Write()

H2_rhoZ_5.Write()
H2_nhitQPE_5.Write()
H2_nhitNSCB_5.Write()
H1_fmaxpe_5.Write()
H1_scintlike_5.Write()
H1_llneutron_5.Write()

H2_rhoZ_6.Write()
H2_nhitQPE_6.Write()
H2_nhitNSCB_6.Write()
H1_fmaxpe_6.Write()
H1_scintlike_6.Write()
H1_llneutron_6.Write()

H2_rhoZ_7.Write()
H2_nhitQPE_7.Write()
H2_nhitNSCB_7.Write()
H1_fmaxpe_7.Write()
H1_scintlike_7.Write()
H1_llneutron_7.Write()

H2_rhoZ_7_neck.Write()
H2_nhitQPE_7_neck.Write()
H2_nhitNSCB_7_neck.Write()
H1_fmaxpe_7_neck.Write()
H1_scintlike_7_neck.Write()
H1_llneutron_7_neck.Write()

H2_rhoZ_7_neck_fmaxpe.Write()
H2_nhitQPE_7_neck_fmaxpe.Write()
H2_nhitNSCB_7_neck_fmaxpe.Write()
H1_fmaxpe_7_neck_fmaxpe.Write()
H1_scintlike_7_neck_fmaxpe.Write()
H1_llneutron_7_neck_fmaxpe.Write()

H2_rhoZ_7_neck_fmaxpe_mbR.Write()
H2_nhitQPE_7_neck_fmaxpe_mbR.Write()
H2_nhitNSCB_7_neck_fmaxpe_mbR.Write()
H1_fmaxpe_7_neck_fmaxpe_mbR.Write()
H1_scintlike_7_neck_fmaxpe_mbR.Write()
H1_llneutron_7_neck_fmaxpe_mbR.Write()

H2_rhoZ_7_str_qRatio.Write()
H2_nhitQPE_7_str_qRatio.Write()
H2_nhitNSCB_7_str_qRatio.Write()
H1_fmaxpe_7_str_qRatio.Write()
H1_scintlike_7_str_qRatio.Write()
H1_llneutron_7_str_qRatio.Write()

H2_rhoZ_7_str_qRatio_mbZ.Write()
H2_nhitQPE_7_str_qRatio_mbZ.Write()
H2_nhitNSCB_7_str_qRatio_mbZ.Write()
H1_fmaxpe_7_str_qRatio_mbZ.Write()
H1_scintlike_7_str_qRatio_mbZ.Write()
H1_llneutron_7_str_qRatio_mbZ.Write()

H2_rhoZ_7_str_qRatio_mbZ_pulseG.Write()
H2_nhitQPE_7_str_qRatio_mbZ_pulseG.Write()
H2_nhitNSCB_7_str_qRatio_mbZ_pulseG.Write()
H1_fmaxpe_7_str_qRatio_mbZ_pulseG.Write()
H1_scintlike_7_str_qRatio_mbZ_pulseG.Write()
H1_llneutron_7_str_qRatio_mbZ_pulseG.Write()

H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon.Write()
H2_qpe_fprompt_7_str_LAr_zCon.Write()
H2_rhoZ_7_str_LAr_zCon.Write()
H2_nhitQPE_7_str_LAr_zCon.Write()
H2_nhitNSCB_7_str_LAr_zCon.Write()
H1_fmaxpe_7_str_LAr_zCon.Write()
H1_scintlike_7_str_LAr_zCon.Write()
H1_llneutron_7_str_LAr_zCon.Write()

H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon.Write()
H2_qpe_fprompt_7_str_LAr_zrCon.Write()
H2_rhoZ_7_str_LAr_zrCon.Write()
H2_nhitQPE_7_str_LAr_zrCon.Write()
H2_nhitNSCB_7_str_LAr_zrCon.Write()
H1_fmaxpe_7_str_LAr_zrCon.Write()
H1_scintlike_7_str_LAr_zrCon.Write()
H1_llneutron_7_str_LAr_zrCon.Write()

H2_nSCBayes_rprompt60Bayes_7_str_LAr_zCon_tf2Cerenkov.Write()
H2_qpe_fprompt_7_str_LAr_zCon_tf2Cerenkov.Write()
H2_rhoZ_7_str_LAr_zCon_tf2Cerenkov.Write()
H2_nhitQPE_7_str_LAr_zCon_tf2Cerenkov.Write()
H2_nhitNSCB_7_str_LAr_zCon_tf2Cerenkov.Write()
H1_fmaxpe_7_str_LAr_zCon_tf2Cerenkov.Write()
H1_scintlike_7_str_LAr_zCon_tf2Cerenkov.Write()
H1_llneutron_7_str_LAr_zCon_tf2Cerenkov.Write()

H2_nSCBayes_rprompt60Bayes_7_str_LAr_zrCon_tf2Cerenkov.Write()
H2_qpe_fprompt_7_str_LAr_zrCon_tf2Cerenkov.Write()
H2_rhoZ_7_str_LAr_zrCon_tf2Cerenkov.Write()
H2_nhitQPE_7_str_LAr_zrCon_tf2Cerenkov.Write()
H2_nhitNSCB_7_str_LAr_zrCon_tf2Cerenkov.Write()
H1_fmaxpe_7_str_LAr_zrCon_tf2Cerenkov.Write()
H1_scintlike_7_str_LAr_zrCon_tf2Cerenkov.Write()
H1_llneutron_7_str_LAr_zrCon_tf2Cerenkov.Write()

H_countTotal.Write()
H_countROI.Write()
H_countSideBand.Write()
fout.Close()
