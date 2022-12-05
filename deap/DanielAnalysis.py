import sys
import glob
from rat import *
from ROOT import *
import numpy as np
import struct
#from run_numbers import runs
from datetime import datetime, time as datetime_time, timedelta

#runs_name = ["U232Source_SLB005_2020_CalC_L0", "U232Source_SLB005_2020_CalC_RT465_L0", "U232Source_SLB005_2020_CalC_RT480_L0", "U232Source_SLB005_2020_CalF_L0", "U232Source_SLB005_2020_CalFpos6_RT465_L0", "U232Source_SLB005_2020_CalFpos6_RT480_L0", "U232Source_SLB005_2020_CalFpos7_RT465_L0", "U232Source_SLB005_2020_CalFpos7_RT480_L0", "U232Source_SLB005_2021_CalFpos6_RT480_L0"]
#
#RunList = runs
#
#Run            = sys.argv[1]
file0 = str(sys.argv[1])

###############################################################################################################################
H2_nSCBayes_rprompt60Bayes_0          = TH2D( "H2_nSCBayes_rprompt60Bayes_0"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_0                      = TH2D( "H2_qpe_fprompt_0"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_1          = TH2D( "H2_nSCBayes_rprompt60Bayes_1"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_1                      = TH2D( "H2_qpe_fprompt_1"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_2          = TH2D( "H2_nSCBayes_rprompt60Bayes_2"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_2                      = TH2D( "H2_qpe_fprompt_2"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_3          = TH2D( "H2_nSCBayes_rprompt60Bayes_3"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_3                      = TH2D( "H2_qpe_fprompt_3"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_4          = TH2D( "H2_nSCBayes_rprompt60Bayes_4"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_4                      = TH2D( "H2_qpe_fprompt_4"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_5          = TH2D( "H2_nSCBayes_rprompt60Bayes_5"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_5                      = TH2D( "H2_qpe_fprompt_5"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_6          = TH2D( "H2_nSCBayes_rprompt60Bayes_6"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_6                      = TH2D( "H2_qpe_fprompt_6"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )

H2_nSCBayes_rprompt60Bayes_7          = TH2D( "H2_nSCBayes_rprompt60Bayes_7"         , "; nSCBayes [1/bin]; rprompt60Bayes [0.01/bin]"        , 2000, 0, 2000, 100, 0.0, 1.0 )
H2_qpe_fprompt_7                      = TH2D( "H2_qpe_fprompt_7"                     , "; qpe [1/bin]; fprompt [0.01/bin]"                    , 2000, 0, 2000, 100, 0.0, 1.0 )
###########################################################################################################################
#Files = AllFiles[subrun];
File = TFile(file0)
print ("\n********************************** File " , File , " loaded... **********************************\n")
data = File.Get("data_satCorr");#TTree
nentries = data.GetEntries();           update = int(nentries/10);
fout = TFile("Daniel_"+file0, "RECREATE")
###########################################################################################################################
#foutfile="/project/6004969/dpapi/cherenkov/skim/{}/{}".format(runs_name[int(Run)-1], value.split("/")[-1])
#fout = TFile(foutfile,"RECREATE")
dstreeclone = data.CloneTree(0);
dstreeclone.SetDirectory(fout);
###########################################################################################################################
for event in range (nentries):
    if (event+1)%update == 0:
        print (event+1), "Analyzed..."
    data.GetEntry(event)
    flag=0
    H2_nSCBayes_rprompt60Bayes_0  .Fill(data.nSCBayes, data.rprompt60Bayes)
    H2_qpe_fprompt_0              .Fill(data.qPE     , data.fprompt       )

    if data.dtmTrigSrc&0x82 == 0:
        H2_nSCBayes_rprompt60Bayes_1  .Fill(data.nSCBayes, data.rprompt60Bayes)
        H2_qpe_fprompt_1              .Fill(data.qPE     , data.fprompt       )
        if data.calcut&0x31f8 == 0:
            H2_nSCBayes_rprompt60Bayes_2  .Fill(data.nSCBayes, data.rprompt60Bayes)
            H2_qpe_fprompt_2              .Fill(data.qPE     , data.fprompt       )
            if data.deltat > 20000:
                H2_nSCBayes_rprompt60Bayes_3  .Fill(data.nSCBayes, data.rprompt60Bayes)
                H2_qpe_fprompt_3              .Fill(data.qPE     , data.fprompt       )
                if data.numEarlyPulses <= 3:
                    H2_nSCBayes_rprompt60Bayes_4  .Fill(data.nSCBayes, data.rprompt60Bayes)
                    H2_qpe_fprompt_4              .Fill(data.qPE     , data.fprompt       )
                    if data.subeventN==1:
                        H2_nSCBayes_rprompt60Bayes_5  .Fill(data.nSCBayes, data.rprompt60Bayes)
                        H2_qpe_fprompt_5              .Fill(data.qPE     , data.fprompt       )
                        if data.eventTime > 2250  and data.eventTime < 2700:
                            H2_nSCBayes_rprompt60Bayes_6  .Fill(data.nSCBayes, data.rprompt60Bayes)
                            H2_qpe_fprompt_6              .Fill(data.qPE     , data.fprompt       )
                            if data.qPE>60:
                                H2_nSCBayes_rprompt60Bayes_7  .Fill(data.nSCBayes, data.rprompt60Bayes)
                                H2_qpe_fprompt_7              .Fill(data.qPE     , data.fprompt       )
                                flag=1
    if flag:
        dstreeclone.Fill()

#################################################################################################################################

#dstreeclone.Write()
#fout.Close()
#File.Close()

fout.cd()
H2_nSCBayes_rprompt60Bayes_0    .Write( "H2_nSCBayes_rprompt60Bayes_0")
H2_qpe_fprompt_0                .Write( "H2_qpe_fprompt_0"            )

H2_nSCBayes_rprompt60Bayes_1    .Write( "H2_nSCBayes_rprompt60Bayes_1")
H2_qpe_fprompt_1                .Write( "H2_qpe_fprompt_1"            )

H2_nSCBayes_rprompt60Bayes_2    .Write( "H2_nSCBayes_rprompt60Bayes_2")
H2_qpe_fprompt_2                .Write( "H2_qpe_fprompt_2"            )

H2_nSCBayes_rprompt60Bayes_3    .Write( "H2_nSCBayes_rprompt60Bayes_3")
H2_qpe_fprompt_3                .Write( "H2_qpe_fprompt_3"            )

H2_nSCBayes_rprompt60Bayes_4    .Write( "H2_nSCBayes_rprompt60Bayes_4")
H2_qpe_fprompt_4                .Write( "H2_qpe_fprompt_4"            )

H2_nSCBayes_rprompt60Bayes_5    .Write( "H2_nSCBayes_rprompt60Bayes_5")
H2_qpe_fprompt_5                .Write( "H2_qpe_fprompt_5"            )

H2_nSCBayes_rprompt60Bayes_6    .Write( "H2_nSCBayes_rprompt60Bayes_6")
H2_qpe_fprompt_6                .Write( "H2_qpe_fprompt_6"            )

H2_nSCBayes_rprompt60Bayes_7    .Write( "H2_nSCBayes_rprompt60Bayes_7")
H2_qpe_fprompt_7                .Write( "H2_qpe_fprompt_7"            )

fout.Close();
