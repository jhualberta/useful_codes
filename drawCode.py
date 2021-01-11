i = 0.1
for k in range(100):
    ss = str(i)
    ss1 = ss.replace('.','p')
    s1 = "TFile *f" +str(k+1)+"= new TFile(\"MC_Water_"+ss1+"_beta_IsoCenter.root\");"
    s2 = "TH1F *hNCherPhoton"+ss1+"MeV = new TH1F(\"hNCherPhoton"+ss1+"MeV\",\"\",240,0,6000);"
    s3 = "TTree* t"+str(k+1)+" = (TTree*)f"+str(k+1)+"->Get(\"T\");"
    s4 = "t"+str(k+1)+"->Project(\"hNCherPhoton"+ss1+"MeV\",\"ds.mc.nCherPhotons\");"
    i = i + 0.1
    print s1
    print s2
    print s3
    print s4

i = 0.1
for k in range(100):
    ss = str(i)
    ss1 = ss.replace('.','p')
    s2 = "hNCherPhoton"+ss1+"MeV->Write();"
    print s2
    i = i + 0.1
