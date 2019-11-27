{
    gStyle->SetOptStat(0);
    TH1F *hTotal = new TH1F("hTotal",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBackTotal = new TH1F("hBackTotal",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hSigTotal = new TH1F("hSigTotal",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBackExt = new TH1F("hBackExt",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBackOther = new TH1F("hBackOther",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAc228 = new TH1F("hAc228",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAlphan_Lab_13c = new TH1F("hAlphan_Lab_13c",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAlphan_Lab_Avin_Av_13c = new TH1F("hAlphan_Lab_Avin_Av_13c",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAlphan_Lab_Avin_Av_18o = new TH1F("hAlphan_Lab_Avin_Av_18o",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAlphan_Lab_Avin_Ls_13c = new TH1F("hAlphan_Lab_Avin_Ls_13c",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAlphan_Lab_Avout_Av_13c = new TH1F("hAlphan_Lab_Avout_Av_13c",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAlphan_Lab_Avout_Av_18o = new TH1F("hAlphan_Lab_Avout_Av_18o",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hAr39 = new TH1F("hAr39",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi210 = new TH1F("hBi210",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi210_Av_Innerdust = new TH1F("hBi210_Av_Innerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi210_Av_Outerdust = new TH1F("hBi210_Av_Outerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi212 = new TH1F("hBi212",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi212po212 = new TH1F("hBi212po212",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214 = new TH1F("hBi214",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Av = new TH1F("hBi214_Av",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Av_Innerdust = new TH1F("hBi214_Av_Innerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Av_Outerdust = new TH1F("hBi214_Av_Outerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Exwater = new TH1F("hBi214_Exwater",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Exwater_r = new TH1F("hBi214_Exwater_r",";r;Events (mm^{-1} )",10000,0,10000);
    TH1F *hBi214_Hdropes = new TH1F("hBi214_Hdropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Hupropes = new TH1F("hBi214_Hupropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214_Intropes = new TH1F("hBi214_Intropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBi214po214 = new TH1F("hBi214po214",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBipo212_500 = new TH1F("hBipo212_500",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBipo212_500_Av_Innerdust = new TH1F("hBipo212_500_Av_Innerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBipo212_500_Intropes = new TH1F("hBipo212_500_Intropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBipo214_500 = new TH1F("hBipo214_500",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBipo214_500_Av_Innerdust = new TH1F("hBipo214_500_Av_Innerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hBipo214_500_Intropes = new TH1F("hBipo214_500_Intropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hC11 = new TH1F("hC11",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hC14 = new TH1F("hC14",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hK40 = new TH1F("hK40",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hK40_Intropes = new TH1F("hK40_Intropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hKr85 = new TH1F("hKr85",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPa234m = new TH1F("hPa234m",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPb210 = new TH1F("hPb210",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPb212 = new TH1F("hPb212",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPb214 = new TH1F("hPb214",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPmt_Betagammas = new TH1F("hPmt_Betagammas",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo210 = new TH1F("hPo210",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo210_Av_Innerdust = new TH1F("hPo210_Av_Innerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo210_Av_Outerdust = new TH1F("hPo210_Av_Outerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo212 = new TH1F("hPo212",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo214 = new TH1F("hPo214",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo216 = new TH1F("hPo216",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPo218 = new TH1F("hPo218",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hRa224 = new TH1F("hRa224",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hRa226 = new TH1F("hRa226",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hRa228 = new TH1F("hRa228",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hRn220 = new TH1F("hRn220",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hRn222 = new TH1F("hRn222",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTh228 = new TH1F("hTh228",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTh230 = new TH1F("hTh230",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTh232 = new TH1F("hTh232",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTh234 = new TH1F("hTh234",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208 = new TH1F("hTl208",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Av = new TH1F("hTl208_Av",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Av_Innerdust = new TH1F("hTl208_Av_Innerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Av_Outerdust = new TH1F("hTl208_Av_Outerdust",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Exwater = new TH1F("hTl208_Exwater",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Hdropes = new TH1F("hTl208_Hdropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Hupropes = new TH1F("hTl208_Hupropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl208_Intropes = new TH1F("hTl208_Intropes",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hTl210 = new TH1F("hTl210",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hU234 = new TH1F("hU234",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hU238 = new TH1F("hU238",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hUChain = new TH1F("hUChain",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hThChain = new TH1F("hThChain",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hnewTl = new TH1F("hnewTl",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hnew2Tl = new TH1F("hnew2Tl",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hnew3Tl = new TH1F("hnew3Tl",";nHits;Events (nHit^{-1})",2000,0,2000);
    
    TH1F *hBe7_Nue = new TH1F("hBe7_Nue",";nHits;Events (nHit^{-1})",2000,0,2000);
    //TH1F *hBe7_Numu = new TH1F("hBe7_Numu",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hB8_Nue = new TH1F("hB8_Nue",";nHits;Events (nHit^{-1})",2000,0,2000);
    //TH1F *hB8_Numu = new TH1F("hB8_Numu",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hCno_Nue = new TH1F("hCno_Nue",";nHits;Events (nHit^{-1})",2000,0,2000);
    //TH1F *hCno_Numu = new TH1F("hCno_Numu",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPep_Nue = new TH1F("hPep_Nue",";nHits;Events (nHit^{-1})",2000,0,2000);
    //TH1F *hPep_Numu = new TH1F("hPep_Numu",";nHits;Events (nHit^{-1})",2000,0,2000);
    TH1F *hPp_Nue = new TH1F("hPp_Nue",";nHits;Events (nHit^{-1})",2000,0,2000);
    
    TH1F *hExWatermine = new TH1F("hExWatermine",";nHits;Events (nHit^{-1})",2000,0,2000);
    
    //  ifstream finn("in_SolarBackgrounds.txt");
    string inPath = "/home/djauty/snoing/install/rattools-dev/GridTools/downloaded/";
    string str;
    double posx, posy, posz, posr,scale,energy;
    int nhits;
    int nhits_high = 10000;
    int nhits_low = 1;
    double radius_cut = 10000;
    TString inBackgrounds = inPath;
    /*  inBackgrounds += "SolarAc228_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAc228->Fill(nhits);
     }
     scale = 0.09741476636;
     cout<<"Ac228 "<<endl;
     hAc228->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAlphan_Lab_13c_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAlphan_Lab_13c->Fill(nhits);
     }
     scale = 0.006543290409;
     cout<<"Alphan_Lab_13c"<<endl;
     hAlphan_Lab_13c->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAlphan_Lab_Avin_Av_13c_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAlphan_Lab_Avin_Av_13c->Fill(nhits);
     }
     scale = 0.09138412493;
     cout<<"SolarAlphan_Lab_Avin_Av_13c"<<endl;
     //  if(scale>=1){cout<<"Okay SolarAlphan_Lab_Avin_Av_13c "<<nBackEntries<<" " <<366<<endl;}
     hAlphan_Lab_Avin_Av_13c->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAlphan_Lab_Avin_Av_18o_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAlphan_Lab_Avin_Av_18o->Fill(nhits);
     }
     scale = 0.0225358899;
     cout<<"SolarAlphan_Lab_Avin_Av_18o"<<endl;
     hAlphan_Lab_Avin_Av_18o->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAlphan_Lab_Avin_Ls_13c_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAlphan_Lab_Avin_Ls_13c->Fill(nhits);
     }
     scale = 0.1185996157;
     cout<<"SolarAlphan_Lab_Avin_Ls_13c"<<endl;
     hAlphan_Lab_Avin_Ls_13c->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAlphan_Lab_Avout_Av_13c_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAlphan_Lab_Avout_Av_13c->Fill(nhits);
     }
     scale = 0.05284965331;
     cout<<"SolarAlphan_Lab_Avout_Av_13c"<<endl;
     hAlphan_Lab_Avout_Av_13c->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAlphan_Lab_Avout_Av_18o_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAlphan_Lab_Avout_Av_18o->Fill(nhits);
     }
     scale = 0.01628340669;
     cout<<"SolarAlphan_Lab_Avout_Av_18o"<<endl;
     hAlphan_Lab_Avin_Av_18o->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarAr39.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     size_t nAr39Entries = nBackEntries;
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hAr39->Fill(nhits);
     }
     scale = 0.5344337832;
     cout<<"SolarAr39 = "<<endl;
     hAr39->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi210_Av_Innerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi210_Av_Innerdust->Fill(nhits);
     }
     
     scale = 1.06E+04;
     cout<<"SolarBi210_Av_Innerdust"<<endl;
     hBi210->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi210_Av_Outerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi210_Av_Outerdust->Fill(nhits);
     }
     
     scale = 3.78E+03;
     if(scale<1){cout<<"SolarBi210_av_Outerdust"<<endl;}
     hBi210_Av_Outerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi210_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi210->Fill(nhits);
     }
     scale = 3.78E+03;
     cout<<"SolarBi210"<<endl;
     hBi210_Av_Outerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi212po212_r1_s0_p2.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi212po212->Fill(nhits);
     }
     scale = 0.1138214639;
     cout<<"SolarBi212_po212"<<endl;
     hBi212po212->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi212_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi212->Fill(nhits);
     }
     scale = 0.06822106667;
     cout<<"SolarBi212"<<endl;
     hBi212->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi214_Av_Innerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     size_t nBi214_av_innerEntries = nBackEntries;
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi214_Av_Innerdust->Fill(nhits);
     }
     
     scale = 0.576;
     cout<<"SolarBi214_AV_Innerdust "<<endl;
     hBi214_Av_Innerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi214_Av_Outerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi214_Av_Outerdust->Fill(nhits);
     }
     scale = 1.23E+00;
     cout<<"SolarBi214_AV_Outerdust"<<endl;
     hBi214_Av_Outerdust->Scale(scale);
     
    inBackgrounds = inPath;
    inBackgrounds += "SolarBi214_Av.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hBi214_Av->Fill(nhits);
    }
    scale = 1.18;
    cout<<"SolarBi214_AV"<<endl;
    hBi214_Av->Scale(scale);
    */
    cout<<"here"<<endl;
    inBackgrounds = inPath;
    inBackgrounds += "SolarBi214_Exwater.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hBi214_Exwater->Fill(nhits);
        hBi214_Exwater_r->Fill(posr);
    }
    scale = 1.14E+00;
    cout<<"SolarBi214_Exwater"<<endl;
    hBi214_Exwater->Scale(scale);
    hBi214_Exwater_r->Scale(scale);
    /*
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi214_Hdropes.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi214_Hdropes->Fill(nhits);
     }
     scale = 1.17;
     cout<<"SolarBi214_Hdropes"<<endl;
     hBi214_Hdropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi214_Hupropes_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi214_Hupropes->Fill(nhits);
     }
     scale = 7.58E-01;
     cout<<"SolarBi214_Hupropes"<<endl;
     hBi214_Hupropes->Scale(scale);
     
    inBackgrounds = inPath;
    inBackgrounds += "SolarBi214_Intropes_r1_s0_p1.ntuple.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hBi214_Intropes->Fill(nhits);
    }
    scale = 0.5833427288;
    cout<<"SolarBi214_Intropes"<<endl;
    hBi214_Intropes->Scale(scale);
    
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi214po214_r1_s0_p2.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi214po214->Fill(nhits);
     }
     scale = 0.8174450586;
     cout<<"SolarBi214po214"<<endl;
     hBi214po214->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBi214_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBi214->Fill(nhits);
     }
     scale = 0.9768211392;
     cout<<"SolarBi214"<<endl;
     hBi214->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBipo212-500_Av_Innerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBipo212_500_Av_Innerdust->Fill(nhits);
     }
     scale = 0.9428201291;
     cout<<"SolarBipo212-500_Av_Innerdust"<<endl;
     hBipo212_500_Av_Innerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBipo212-500_Intropes_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBipo212_500_Intropes->Fill(nhits);
     }
     scale = 0.02545561069;
     cout<<"SolarBipo212-500_Intropes"<<endl;
     hBipo212_500_Intropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBipo212-500_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBipo212_500->Fill(nhits);
     }
     scale = 0.05956255644;
     cout<<"SolarBipo212-500"<<endl;
     hBipo212_500->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBipo214-500_Av_Innerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBipo214_500_Av_Innerdust->Fill(nhits);
     }
     scale = 0.01244823591;
     cout<<"solarBipo214-500_Av_Innerdust"<<endl;
     hBipo214_500_Av_Innerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBipo214-500_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBipo214_500->Fill(nhits);
     }
     scale = 0.05956255644;
     cout<<"solarBipo214-500"<<endl;
     hBipo214_500->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarBipo214-500_Intropes_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBipo214_500_Intropes->Fill(nhits);
     }
     scale = 0.5833427288;
     cout<<"solarBipo214-500-Intropes"<<endl;
     hBipo214_500_Intropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarC11_r1_s0_p0.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hC11->Fill(nhits);
     }
     scale = 0.201;
     cout<<"solarC11"<<endl;
     hC11->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarC14.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hC14->Fill(nhits);
     }
     scale = 1.93E+02;
     cout<<"solarC14"<<endl;
     hC14->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarK40_Intropes_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hK40_Intropes->Fill(nhits);
     }
     scale = 0.8481619124;
     cout<<"solarK40_Intropes"<<endl;
     hK40_Intropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarK40_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hK40->Fill(nhits);
     }
     scale = 0.8027433827;
     cout<<"solarK40"<<endl;
     hK40->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarKr85.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hKr85->Fill(nhits);
     }
     scale = 0.6128917159;
     cout<<"solarKr85"<<endl;
     hKr85->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPa234m_r0_s0_p2.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPa234m->Fill(nhits);
     }
     scale = 0.489951297;
     cout<<"solarPa234m"<<endl;
     hPa234m->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPb210.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPb210->Fill(nhits);
     }
     scale = 3.86E+02;
     cout<<"solarPb210"<<endl;
     hPb210->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPb212_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPb212->Fill(nhits);
     }
     scale = 0.02727122085;
     cout<<"solarPb212"<<endl;
     hPb212->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPb214_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPb214->Fill(nhits);
     }
     scale = 0.3262910948;
     cout<<"solarPb214"<<endl;
     hPb214->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPmt_Betagammas.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPmt_Betagammas->Fill(nhits);
     }
     scale = 8.41E-01;
     cout<<"solarPmt_Betagammas"<<endl;
     hPmt_Betagammas->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo210_Av_Innerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPo210_Av_Innerdust->Fill(nhits);
     }
     scale = 1.16E+04;
     cout<<"solarPo210_Av_Innerdust"<<endl;
     hPo210_Av_Innerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo210_Av_Outerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPo210_Av_Outerdust->Fill(nhits);
     }
     scale = 3.90E+03;
     cout<<"solarPo210_Av_Outerdust"<<endl;
     hPo210_Av_Outerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo210.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPb210->Fill(nhits);
     }
     scale = 7.30E+01;
     cout<<"solarPo210"<<endl;
     hPo210->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo212_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPo212->Fill(nhits);
     }
     scale = 0.13601307;
     cout<<"solarPo212"<<endl;
     hPo212->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo214_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPo214->Fill(nhits);
     }
     scale = 0.489951297;
     cout<<"solarPo214"<<endl;
     hPo214->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo216_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPo216->Fill(nhits);
     }
     scale = 0.06822106667;
     cout<<"solarPo216"<<endl;
     hPo216->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarPo218_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hPo218->Fill(nhits);
     }
     scale = 0.3262910948;
     cout<<"solarPo218"<<endl;
     hPo218->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarRa224_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hRa224->Fill(nhits);
     }
     scale = 0.04543293725;
     cout<<"solarRa224"<<endl;
     hRa224->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarRa226_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hRa226->Fill(nhits);
     }
     scale = 0.2449756485;
     cout<<"solarRa226"<<endl;
     hRa226->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarRa228_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hRa228->Fill(nhits);
     }
     scale = 0.006666485244;
     cout<<"solarRa228"<<endl;
     hRa228->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarRn220_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hRn220->Fill(nhits);
     }
     scale = 0.04543293725;
     cout<<"solarRn220"<<endl;
     hRn220->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarRn222_r1_s0_p2.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hRn222->Fill(nhits);
     }
     scale = 0.3262910948;
     cout<<"solarRn222"<<endl;
     hRn222->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTh228_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTh228->Fill(nhits);
     }
     scale = 0.03411053333;
     cout<<"solarTh228"<<endl;
     hTh228->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTh230_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTh230->Fill(nhits);
     }
     scale = 0.2449756485;
     cout<<"solarTh230"<<endl;
     hTh230->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTh232_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTh232->Fill(nhits);
     }
     scale = 0.02727122085;
     cout<<"solarTh232"<<endl;
     hTh232->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTh234_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTh234->Fill(nhits);
     }
     scale = 0.04897967869;
     cout<<"solarTh234"<<endl;
     hTh234->Scale(scale);
     
    inBackgrounds = inPath;
    inBackgrounds += "SolarTl208_Av.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hTl208_Av_Innerdust->Fill(nhits);
    }
    scale = 5.88E-01;
    cout<<"solarTl208_Av"<<endl;
    hTl208_Av->Scale(scale);
    
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl208_Av_Innerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl208_Av_Innerdust->Fill(nhits);
     }
     scale = 2.75E-01;
     cout<<"solarTl208_innerdust"<<endl;
     hTl208_Av_Innerdust->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl208_Av_Outerdust.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl208_Av_Outerdust->Fill(nhits);
     }
     scale = 1.16E+00;
     cout<<"solarTl208_outerdust"<<endl;
     hTl208_Av_Outerdust->Scale(scale);
     */
    inBackgrounds = inPath;
    inBackgrounds += "SolarTl208_Exwater.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hTl208_Exwater->Fill(nhits);
    }
    scale = 1.06E+00;
    cout<<"solarTl208_Exwater"<<endl;
    hTl208_Exwater->Scale(scale);
    
    inBackgrounds = inPath;
    inBackgrounds += "Tl208_exwater_full.ntuple.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hnewTl->Fill(nhits);
    }
    scale = 1.00E+00;
    cout<<"Tl208_new"<<endl;
    hnewTl->Scale(scale);
    
    inBackgrounds = inPath;
    inBackgrounds += "Tl208_exwater2_full.ntuple.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hnew2Tl->Fill(nhits);
    }
    scale = 1.00E+00;
    cout<<"Tl208_new2"<<endl;
    hnew2Tl->Scale(scale);
    
    inBackgrounds = inPath;
    inBackgrounds += "Tl208_exwater3_full.ntuple.root";
    TFile *f = TFile::Open(inBackgrounds, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nBackEntries = t->GetEntries();
    for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hnew3Tl->Fill(nhits);
    }
    scale = 250000/249500.;
    cout<<"Tl208_new3"<<endl;
    hnew3Tl->Scale(scale);
    
    /*
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl208_Hdropes.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl208_Hdropes->Fill(nhits);
     }
     scale = 1.24;
     cout<<"solarTl208_Hdropes"<<endl;
     hTl208_Hdropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl208_Hupropes.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl208_Hupropes->Fill(nhits);
     }
     scale = 0.352;
     cout<<"solarTl208_Hupropes"<<endl;
     hTl208_Hupropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl208_Intropes_r1_s0_p1.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl208_Intropes->Fill(nhits);
     }
     scale = 0.08336284937;
     cout<<"solarTl208_Intropes"<<endl;
     hTl208_Intropes->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl208_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl208->Fill(nhits);
     }
     scale = 0.08211167187;
     cout<<"solarTl208"<<endl;
     hTl208->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarTl210_r1_s0_p2.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hTl210->Fill(nhits);
     }
     scale = 0.0002856738017;
     cout<<"solarTl210"<<endl;
     hTl210->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarU234_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hU234->Fill(nhits);
     }
     scale = 0.2449756485;
     cout<<"solarU234"<<endl;
     hU234->Scale(scale);
     
     inBackgrounds = inPath;
     inBackgrounds += "SolarU238_r1_s0_p3.ntuple.root";
     TFile *f = TFile::Open(inBackgrounds, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nBackEntries = t->GetEntries();
     for(size_t iEntry = 0;iEntry<nBackEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hU238->Fill(nhits);
     }
     scale = 0.1958569497;
     cout<<"solarU238"<<endl;
     hU238->Scale(scale);
     
     
     ///////signal
     TString inSignal = inPath;
     inSignal += "SolarBe7.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     size_t Be7Enties = nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hBe7_Nue->Fill(nhits);
     }
     /*    scale = nSigEntries/nSigEntries;
     if(scale<1){cout<<"solarBe7nue"<<endl;}
     hBe7_Nue->Scale(1/scale);
     TString inSignal = inPath;
     inSignal += "SolarBe7_Solar_Numu_r1_s0_p0.ntuple.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     Be7Enties += nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     // hBe7_Numu->Fill(nhits);
     hBe7_Nue->Fill(nhits);
     }
     scale = nSigEntries/nSigEntries;
     if(scale<1){cout<<"solarBe7num"<<endl;}
     hBe7_Numu->Scale(1/scale);*/
    /*    scale = nSigEntries/(309452.*18);
     if(scale<1){cout<<"solarBe7nue"<<endl;}
     hBe7_Nue->Scale(1/scale);
     scale = nSigEntries/(220805+49258.6);
     if(scale<1){cout<<"solarBe7nue"<<endl;}
     hBe7_Nue->Scale(1/scale);
     
     TString inSignal = inPath;
     inSignal += "SolarB8.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     size_t B8Entries = nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hB8_Nue->Fill(nhits);
     }
     /*   scale = nSigEntries/3871.;
     if(scale<1){cout<<"solarB8nue"<<endl;}
     hB8_Nue->Scale(1/scale);
     
     TString inSignal = inPath;
     inSignal += "SolarB8_Solar_Numu_r1_s0_p0.ntuple.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     B8Entries += nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     //        hB8_Numu->Fill(nhits);
     hB8_Nue->Fill(nhits);
     }
     /*   scale = nSigEntries/nSigEntries;
     if(scale<1){cout<<"solarB8num"<<endl;}
     hB8_Numu->Scale(1/scale);*/
    /*    scale = B8Entries/(3871.*4);
     if(scale<1){cout<<"solarB8"<<endl;}
     hB8_Nue->Scale(1/scale);
     scale = B8Entries/(2900.13+517.59);
     if(scale<1){cout<<"solarB8"<<endl;}
     hB8_Nue->Scale(1/scale);
     */
    TString inSignal = inPath;
    inSignal += "SolarCno.root";
    TFile *f = TFile::Open(inSignal, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nSigEntries = t->GetEntries();
    size_t CnoEntries = nSigEntries;
    for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hCno_Nue->Fill(nhits);
    }
    /*    scale = nSigEntries/36503.;
     if(scale<1){cout<<"solarCnonue"<<endl;}
     hCno_Nue->Scale(1/scale);
     
     TString inSignal = inPath;
     inSignal += "SolarCno_Solar_Numu_r2_s0_p0.ntuple.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     CnoEntries += nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     //    hCno_Numu->Fill(nhits);
     hCno_Nue->Fill(nhits);
     }
     /*   scale = nSigEntries/nSigEntries;
     if(scale<1){cout<<"solarCnonum"<<endl;}
     hCno_Numu->Scale(1/scale);*/
/*    scale = nSigEntries/(36503.*3);
    if(scale<1){cout<<"solarCno"<<endl;}
    hCno_Nue->Scale(1/scale);*/
    
    scale = nSigEntries/(13866.5+2957.31+11859.8+2690.38+348.881+74.3682);
    if(scale<1){cout<<"solarCno"<<endl;}
    hCno_Nue->Scale(1/scale);
    
    TString inSignal = inPath;
    inSignal += "SolarPep.root";
    TFile *f = TFile::Open(inSignal, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nSigEntries = t->GetEntries();
    size_t pepEntries =nSigEntries;
    for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hPep_Nue->Fill(nhits);
    }
    /*  scale = nSigEntries/18021;
     if(scale<1){cout<<"solarPepnue"<<endl;}
     hPep_Nue->Scale(1/scale);
     
     TString inSignal = inPath;
     inSignal += "SolarPep_Solar_Numu_r1_s0_p0.ntuple.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     pepEntries +=nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     // hPep_Numu->Fill(nhits);
     hPep_Nue->Fill(nhits);
     }
     /*    scale = nSigEntries/nSigEntries;
     if(scale<1){cout<<"solarPepnum"<<endl;}
     hPep_Numu->Scale(1/scale);*/
   /* scale = pepEntries/(18021.*7);
    if(scale<1){cout<<"solarPep"<<endl;}
    hPep_Nue->Scale(1/scale);*/
    scale = pepEntries/(13081.3+2659.4);
    if(scale<1){cout<<"solarPep"<<endl;}
    hPep_Nue->Scale(1/scale);
    
    
   /* TString inSignal = inPath;
    inSignal += "SolarPp.root";
    TFile *f = TFile::Open(inSignal, "READ");
    TTree *t = (TTree*)f->Get("output");
    t->SetBranchAddress("nhits",&nhits);
    t->SetBranchAddress("posx",&posx);
    t->SetBranchAddress("posy",&posy);
    t->SetBranchAddress("posz",&posz);
    t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
    size_t nSigEntries = t->GetEntries();
    size_t ppEntries =nSigEntries;
    for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
        t->GetEntry(iEntry);
        if(nhits>nhits_high){continue;}
        if(nhits<nhits_low){continue;}
        if(posr>radius_cut){continue;}
        hPp_Nue->Fill(nhits);
    }
    /*  scale = nSigEntries/18021;
     if(scale<1){cout<<"solarPepnue"<<endl;}
     hPep_Nue->Scale(1/scale);
     
     TString inSignal = inPath;
     inSignal += "SolarPep_Solar_Numu_r1_s0_p0.ntuple.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nSigEntries = t->GetEntries();
     pepEntries +=nSigEntries;
     for(size_t iEntry = 0;iEntry<nSigEntries ;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     // hPep_Numu->Fill(nhits);
     hPep_Nue->Fill(nhits);
     }
     /*    scale = nSigEntries/nSigEntries;
     if(scale<1){cout<<"solarPepnum"<<endl;}
     hPep_Numu->Scale(1/scale);*/
    /* scale = pepEntries/(18021.*7);
     if(scale<1){cout<<"solarPep"<<endl;}
     hPep_Nue->Scale(1/scale);*/
   /* scale = ppEntries/(588269+166375.);
    if(scale<1){cout<<"solarPp"<<endl;}
    hPp_Nue->Scale(1/scale);
    
    /*
     TString inSignal = "/home/djauty/codes/mac_files/kalpana/";
     inSignal += "Bi214exwater.root";
     TFile *f = TFile::Open(inSignal, "READ");
     TTree *t = (TTree*)f->Get("output");
     t->SetBranchAddress("nhits",&nhits);
     t->SetBranchAddress("posx",&posx);
     t->SetBranchAddress("posy",&posy);
     t->SetBranchAddress("posz",&posz);
     t->SetBranchAddress("posr",&posr);t->SetBranchAddress("energy",&energy);
     size_t nExWaterEntries = t->GetEntries();
     size_t CnoEntries = nSigEntries;
     for(size_t iEntry = 0;iEntry<nExWaterEntries ;++iEntry){
     t->GetEntry(iEntry);
     if(nhits>nhits_high){continue;}
     if(nhits<nhits_low){continue;}
     if(posr>radius_cut){continue;}
     hExWatermine->Fill(nhits);
     }
     scale = 1.32e8/50000;
     cout<<"Biwater mine"<<endl;
     hExWatermine->Scale(scale);
     */
    hBe7_Nue ->SetLineColor(kRed);
    //  hBe7_Numu->SetLineColor(kBlue);
    hB8_Nue ->SetLineColor(kOrange);
    //  hB8_Numu->SetLineColor(kBlue);
    hCno_Nue ->SetLineColor(kGreen);
    //  hCno_Numu ->SetLineColor(kBlue);
    hPep_Nue ->SetLineColor(kBlue);
    //  hPep_Numu->SetLineColor(kBlue);
    hSigTotal->SetLineColor(kBlue);
    hSigTotal->SetLineWidth(3);
    
    hSigTotal->Add(hBe7_Nue );
    // hSigTotal->Add(hBe7_Numu);
    hSigTotal->Add(hB8_Nue );
    // hSigTotal->Add(hB8_Numu);
    hSigTotal->Add(hCno_Nue );
    // hSigTotal->Add(hCno_Numu );
    hSigTotal->Add(hPep_Nue );
    // hSigTotal->Add(hPep_Numu);
    
    hThChain->Add(hTh232);
    hThChain->Add(hTh228);
    hThChain->Add(hAc228);
    hThChain->Add(hRa228);
    hThChain->Add(hRa224);
    hThChain->Add(hRn220);
    hThChain->Add(hPo216);
    hThChain->Add(hPb212);
    hThChain->Add(hBi212);
    hThChain->Add(hPo212);
    hThChain->Add(hTl208);
    hThChain->SetLineColor(kBlue);
    hThChain->SetLineStyle(3);
    hThChain->SetLineWidth(2);
    
    hAr39->SetLineColor(14);
    hAr39->SetLineStyle(3);
    hAr39->SetLineWidth(2);
    
    hBi210->SetLineColor(3);
    hBi210->SetLineStyle(3);
    hBi210->SetLineWidth(2);
    
    hUChain->Add(hU238);
    hUChain->Add(hTh234);
    hUChain->Add(hPa234m);
    hUChain->Add(hU234);
    hUChain->Add(hTh230);
    hUChain->Add(hRa226);
    hUChain->Add(hRn222);
    hUChain->Add(hPo218);
    hUChain->Add(hPb214);
    hUChain->Add(hBi214);
    hUChain->Add(hPo214);
    hUChain->Add(hTl210);
    hUChain->Add(hPb210);
    hUChain->Add(hBi210);
    // hUChain->Add(hPo210);
    hUChain->SetLineColor(kRed);
    hUChain->SetLineStyle(3);
    hUChain->SetLineWidth(2);
    
    hC11->SetLineColor(8);
    hC11->SetLineStyle(3);
    hC11->SetLineWidth(2);
    
    hKr85->SetLineColor(kOrange);
    hKr85->SetLineStyle(3);
    hKr85->SetLineWidth(2);
    
    hPo210->SetLineColor(kCyan);
    hPo210->SetLineStyle(3);
    hPo210->SetLineWidth(2);
    
    /* hAc228->SetLineColor(kRed);
     hAlphan_Lab_13c->SetLineColor(kRed);
     hAlphan_Lab_Avin_Av_13c->SetLineColor(kRed);
     hAlphan_Lab_Avin_Av_18o->SetLineColor(kRed);
     hAlphan_Lab_Avin_Ls_13c->SetLineColor(kRed);
     hAlphan_Lab_Avout_Av_13c->SetLineColor(kRed);
     hAlphan_Lab_Avout_Av_18o->SetLineColor(kRed);
     
     //  hBi210->SetLineColor(kRed);
     hBi210_Av_Innerdust->SetLineColor(kRed);
     hBi210_Av_Outerdust->SetLineColor(kRed);
     hBi212->SetLineColor(kRed);
     hBi212po212->SetLineColor(kRed);
     hBi214->SetLineColor(kRed);
     hBi214_Av->SetLineColor(kRed);
     hBi214_Av_Innerdust->SetLineColor(kRed);
     hBi214_Av_Outerdust->SetLineColor(kRed);
     hBi214_Exwater->SetLineColor(kRed);
     hBi214_Hdropes->SetLineColor(kRed);
     hBi214_Hupropes->SetLineColor(kRed);
     hBi214_Intropes->SetLineColor(kRed);
     hBi214po214->SetLineColor(kRed);
     hBipo212_500->SetLineColor(kRed);
     hBipo212_500_Av_Innerdust->SetLineColor(kRed);
     hBipo212_500_Intropes->SetLineColor(kRed);
     hBipo214_500->SetLineColor(kRed);
     hBipo214_500_Av_Innerdust->SetLineColor(kRed);
     hBipo214_500_Intropes->SetLineColor(kRed);
     // hC11->SetLineColor(kRed);
     hC14->SetLineColor(kRed);
     hK40->SetLineColor(kRed);
     hK40_Intropes->SetLineColor(kRed);
     // hKr85->SetLineColor(kRed);
     hPa234m->SetLineColor(kRed);
     hPb210->SetLineColor(kRed);
     hPb212->SetLineColor(kRed);
     hPb214->SetLineColor(kRed);
     hPmt_Betagammas->SetLineColor(kRed);
     // hPo210->SetLineColor(kRed);
     hPo210_Av_Innerdust->SetLineColor(kRed);
     hPo210_Av_Outerdust->SetLineColor(kRed);
     hPo212->SetLineColor(kRed);
     hPo214->SetLineColor(kRed);
     hPo216->SetLineColor(kRed);
     hPo218->SetLineColor(kRed);
     hRa224->SetLineColor(kRed);
     hRa226->SetLineColor(kRed);
     hRa228->SetLineColor(kRed);
     hRn220->SetLineColor(kRed);
     hRn222->SetLineColor(kRed);
     hTh228->SetLineColor(kRed);
     hTh230->SetLineColor(kRed);
     hTh232->SetLineColor(kRed);
     hTh234->SetLineColor(kRed);
     hTl208->SetLineColor(kRed);
     hTl208_Av->SetLineColor(kRed);
     hTl208_Av_Innerdust->SetLineColor(kRed);
     hTl208_Av_Outerdust->SetLineColor(kRed);
     hTl208_Exwater->SetLineColor(kRed);
     hTl208_Hdropes->SetLineColor(kRed);
     hTl208_Hupropes->SetLineColor(kRed);
     hTl208_Intropes->SetLineColor(kRed);
     hTl210->SetLineColor(kRed);
     hU234->SetLineColor(kRed);
     hU238->SetLineColor(kRed);*/
    hBackTotal->SetLineColor(kRed);
    hBackTotal->SetLineWidth(3);
    
    hBackTotal->Add(hAc228);
    hBackTotal->Add(hAlphan_Lab_13c);
    hBackTotal->Add(hAlphan_Lab_Avin_Av_13c);
    hBackTotal->Add(hAlphan_Lab_Avin_Av_18o);
    hBackTotal->Add(hAlphan_Lab_Avin_Ls_13c);
    hBackTotal->Add(hAlphan_Lab_Avout_Av_13c);
    hBackTotal->Add(hAlphan_Lab_Avout_Av_18o);
    hBackTotal->Add(hAr39);
    hBackTotal->Add(hBi210);
    hBackTotal->Add(hBi210_Av_Innerdust);
    hBackTotal->Add(hBi210_Av_Outerdust);
    hBackTotal->Add(hBi212);
    hBackTotal->Add(hBi212po212);
    hBackTotal->Add(hBi214);
    hBackTotal->Add(hBi214_Av);
    hBackTotal->Add(hBi214_Av_Innerdust);
    hBackTotal->Add(hBi214_Av_Outerdust);
    hBackTotal->Add(hBi214_Exwater);
    hBackTotal->Add(hBi214_Hdropes);
    hBackTotal->Add(hBi214_Hupropes);
    hBackTotal->Add(hBi214_Intropes);
    hBackTotal->Add(hBi214po214);
    hBackTotal->Add(hBipo212_500);
    hBackTotal->Add(hBipo212_500_Av_Innerdust);
    hBackTotal->Add(hBipo212_500_Intropes);
    hBackTotal->Add(hBipo214_500);
    hBackTotal->Add(hBipo214_500_Av_Innerdust);
    hBackTotal->Add(hBipo214_500_Intropes);
    hBackTotal->Add(hC11);
    hBackTotal->Add(hC14);
    hBackTotal->Add(hK40);
    hBackTotal->Add(hK40_Intropes);
    hBackTotal->Add(hKr85);
    hBackTotal->Add(hPa234m);
    hBackTotal->Add(hPb210);
    hBackTotal->Add(hPb212);
    hBackTotal->Add(hPb214);
    hBackTotal->Add(hPmt_Betagammas);
    hBackTotal->Add(hPo210);
    hBackTotal->Add(hPo210_Av_Innerdust);
    hBackTotal->Add(hPo210_Av_Outerdust);
    hBackTotal->Add(hPo212);
    hBackTotal->Add(hPo214);
    hBackTotal->Add(hPo216);
    hBackTotal->Add(hPo218);
    hBackTotal->Add(hRa224);
    hBackTotal->Add(hRa226);
    hBackTotal->Add(hRa228);
    hBackTotal->Add(hRn220);
    hBackTotal->Add(hRn222);
    hBackTotal->Add(hTh228);
    hBackTotal->Add(hTh230);
    hBackTotal->Add(hTh232);
    hBackTotal->Add(hTh234);
    hBackTotal->Add(hTl208);
    hBackTotal->Add(hTl208_Av);
    hBackTotal->Add(hTl208_Av_Innerdust);
    hBackTotal->Add(hTl208_Av_Outerdust);
    hBackTotal->Add(hTl208_Exwater);
    hBackTotal->Add(hTl208_Hdropes);
    hBackTotal->Add(hTl208_Hupropes);
    hBackTotal->Add(hTl208_Intropes);
    hBackTotal->Add(hTl210);
    hBackTotal->Add(hU234);
    hBackTotal->Add(hU238);
    
    
    hBackExt->Add(hAlphan_Lab_Avout_Av_13c);
    hBackExt->Add(hAlphan_Lab_Avout_Av_18o);
    hBackExt->Add(hBi210_Av_Outerdust);
    hBackExt->Add(hBi214_Av);
    hBackExt->Add(hBi214_Av_Outerdust);
    hBackExt->Add(hBi214_Exwater);
    hBackExt->Add(hBi214_Hdropes);
    hBackExt->Add(hBi214_Hupropes);
    hBackExt->Add(hBipo212_500_Intropes);
    hBackExt->Add(hBipo214_500_Intropes);
    hBackExt->Add(hPmt_Betagammas);
    hBackExt->Add(hPo210_Av_Outerdust);
    hBackExt->Add(hTl208_Av);
    hBackExt->Add(hTl208_Av_Innerdust);
    hBackExt->Add(hTl208_Av_Outerdust);
    hBackExt->Add(hTl208_Exwater);
    hBackExt->Add(hTl208_Hdropes);
    hBackExt->Add(hTl208_Hupropes);
    hBackExt->SetLineColor(kMagenta);
    hBackExt->SetLineWidth(3);
    
    hBackOther->Add(hAlphan_Lab_Avin_Ls_13c);
    hBackOther->Add(hBi210_Av_Innerdust);
    hBackOther->Add(hBi212po212);
    hBackOther->Add(hBi214_Av_Innerdust);
    hBackOther->Add(hBi214_Intropes);
    hBackOther->Add(hBi214po214);
    hBackOther->Add(hBipo212_500);
    hBackOther->Add(hBipo212_500_Av_Innerdust);
    hBackOther->Add(hBipo214_500);
    hBackOther->Add(hBipo214_500_Av_Innerdust);
    hBackOther->Add(hC14);
    hBackOther->Add(hK40);
    hBackOther->Add(hK40_Intropes);
    hBackOther->Add(hPmt_Betagammas);
    hBackOther->Add(hPo210_Av_Innerdust);
    hBackOther->Add(hTl208_Intropes);
    hBackOther->SetLineWidth(3);
    hBackOther->SetLineColor(kCyan);
    
    hTotal->Add(hAc228);//Thchain
    hTotal->Add(hAlphan_Lab_13c);//BackEx
    hTotal->Add(hAlphan_Lab_Avin_Av_13c);//BackEx
    hTotal->Add(hAlphan_Lab_Avin_Av_18o);//BackEx
    hTotal->Add(hAlphan_Lab_Avin_Ls_13c);//other
    hTotal->Add(hAlphan_Lab_Avout_Av_13c);//BackEx
    hTotal->Add(hAlphan_Lab_Avout_Av_18o);//BackEx
    hTotal->Add(hAr39);//Ar39
    hTotal->Add(hB8_Nue);//signal
    // hTotal->Add(hB8_Numu);
    hTotal->Add(hBe7_Nue);//signal
    // hTotal->Add(hBe7_Numu);
    hTotal->Add(hBi210);//Uchain
    hTotal->Add(hBi210_Av_Innerdust);//other
    hTotal->Add(hBi210_Av_Outerdust);//BackEx
    hTotal->Add(hBi212);//Thchain
    hTotal->Add(hBi212po212);//other
    hTotal->Add(hBi214);//Uchain
    hTotal->Add(hBi214_Av);//BackEx
    hTotal->Add(hBi214_Av_Innerdust);//other
    hTotal->Add(hBi214_Av_Outerdust);//BackEx
    hTotal->Add(hBi214_Exwater);//BackEx
    hTotal->Add(hBi214_Hdropes);//BackEx
    hTotal->Add(hBi214_Hupropes);//BackEx
    hTotal->Add(hBi214_Intropes);//other
    hTotal->Add(hBi214po214);//other
    hTotal->Add(hBipo212_500);//other
    hTotal->Add(hBipo212_500_Av_Innerdust);//other
    hTotal->Add(hBipo212_500_Intropes);//BackEx
    hTotal->Add(hBipo214_500);//other
    hTotal->Add(hBipo214_500_Av_Innerdust);//other
    hTotal->Add(hBipo214_500_Intropes);//BackEx
    hTotal->Add(hC11);//C11
    hTotal->Add(hC14);//other
    hTotal->Add(hCno_Nue);//sig
    //hTotal->Add(hCno_Numu);
    hTotal->Add(hK40);//other
    hTotal->Add(hK40_Intropes);//other
    hTotal->Add(hKr85);//Kr85
    hTotal->Add(hPa234m);//U chain
    hTotal->Add(hPb210);//U chain
    hTotal->Add(hPb212);//Th chain
    hTotal->Add(hPb214);//U chain
    hTotal->Add(hPep_Nue);//sig
    //hTotal->Add(hPep_Numu);
    hTotal->Add(hPmt_Betagammas);//other
    hTotal->Add(hPo210);//Po210
    hTotal->Add(hPo210_Av_Innerdust);//other
    hTotal->Add(hPo210_Av_Outerdust);//BackEx
    hTotal->Add(hPo212);//Th chain
    hTotal->Add(hPo214);//Uchain
    hTotal->Add(hPo216);//Thchain
    hTotal->Add(hPo218);//Uchain
    hTotal->Add(hRa224);//Thchain
    hTotal->Add(hRa226);//Uchain
    hTotal->Add(hRa228);//Thchain
    hTotal->Add(hRn220);//Thchain
    hTotal->Add(hRn222);//Uchain
    hTotal->Add(hTh228);//Thchain
    hTotal->Add(hTh230);//Uchain
    hTotal->Add(hTh232);//Thchain
    hTotal->Add(hTh234);//Uchain
    hTotal->Add(hTl208);//Thchain
    hTotal->Add(hTl208_Av);//BackExt
    hTotal->Add(hTl208_Av_Innerdust);//BackEx
    hTotal->Add(hTl208_Av_Outerdust);//BackEx
    hTotal->Add(hTl208_Exwater);//BackEx
    hTotal->Add(hTl208_Hdropes);//BackEx
    hTotal->Add(hTl208_Hupropes);//BackEx
    hTotal->Add(hTl208_Intropes);//other
    hTotal->Add(hTl210);//Uchain
    hTotal->Add(hU234);//Uchian
    hTotal->Add(hU238);//Uchain
    
    hTotal->SetLineColor(kBlack);
    hTotal->SetLineWidth(3);
    hTotal->SetMinimum(0.8);
    //  hTotal->SetMaximum(10e10);
    
    TCanvas *c1 = new TCanvas("c1");
    hTotal->Draw();
    // hBackTotal->Draw("same");
    //  hSigTotal->Draw("same");
    hPep_Nue->Draw("same");
    hB8_Nue->Draw("same");
    hBe7_Nue->Draw("same");
    hCno_Nue->Draw("same");
    hThChain->Draw("same");
    hAr39->Draw("sames");
    hBi210->Draw("sames");
    hUChain->Draw("same");
    hC11->Draw("sames");
    hKr85->Draw("sames");
    hPo210->Draw("sames");
    hBackExt->Draw("same");
    hBackOther->Draw("same");
    
    hExWatermine->SetLineColor(kYellow);
    hExWatermine->SetLineWidth(2);
    hExWatermine->SetLineStyle(5);
    
    /* hAc228->Draw("same");
     hAlphan_Lab_13c->Draw("same");
     hAlphan_Lab_Avin_Av_13c->Draw("same");
     hAlphan_Lab_Avin_Av_18o->Draw("same");
     hAlphan_Lab_Avin_Ls_13c->Draw("same");
     hAlphan_Lab_Avout_Av_13c->Draw("same");
     hAlphan_Lab_Avout_Av_18o->Draw("same");
     hAr39->Draw("same");
     // hB8_Numu->Draw("same");
     // hBe7_Numu->Draw("same");
     hBi210->Draw("same");
     hBi210_Av_Innerdust->Draw("same");
     hBi210_Av_Outerdust->Draw("same");
     hBi212->Draw("same");
     hBi212po212->Draw("same");
     hBi214->Draw("same");
     */
    hBi214_Av->SetLineStyle(3);
    hBi214_Av->SetLineColor(kRed);
    hBi214_Av->Draw("same");
    /*    hBi214_Av_Innerdust->Draw("same");
     hBi214_Av_Outerdust->Draw("same");
     */hBi214_Exwater->SetLineColor(kBlue);
    hBi214_Exwater->SetLineStyle(3);
    hBi214_Exwater->Draw("same");
    /*   hBi214_Hdropes->Draw("same");
     hBi214_Hupropes->Draw("same");
     */hBi214_Intropes->SetLineColor(kOrange);
    hBi214_Intropes->SetLineStyle(3);
    hBi214_Intropes->Draw("same");
    /*    hBi214po214->Draw("same");
     hBipo212_500->Draw("same");
     hBipo212_500_Av_Innerdust->Draw("same");
     hBipo212_500_Intropes->Draw("same");
     hBipo214_500->Draw("same");
     hBipo214_500_Av_Innerdust->Draw("same");
     hBipo214_500_Intropes->Draw("same");
     hC11->Draw("same");
     hC14->Draw("same");
     hCno_Numu->Draw("same");
     hK40->Draw("same");
     hK40_Intropes->Draw("same");
     hKr85->Draw("same");
     hPa234m->Draw("same");
     hPb210->Draw("same");
     hPb212->Draw("same");
     hPb214->Draw("same");
     // hPep_Numu->Draw("same");
     hPmt_Betagammas->Draw("same");
     hPo210->Draw("same");
     hPo210_Av_Innerdust->Draw("same");
     hPo210_Av_Outerdust->Draw("same");
     hPo212->Draw("same");
     hPo214->Draw("same");
     hPo216->Draw("same");
     hPo218->Draw("same");
     hRa224->Draw("same");
     hRa226->Draw("same");
     hRa228->Draw("same");
     hRn220->Draw("same");
     hRn222->Draw("same");
     hTh228->Draw("same");
     hTh230->Draw("same");
     hTh232->Draw("same");
     hTh234->Draw("same");
     hTl208->Draw("same");
     hTl208_Av->SetLineColor(kRed);
    hTl208_Av->SetLineStyle(5);
    hTl208_Av->Draw("same");
        hTl208_Av_Innerdust->Draw("same");
     hTl208_Av_Outerdust->Draw("same");
     */hTl208_Exwater->SetLineColor(kBlue);
    hTl208_Exwater->SetLineStyle(5);
    hTl208_Exwater->Draw("same");
    hnewTl->SetLineColor(kGreen);
    hnewTl->SetLineWidth(3);
    hnewTl->Draw("same");
    hnew2Tl->SetLineColor(kCyan);
    hnew2Tl->SetLineWidth(3);
    hnew2Tl->Draw("same");
    hnew3Tl->SetLineColor(kOrange);
    hnew3Tl->SetLineWidth(3);
    hnew3Tl->Draw("same");
    
    cout<<"production = "<<hTl208_Exwater->Integral()<<endl;
    cout<<"mine 1 = "<<hnewTl->Integral()<<endl;
    cout<<"mine 2 = "<<hnew2Tl->Integral()<<endl;
    cout<<"mine 3 = "<<hnew3Tl->Integral()<<endl;
    /* hTl208_Hdropes->Draw("same");
     hTl208_Hupropes->Draw("same");
     hTl208_Intropes->SetLineColor(kOrange);
    hTl208_Intropes->SetLineStyle(5);
    hTl208_Intropes->Draw("same");
        hTl210->Draw("same");
     hU234->Draw("same");
     hU238->Draw("same");*/
    /*
     TLegend *leg = new TLegend(0.6,0.4,0.9,0.9);
     leg->AddEntry(hThChain,"Th chain","l");
     leg->AddEntry(hAr39,"^{39}Ar","l");
     leg->AddEntry(hBi210, "^{210}Bi","l");
     leg->AddEntry(hUChain, "U chain","l");
     leg->AddEntry(hC11, "^{11}C","l");
     leg->AddEntry(hKr85, "^{85}Kr","l");
     leg->AddEntry(hPo210, "^{210}Po","l");
     leg->SetBorderSize(0);
     leg->Draw();
     
     TLegend *legsig = new TLegend(0.4,0.6,0.6,0.9);
     legsig->AddEntry(hPep_Nue,"pep","l");
     legsig->AddEntry(hB8_Nue,"B8","l");
     legsig->AddEntry(hBe7_Nue,"Be7","l");
     legsig->AddEntry(hCno_Nue,"CNO","l");
     legsig->SetBorderSize(0);
     legsig->Draw();
     
     TLegend *legoth = new TLegend(0.2,0.6,0.4,0.9);
     legoth->AddEntry(hTotal,"Total spectrum","l");
     legoth->AddEntry(hBackExt,"hBackExt","l");
     legoth->AddEntry(hTotal,"Total","l");
     legoth->AddEntry(hExWatermine,"Biwater","l");
     
     legoth->SetBorderSize(0);
     legoth->Draw();*/
    
    TCanvas *c2 = new TCanvas("c2");
    TH1F *hratio1 = new TH1F("hratio1",";nHits;mine/production",2000,0,2000);
    hratio1->Divide(hnewTl,hTl208_Exwater);
    hratio1->Draw();
    
    TCanvas *c3 = new TCanvas("c3");
    TH1F *hratio2 = new TH1F("hratio2",";nHits;mine2/production",2000,0,2000);
    hratio2->Divide(hnew2Tl,hTl208_Exwater);
    hratio2->Draw();
    
    TCanvas *c4 = new TCanvas("c4");
    TH1F *hratio3 = new TH1F("hratio3",";nHits;mine3/production",2000,0,2000);
    hratio3->Divide(hnew3Tl,hTl208_Exwater);
    hratio3->Draw();
    
    TCanvas *c5 = new TCanvas("c5");
    hBi214_Exwater_r->Draw();
}
