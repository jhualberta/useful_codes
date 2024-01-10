//CPP
#include <iostream>
//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
//RAT
#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/TrackNav.hh>
#include <RAT/TrackCursor.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <RAT/DS/UniversalTime.hh>
//定义想要记录的数据结构
struct Result
{
    void ClearMeta()
    {
        C_MCID = 0;
        C_RunID = 0;
        C_SubRunID = 0;
    };
    void ClearParent()
    {
        //C_Parent_Time = 0;
        C_Parent_PosX = 0;
        C_Parent_PosY = 0;
        C_Parent_PosZ = 0;

        C_Parent_KineticE = 0;
        C_Parent_MomX = 0;
        C_Parent_MomY = 0;
        C_Parent_MomZ = 0;

        C_Parent_PDG = 0;
        C_Parent_VertexID = 0;
    };
    void ClearParticle1()
    {
        C_Particle1_Time = 0;
        C_Particle1_PosX = 0;
        C_Particle1_PosY = 0;
        C_Particle1_PosZ = 0;

        C_Particle1_KineticE = 0;
        C_Particle1_MomX = 0;
        C_Particle1_MomY = 0;
        C_Particle1_MomZ = 0;

        C_Particle1_PDG = 0;

        C_Particle1_TrackID = 0;
        C_Particle1_StepID = 0;
        C_Particle1_IsTrackStart = 0;
        C_Particle1_IsTrackEnd = 0;



        C_Particle1_An_Time = 0;
        C_Particle1_An_PosX = 0;
        C_Particle1_An_PosY = 0;
        C_Particle1_An_PosZ = 0;

        C_Particle1_An_KineticE = 0;
        C_Particle1_An_MomX = 0;
        C_Particle1_An_MomY = 0;
        C_Particle1_An_MomZ = 0;

        C_Particle1_An_PDG = 0;
        
        C_Particle1_An_TrackID = 0;
        C_Particle1_An_StepID = 0;
        C_Particle1_An_IsTrackStart = 0;
        C_Particle1_An_IsTrackEnd = 0;
    };
    void ClearParticle2()
    {
        C_Particle2_Time = 0;
        C_Particle2_PosX = 0;
        C_Particle2_PosY = 0;
        C_Particle2_PosZ = 0;

        C_Particle2_KineticE = 0;
        C_Particle2_MomX = 0;
        C_Particle2_MomY = 0;
        C_Particle2_MomZ = 0;

        C_Particle2_PDG = 0;

        C_Particle2_TrackID = 0;
        C_Particle2_StepID = 0;
        C_Particle2_IsTrackStart = 0;
        C_Particle2_IsTrackEnd = 0;

        C_Particle2_pCap_Time = 0;
        C_Particle2_pCap_PosX = 0;
        C_Particle2_pCap_PosY = 0;
        C_Particle2_pCap_PosZ = 0;

        C_Particle2_pCap_KineticE = 0;
        C_Particle2_pCap_MomX = 0;
        C_Particle2_pCap_MomY = 0;
        C_Particle2_pCap_MomZ = 0;

        C_Particle2_pCap_PDG = 0;

        C_Particle2_pCap_TrackID = 0;
        C_Particle2_pCap_StepID = 0;
        C_Particle2_pCap_IsTrackStart = 0;
        C_Particle2_pCap_IsTrackEnd = 0;
    };
//Run等信息
    Int_t C_MCID;
    UInt_t C_RunID, C_SubRunID;
//Parent信息
    //RAT::UniversalTime C_parent_Time;
    Double_t C_Parent_PosX, C_Parent_PosY, C_Parent_PosZ;
    Double_t C_Parent_KineticE, C_Parent_MomX, C_Parent_MomY, C_Parent_MomZ;
    Int_t C_Parent_PDG, C_Parent_VertexID;
//Particle1------e+信息
    //初始信息
    Double_t C_Particle1_Time, C_Particle1_PosX, C_Particle1_PosY, C_Particle1_PosZ;
    Double_t C_Particle1_KineticE, C_Particle1_MomX, C_Particle1_MomY, C_Particle1_MomZ;
    Int_t C_Particle1_PDG;
    int C_Particle1_TrackID, C_Particle1_StepID, C_Particle1_IsTrackStart, C_Particle1_IsTrackEnd;
    //湮灭信息
    Double_t C_Particle1_An_Time, C_Particle1_An_PosX, C_Particle1_An_PosY, C_Particle1_An_PosZ;
    Double_t C_Particle1_An_KineticE, C_Particle1_An_MomX, C_Particle1_An_MomY, C_Particle1_An_MomZ;
    Int_t C_Particle1_An_PDG;
    int C_Particle1_An_TrackID, C_Particle1_An_StepID, C_Particle1_An_IsTrackStart, C_Particle1_An_IsTrackEnd;
    
//Particle2--------n信息
    //初始信息
    Double_t C_Particle2_Time, C_Particle2_PosX, C_Particle2_PosY, C_Particle2_PosZ;
    Double_t C_Particle2_KineticE, C_Particle2_MomX, C_Particle2_MomY, C_Particle2_MomZ;
    Int_t C_Particle2_PDG;
    int C_Particle2_TrackID, C_Particle2_StepID, C_Particle2_IsTrackStart, C_Particle2_IsTrackEnd;
    //湮灭信息
    Double_t C_Particle2_pCap_Time, C_Particle2_pCap_PosX, C_Particle2_pCap_PosY, C_Particle2_pCap_PosZ;
    Double_t C_Particle2_pCap_KineticE, C_Particle2_pCap_MomX, C_Particle2_pCap_MomY, C_Particle2_pCap_MomZ;
    Int_t C_Particle2_pCap_PDG;
    int C_Particle2_pCap_TrackID, C_Particle2_pCap_StepID, C_Particle2_pCap_IsTrackStart, C_Particle2_pCap_IsTrackEnd;
};

int main(int argc, char** argv)
{
    // std::string InFile = "/home/shuaioy/scratch/Geo/Reactor/XDireaction/ratds/temp/ScintFit_2p2ReactoribdRun_r306499_s0_p0.root";
    std::string InFile = argv[1];
    std::string OutFile = argv[2];
    // std::string InFile = "/home/shuaioy/scratch/Geo/Reactor/XDireaction/ratds/ROOT/306499.root";
    // std::string OutFile = "/home/shuaioy/scratch/Geo/Reactor/XDireaction/ratds/ROOT/ReactorIBD_OutPut_306499.root";

    TFile *outfile = new TFile(OutFile.c_str(), "recreate");
    TTree *tree = new TTree("output", "output");

    Result res;
//Meta
    tree->Branch("MCID", &res.C_MCID, "MCID/I");
    tree->Branch("RunID", &res.C_RunID, "RunID/I");
    tree->Branch("SubRunID", &res.C_SubRunID, "SubRunID/I");
//Parent
    tree->Branch("ParentPosX", &res.C_Parent_PosX, "ParentPosX/D");
    tree->Branch("ParentPosY", &res.C_Parent_PosY, "ParentPosY/D");
    tree->Branch("ParentPosZ", &res.C_Parent_PosZ, "ParentPosZ/D");
    tree->Branch("ParentKineticE", &res.C_Parent_KineticE, "ParentKineticE/D");
    tree->Branch("ParentMomX", &res.C_Parent_MomX, "ParentMomX/D");
    tree->Branch("ParentMomY", &res.C_Parent_MomY, "ParentMomY/D");
    tree->Branch("ParentMomZ", &res.C_Parent_MomZ, "ParentMomZ/D");
    tree->Branch("ParentPDG", &res.C_Parent_PDG, "ParentPDG/I");
    tree->Branch("ParentVertexID", &res.C_Parent_VertexID, "ParentVertexID/I");
//Particle1
    //初始信息
    tree->Branch("Particle1Time", &res.C_Particle1_Time, "Particle1Time/D");
    tree->Branch("Particle1PosX", &res.C_Particle1_PosX, "Particle1PosX/D");
    tree->Branch("Particle1PosY", &res.C_Particle1_PosY, "Particle1PosY/D");
    tree->Branch("Particle1PosZ", &res.C_Particle1_PosZ, "Particle1PosZ/D");
    tree->Branch("Particle1KineticE", &res.C_Particle1_KineticE, "Particle1KineticE/D");
    tree->Branch("Particle1MomX", &res.C_Particle1_MomX, "Particle1MomX/D");
    tree->Branch("Particle1MomY", &res.C_Particle1_MomY, "Particle1MomY/D");
    tree->Branch("Particle1MomZ", &res.C_Particle1_MomZ, "Particle1MomZ/D");
    tree->Branch("Particle1PDG", &res.C_Particle1_PDG, "Particle1PDG/I");
    tree->Branch("Particle1TrackID", &res.C_Particle1_TrackID, "Particle1TrackID/I");
    tree->Branch("Particle1StepID", &res.C_Particle1_StepID, "Particle1StepID/I");
    tree->Branch("Particle1IsTrackStart", &res.C_Particle1_IsTrackStart, "Particle1IsTrackStart/I");
    tree->Branch("Particle1IsTrackEnd", &res.C_Particle1_IsTrackEnd, "Particle1IsTrackEnd/I");
    //湮灭信息
    tree->Branch("Particle1_An_Time", &res.C_Particle1_An_Time, "Particle1_An_Time/D");
    tree->Branch("Particle1_An_PosX", &res.C_Particle1_An_PosX, "Particle1_An_PosX/D");
    tree->Branch("Particle1_An_PosY", &res.C_Particle1_An_PosY, "Particle1_An_PosY/D");
    tree->Branch("Particle1_An_PosZ", &res.C_Particle1_An_PosZ, "Particle1_An_PosZ/D");
    tree->Branch("Particle1_An_KineticE", &res.C_Particle1_An_KineticE, "Particle1_An_KineticE/D");
    tree->Branch("Particle1_An_MomX", &res.C_Particle1_An_MomX, "Particle1_An_MomX/D");
    tree->Branch("Particle1_An_MomY", &res.C_Particle1_An_MomY, "Particle1_An_MomY/D");
    tree->Branch("Particle1_An_MomZ", &res.C_Particle1_An_MomZ, "Particle1_An_MomZ/D");
    tree->Branch("Particle1_An_PDG", &res.C_Particle1_An_PDG, "Particle1_An_PDG/I");
    tree->Branch("Particle1_An_TrackID", &res.C_Particle1_An_TrackID, "Particle1_An_TrackID/I");
    tree->Branch("Particle1_An_StepID", &res.C_Particle1_An_StepID, "Particle1_An_StepID/I");
    tree->Branch("Particle1_An_IsTrackStart", &res.C_Particle1_An_IsTrackStart, "Particle1_An_IsTrackStart/I");
    tree->Branch("Particle1_An_IsTrackEnd", &res.C_Particle1_An_IsTrackEnd, "Particle1_An_IsTrackEnd/I");
//Particle2
    //中子初始信息
    tree->Branch("Particle2Time", &res.C_Particle2_Time, "Particle2Time/D");
    tree->Branch("Particle2PosX", &res.C_Particle2_PosX, "Particle2PosX/D");
    tree->Branch("Particle2PosY", &res.C_Particle2_PosY, "Particle2PosY/D");
    tree->Branch("Particle2PosZ", &res.C_Particle2_PosZ, "Particle2PosZ/D");
    tree->Branch("Particle2KineticE", &res.C_Particle2_KineticE, "Particle2KineticE/D");
    tree->Branch("Particle2MomX", &res.C_Particle2_MomX, "Particle2MomX/D");
    tree->Branch("Particle2MomY", &res.C_Particle2_MomY, "Particle2MomY/D");
    tree->Branch("Particle2MomZ", &res.C_Particle2_MomZ, "Particle2MomZ/D");
    tree->Branch("Particle2PDG", &res.C_Particle2_PDG, "Particle2PDG/I");
    tree->Branch("Particle2TrackID", &res.C_Particle2_TrackID, "Particle2TrackID/I");
    tree->Branch("Particle2StepID", &res.C_Particle2_StepID, "Particle2StepID/I");
    tree->Branch("Particle2IsTrackStart", &res.C_Particle2_IsTrackStart, "Particle2IsTrackStart/I");
    tree->Branch("Particle2IsTrackEnd", &res.C_Particle2_IsTrackEnd, "Particle2IsTrackEnd/I");
    //质子捕获信息
    tree->Branch("Particle2_pCap_Time", &res.C_Particle2_pCap_Time, "Particle2_pCap_Time/D");
    tree->Branch("Particle2_pCap_PosX", &res.C_Particle2_pCap_PosX, "Particle2_pCap_PosX/D");
    tree->Branch("Particle2_pCap_PosY", &res.C_Particle2_pCap_PosY, "Particle2_pCap_PosY/D");
    tree->Branch("Particle2_pCap_PosZ", &res.C_Particle2_pCap_PosZ, "Particle2_pCap_PosZ/D");
    tree->Branch("Particle2_pCap_KineticE", &res.C_Particle2_pCap_KineticE, "Particle2_pCap_KineticE/D");
    tree->Branch("Particle2_pCap_MomX", &res.C_Particle2_pCap_MomX, "Particle2_pCap_MomX/D");
    tree->Branch("Particle2_pCap_MomY", &res.C_Particle2_pCap_MomY, "Particle2_pCap_MomY/D");
    tree->Branch("Particle2_pCap_MomZ", &res.C_Particle2_pCap_MomZ, "Particle2_pCap_MomZ/D");
    tree->Branch("Particle2_pCap_PDG", &res.C_Particle2_pCap_PDG, "Particle2_pCap_PDG/I");
    tree->Branch("Particle2_pCap_TrackID", &res.C_Particle2_pCap_TrackID, "Particle2_pCap_TrackID/I");
    tree->Branch("Particle2_pCap_StepID", &res.C_Particle2_pCap_StepID, "Particle2_pCap_StepID/I");
    tree->Branch("Particle2_pCap_IsTrackStart", &res.C_Particle2_pCap_IsTrackStart, "Particle2_pCap_IsTrackStart/I");
    tree->Branch("Particle2_pCap_IsTrackEnd", &res.C_Particle2_pCap_IsTrackEnd, "Particle2_pCap_IsTrackEnd/I");

    RAT::DU::DSReader dsReader(InFile);
    const size_t lenEntry = dsReader.GetEntryCount();
    //提前声明，减少loop中重复创建与消除，减少内存的占用
    RAT::DS::Entry entry;
    RAT::DS::MC mc;
    RAT::DS::MCParticle parent;
    TVector3 temp;
    RAT::TrackNode *node;
    RAT::TrackNav *nav = NULL;
    RAT::TrackCursor *cursor = NULL, *particle1 = NULL, *particle2 = NULL;
//开始筛选
    std::cout << "开始筛选" << std::endl;
    for(size_t ii1 = 0; ii1 < lenEntry; ii1++)
    {   
        std::cout << "开始筛选" << ii1 << ", 还剩" << lenEntry - ii1 << std::endl;
        entry = dsReader.GetEntry(ii1);
    //MC Meta
        mc = entry.GetMC();
        //记录数据
        res.ClearMeta();
        res.C_MCID = mc.GetMCID();
        res.C_RunID = entry.GetRunID();
        res.C_SubRunID = entry.GetSubRunID();
        if(ii1 > 60){std::cout << ii1 << ", MC完成" << std::endl;};
    //中微子初始信息：
        parent = mc.GetMCParent(0);
        //记录数据
        res.ClearParent();
        //res.C_Parent_Time = parent.GetTime();
        res.C_Parent_KineticE = parent.GetKineticEnergy();
        res.C_Parent_PDG = parent.GetPDGCode();
        res.C_Parent_VertexID = parent.GetVertexID();
        temp = parent.GetPosition();
        res.C_Parent_PosX = temp.X();
        res.C_Parent_PosY = temp.Y();
        res.C_Parent_PosZ = temp.Z();
        temp = parent.GetMomentum();
        res.C_Parent_MomX = temp.X();
        res.C_Parent_MomY = temp.Y();
        res.C_Parent_MomZ = temp.Z();
        if(ii1 > 60){std::cout << ii1 << ", Parent完成" << std::endl;};
    //Track
        nav = new RAT::TrackNav(&entry);
        if(ii1 > 60){std::cout << ii1 << ", Nav完成" << std::endl;};
        cursor = new RAT::TrackCursor(nav->Cursor(false));
        if(ii1 > 60){std::cout << ii1 << ", Cursor完成" << std::endl;};
    //正电子的信息
        particle1 = new RAT::TrackCursor (cursor->Child(0));//可以直接转成RAT::TrackCurosor
        res.ClearParticle1();
        //记录数据
        //正电子产生
        node = particle1->TrackStart();
        res.C_Particle1_TrackID = node->GetTrackID();
        res.C_Particle1_StepID = node->GetStepID();
        res.C_Particle1_IsTrackStart = node->IsTrackStart();
        res.C_Particle1_IsTrackEnd = node->IsTrackEnd();

        res.C_Particle1_Time = node->GetGlobalTime();
        res.C_Particle1_KineticE = node->GetKineticEnergy();
        res.C_Particle1_PDG = node->GetPDGCode();
        temp = node->GetPosition();
        res.C_Particle1_PosX = temp.X();
        res.C_Particle1_PosY = temp.Y();
        res.C_Particle1_PosZ = temp.Z();
        temp = node->GetMomentum();
        res.C_Particle1_MomX = temp.X();
        res.C_Particle1_MomY = temp.Y();
        res.C_Particle1_MomZ = temp.Z();
        if(ii1 > 60){std::cout << ii1 << ", e+完成" << std::endl;};
        //正电子湮灭
        node = particle1->TrackEnd();
        res.C_Particle1_An_TrackID = node->GetTrackID();
        res.C_Particle1_An_StepID = node->GetStepID();
        res.C_Particle1_An_IsTrackStart = node->IsTrackStart();
        res.C_Particle1_An_IsTrackEnd = node->IsTrackEnd();

        res.C_Particle1_An_Time = node->GetGlobalTime();
        res.C_Particle1_An_KineticE = node->GetKineticEnergy();
        res.C_Particle1_An_PDG = node->GetPDGCode();
        temp = node->GetPosition();
        res.C_Particle1_An_PosX = temp.X();
        res.C_Particle1_An_PosY = temp.Y();
        res.C_Particle1_An_PosZ = temp.Z();
        temp = node->GetMomentum();
        res.C_Particle1_An_MomX = temp.X();
        res.C_Particle1_An_MomY = temp.Y();
        res.C_Particle1_An_MomZ = temp.Z();
        if(ii1 > 60){std::cout << ii1 << ", e+湮灭完成" << std::endl;};
    //中子信息
        particle2 = new RAT::TrackCursor (cursor->Child(1));
        //记录数据
        res.ClearParticle2();
        //中子产生
        node = particle2->TrackStart();
        res.C_Particle2_TrackID = node->GetTrackID();
        res.C_Particle2_StepID = node->GetStepID();
        res.C_Particle2_IsTrackStart = node->IsTrackStart();
        res.C_Particle2_IsTrackEnd = node->IsTrackEnd();

        res.C_Particle2_Time = node->GetGlobalTime();
        res.C_Particle2_KineticE = node->GetKineticEnergy();
        res.C_Particle2_PDG = node->GetPDGCode();
        temp = node->GetPosition();
        res.C_Particle2_PosX = temp.X();
        res.C_Particle2_PosY = temp.Y();
        res.C_Particle2_PosZ = temp.Z();
        temp = node->GetMomentum();
        res.C_Particle2_MomX = temp.X();
        res.C_Particle2_MomY = temp.Y();
        res.C_Particle2_MomZ = temp.Z();
        if(ii1 > 60){std::cout << ii1 << ", 中子完成" << std::endl;};
        //中子捕获
        node = particle1->TrackEnd();
        res.C_Particle2_pCap_TrackID = node->GetTrackID();
        res.C_Particle2_pCap_StepID = node->GetStepID();
        res.C_Particle2_pCap_IsTrackStart = node->IsTrackStart();
        res.C_Particle2_pCap_IsTrackEnd = node->IsTrackEnd();

        res.C_Particle2_pCap_Time = node->GetGlobalTime();
        res.C_Particle2_pCap_KineticE = node->GetKineticEnergy();
        res.C_Particle2_pCap_PDG = node->GetPDGCode();
        temp = node->GetPosition();
        res.C_Particle2_pCap_PosX = temp.X();
        res.C_Particle2_pCap_PosY = temp.Y();
        res.C_Particle2_pCap_PosZ = temp.Z();
        temp = node->GetMomentum();
        res.C_Particle2_pCap_MomX = temp.X();
        res.C_Particle2_pCap_MomY = temp.Y();
        res.C_Particle2_pCap_MomZ = temp.Z();
        if(ii1 > 60){std::cout << ii1 << ", 中子捕获完成" << std::endl;};
        tree->Fill();
        //释放内存
        delete nav; nav = NULL;
        delete cursor; cursor = NULL;
        delete particle1; particle1 = NULL;
        delete particle2; particle2 = NULL;
        if(ii1 > 60){std::cout << ii1 << ", 清除完成" << std::endl;};
    };
    std::cout << "筛选完成" << std::endl;
    outfile->Write();
    outfile->Close();
    
    delete outfile; outfile = NULL;
    delete tree; tree = NULL;

    return 0;
}