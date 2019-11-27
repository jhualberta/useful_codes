//////////////////////////////////////////////////////////////////////
///
/// \class
///
/// \brief
///
/// \author Kalpana Singh <kalpana.singh@ualberta.ca>
///
////
////////////////////////////////////////////////////////////////////////

#include <RAT/MultiPathFitter.hh>
#include "TVector3.h"
#include <CLHEP/Random/RandFlat.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <math.h>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>

namespace RAT {
  //number of steps to do after reaching tolerance
  const int MultiPathFitter::sStepsBelowTolerance = 4;
  //  convergence parameter
  int MultiPathFitter::fMaximumIterations=500 ;
  double MultiPathFitter::sTolerance = 0.001 ; // chisquare tolerance
  TVector3 MultiPathFitter::sVertex; // "best vertex"
  double MultiPathFitter::sTime0;  //best fit time;

  std::vector<double*> MultiPathFitter::sTrialparArray ;
  std::vector<double> MultiPathFitter::sBestparArray ;
  std::vector<double*> MultiPathFitter::sStartparArray ;
  std::vector<double> MultiPathFitter::sparArray ;

  TVector3 MultiPathFitter::sStartVertex; //starting vertex
  double MultiPathFitter::sStartTime0; // seed time
  int MultiPathFitter::nParameter ;
  std::list<MultiPathFitter *> MultiPathFitter::sList;
  std::vector<TH1F> MultiPathFitter::vHists;

  TMatrixT<double> MultiPathFitter::sCovariance( nParameter, nParameter);
  TMatrixT<double> MultiPathFitter::sBeta( nParameter, 1 );
  TMatrixT<double> MultiPathFitter::sAlpha(nParameter,nParameter); //curvature matrix

  //  std::vector<double> MultiPathFitter::fDerivative;
  int MultiPathFitter::sEvent = 0;
  TFile *MultiPathFitter::sFile;

  TH1F MultiPathFitter::LLHFitPosition = TH1F("LLHFitPosition","LLHFitPosition",10000,0,10000);
  TH1F MultiPathFitter::LLHoriginalPosition = TH1F("LLHoriginalPosition","LLHoriginalPosition",10000,0,10000);

  TH1F MultiPathFitter::TresFitPosition = TH1F("TresFitPosition","TresFitPosition",10000,0,10000);
  TH1F MultiPathFitter::TresOriginalPosition = TH1F("TresOriginalPosition","TresOriginalPosition",10000,0,10000);

  double MultiPathFitter::sBestLikelihood; //likelihood of best fit
  double MultiPathFitter::sLikelihood ;//keep
  unsigned int MultiPathFitter::sNHits = 0;
  unsigned int MultiPathFitter::sNHist = 0;
  double MultiPathFitter::fSpeedOfLightWater;


  MultiPathFitter::MultiPathFitter( double aTime, const TVector3 &aPosition, const int aPMTID,  int nPar ): fPosition(aPosition)
  {
    fT = aTime;
    sList.push_front( this );
    if(sList.size()==1){
      nParameter = nPar ;
      sTrialparArray.resize(nParameter) ;
      sBestparArray.resize(nParameter) ;
      sStartparArray.resize(nParameter) ;
      sparArray.reserve(nParameter) ;
      sBeta.ResizeTo(nParameter,1);
      sAlpha.ResizeTo(nParameter,nParameter);
      sCovariance.ResizeTo( nParameter, nParameter);
    }
    fPMTID = aPMTID;
  }

  void MultiPathFitter::DeleteHistograms(){
    vHists.erase(vHists.begin(),vHists.end());
  }


  void MultiPathFitter::WriteHistograms(){
    sFile->cd();
    for(size_t i = 0; i<vHists.size();i++){
      vHists[i].Write();
    }
    LLHoriginalPosition.Write(); TresOriginalPosition.Write();
    LLHFitPosition.Write(); TresFitPosition.Write();
    sFile->Close();
  }


  void MultiPathFitter::SumLikelihoods(){
    for(int j = 0;j < nParameter; j++ ){
      for(int k = 0; k <= j; k++ ) sCovariance[j][k] = 0;
      sBeta[j][0] = 0;
    }
    sLikelihood = 0;
    std::list<MultiPathFitter *>::iterator it;

    for(it = sList.begin(); it != sList.end(); it++ ){
      (*it)->CalculatePMTLikelihood();
    }
    for(int j=1;j<nParameter;j++)for(int k=0;k<j;k++)sCovariance[k][j]=sCovariance[j][k];
  }


  void MultiPathFitter::SetTrialParameters(){
    for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])= (sparArray[i]) ;
  }

  void MultiPathFitter::SetParameters(){
    for(int ii=0 ; ii< nParameter; ii++) (sparArray[ii]) = *(sTrialparArray[ii]) ;
  }


  void MultiPathFitter::ShiftParameters( TMatrixT<double> &aColumn ){
    for(int jj=0; jj< nParameter; jj++ ) *(sTrialparArray[jj]) += aColumn[jj][0] ;
  }

  void MultiPathFitter::Fit(size_t mcevent, Int_t gtID, Int_t check )
  {
    TMatrixT<double> oneda(nParameter,1) ;
    int sGoodFits = 0;
    sBestLikelihood = -1e5;
  
    //fit N times, starting at different locations within the vessel
    for(int sStart = 0 ;sGoodFits < 6 && sStart < 250; sStart++ ){
      if(sList.size()==0) continue ;
      int iStepsBelowTolerance = 0;
      double alambda = 0.001;
      bool oneStep = false; // set to true the first time we find a real step
       std::list<MultiPathFitter *>::iterator it;
      it = sList.begin() ;

      (*it)->SetInitParameters(mcevent);
      for(int i=0; i< nParameter; i++ ) {
	sparArray[i]= *(sStartparArray[i]) ;
      }

      MultiPathFitter::SetTrialParameters();
   
      SumLikelihoods();
  
      sAlpha=sCovariance;
      oneda=sBeta;
      double  oldLikelihood = sLikelihood;

      for(int iterations = 0;iterations < fMaximumIterations; iterations++) {
        //has the fit found the same likelihood value for the set amount of iterations, set alamda to 0
        if( iStepsBelowTolerance == sStepsBelowTolerance ) alambda = 0.0; //last past
        sCovariance = sAlpha;

	//check to see if matrix is singular if so do not continue
        if (sCovariance.Determinant()==0){break;}

        //augment diagonal elements to implement Marquardt
        for(int j = 0; j < nParameter; j++ ) sCovariance[j][j] = (sCovariance[j][j])*(1+alambda);
        sCovariance.InvertFast();

        TMatrixT<double> dA ;
        dA.ResizeTo(nParameter, 1) ;
        dA= sCovariance*sBeta;
        TMatrixT<double> da;
        da.ResizeTo(nParameter, 1);

        for(int j = 0; j < nParameter; j++ ) da[j][0] = dA[j][0];
        //has the fit found the same likelihood value for the set amount of iterations or alamba
        //greater than 1e90 don't proceed further else
        if( iStepsBelowTolerance == sStepsBelowTolerance || alambda > 1.0e90){
            double tmax = 0.;
            for(int j = 0; j < nParameter; j++ ) *(sTrialparArray[j])= sparArray[j] ;

            double maxLikelihood;
            maxLikelihood = oldLikelihood = sLikelihood;
            if(maxLikelihood == oldLikelihood){
              sLikelihood = oldLikelihood;
              break;
            }
            iStepsBelowTolerance = 0;
            alambda = 0.001;
            oneStep = false;
            sparArray[3] += tmax ;
            SetTrialParameters();
            SumLikelihoods();

            // sAlpha always contains last "good" curvature matrix
            for(int j = 0; j < nParameter;j++ )for(int k = 0; k < nParameter; k++ ) sAlpha[j][k] = sCovariance[j][k];

            for(int j = 0; j < nParameter; j++ ) oneda[j][0] = sBeta[j][0];  //save Beta
            oldLikelihood = sLikelihood;
        } else {//if likelihood still not minimised proceed
            ShiftParameters( dA );
            SumLikelihoods();

            //difference between previous and this likelihood should be less than the difined tolerence to
            //proceed
            if( oneStep && fabs( sLikelihood-oldLikelihood )<sTolerance )  iStepsBelowTolerance++;
            //Is this likelihood greater than the previous likelihood? If so, alamda will be multiplied by 0.1
            if( sLikelihood > oldLikelihood ){
                oneStep = true;
                alambda *= 0.1;
                oldLikelihood = sLikelihood;
                SetParameters();
                for(int j = 0; j < nParameter;j++ )for(int k = 0; k < nParameter; k++ ) sAlpha[j][k] = sCovariance[j][k];
                for(int j = 0; j < nParameter; j++ )oneda[j][0] = sBeta[j][0];
            } else {//if not alamda will be multiplied by 10
                alambda *= 10;
                sLikelihood = oldLikelihood;
                for(int j = 0; j < nParameter; j++ ) sBeta[j][0] = oneda[j][0];
                SetTrialParameters();
            }
        }
      } //for iterations

      sVertex.SetXYZ(sparArray[0], sparArray[1],sparArray[2]) ;
      if(sLikelihood > sBestLikelihood-0.5  &&sVertex.Mag()<8400){  //good fit
        if(sBestLikelihood < sLikelihood-0.25)sGoodFits = 0;// resets sGoodfits
        else sGoodFits++;//increments sGoodFits (need 6 good fits to have a valid fit)
        if( sLikelihood > sBestLikelihood ){//reset sBestLikelihood to current likelihood value
          for(int jk=0; jk<nParameter ; jk++) (sBestparArray[jk]) = sparArray[jk] ;
          sBestLikelihood = sLikelihood;
        }
      }
    }  //Start Position iterations
    sEvent++;
    
    if(sGoodFits == 0 || sList.size()==0){//not valid fit set nonphysical values
      for(int jk=0; jk<nParameter ; jk++) (sBestparArray[jk]) = -10000 ;
    }
    if(sList.size()!=0){
    std::list<MultiPathFitter *>::iterator it;
    for( it=sList.begin();it != sList.end(); it++ ){
      (*it)->CalculatePMTOther();
     }
    
    it--;
    sNHits = sList.size();
    (*it)->SummariseOtherOutputs();
    }
    for(int jk=0; jk<nParameter ; jk++) *(sTrialparArray[jk]) = sBestparArray[jk] ;

    if((gtID==1287420|| gtID==1288364|| gtID==1290696|| gtID==1296100|| gtID==1300157|| gtID==1301165|| gtID==1302927|| gtID==1303599|| gtID==1303777|| gtID==1318732|| gtID==1335548|| gtID==1344824|| gtID==1345531|| gtID==1346677)  && check==1) { 
//{ // GTID Loop to dump the likelihoods for event ID 64

    TFile *sFile=new TFile(Form("dumpfit_%d.root",gtID),"RECREATE");

    TH1F *sLLHFitPosition =new TH1F("sLLHFitPosition","LLHFitPosition",10000,0,10000);
    TH1F *sLLHoriginalPosition =new TH1F("sLLHoriginalPosition","LLHoriginalPosition",10000,0,10000);
    TH1F *TresFitPosition =new TH1F("TresFitPosition","TresFitPosition",10000,0,10000);
    TH1F *TresOriginalPosition =new TH1F("TresOriginalPosition","TresOriginalPosition",10000,0,10000);
    TH2F *pmtGEOLLH = new TH2F("pmtGEOLLH","PMT Geometry view of LLH",360,-180,180,180,0,180 );

    TH1F *AngleTres = new TH1F("AngleTres","angle of incidence vs. tres.",90,0,90);

    TH1F *sHx=new TH1F("Hx","Hx",1000,-5000,5000);
    TH1F *sHy=new TH1F("Hy","Hy",1000,-5000,5000);
    TH1F *sHz=new TH1F("Hz","Hz",1000,-5000,5000);
    TH1F *sHt=new TH1F("Ht","Ht",600,-300,300);

    TH1F *sDHx=new TH1F("sDHx","sDHx",1000,-5000,5000);
    TH1F *sDHy=new TH1F("sDHy","sDHy",1000,-5000,5000);
    TH1F *sDHz=new TH1F("sDHz","sDHz",1000,-5000,5000);
    TH1F *sDHt=new TH1F("sDHt","sDHt",600,-300,300);
   
    TH1F *sNDHx=new TH1F("sNDHx","sNDHx",1000,-5000,5000);
    TH1F *sNDHy=new TH1F("sNDHy","sNDHy",1000,-5000,5000);
    TH1F *sNDHz=new TH1F("sNDHz","sNDHz",1000,-5000,5000);
    TH1F *sNDHt=new TH1F("sNDHt","sNDHt",600,-300,300);

    TVector3 AVPosition, fitPosition, pmtUnitPosition ;

    double oldLikelihood=0.0, x=0.0, y=0.0, z=0.0, t=0.0, angle=0.0 ;

    int count =0 ;
    for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
  
    fitPosition.SetXYZ(sBestparArray[0],sBestparArray[1],sBestparArray[2]) ;
    
    std::list<MultiPathFitter *>::iterator it;
    it = sList.begin() ;

    for(it = sList.begin(); it != sList.end(); it++ ){
       sLikelihood= 0.0;
        (*it)->CalculatePMTLikelihood();
        double tres = (*it)->GetResidual() ;
	int PMTid = (*it)->GetPMTID() ;
        pmtUnitPosition = (*it)->GetPosition() ;
        AVPosition= 6.0*pmtUnitPosition.Unit() ;
        angle = TMath::ACos((fitPosition-AVPosition).Unit()*AVPosition.Unit())*180/TMath::Pi() ;
	sLLHFitPosition->Fill(PMTid, sLikelihood);
        TresFitPosition->Fill(PMTid, tres);
        AngleTres->Fill(angle,tres) ;
        pmtGEOLLH->Fill(pmtUnitPosition.Phi()*180/TMath::Pi(),pmtUnitPosition.Theta()*180/TMath::Pi(),sLikelihood) ;
     }

    for(int i=0; i< nParameter; i++ ) { 
    if(i==0) *(sTrialparArray[i])= -0.25 ;
    if(i==1) *(sTrialparArray[i])= 0.15 ;
    if(i==2) *(sTrialparArray[i])= 70.99 ;
    if(i==3) *(sTrialparArray[i])= *(sStartparArray[i]) ;
    }

    for(it = sList.begin(); it != sList.end(); it++ ){
       sLikelihood= 0.0;
        (*it)->CalculatePMTLikelihood();
        double tres = (*it)->GetResidual() ;
	int PMTid = (*it)->GetPMTID() ;
	TresOriginalPosition->Fill(PMTid, tres);
	sLLHoriginalPosition->Fill(PMTid, sLikelihood);
    }

    for(int i=0; i< nParameter; i++ ){ *(sTrialparArray[i])=sBestparArray[i] ;
    std::cout<< " best parameters "<<i<< " is "<<sBestparArray[i] <<"\t" ;}
    SetParameters();
    oldLikelihood=0.0 ;
    count =0 ;
    sHx->Reset(); sDHx-> Reset(); sNDHx-> Reset(); 
    for(x=-5000;x<5000;x+=10.0){        
     *(sTrialparArray[0])=sparArray[0]+x;
      SumLikelihoods();
      if(count ==0){
         oldLikelihood = sLikelihood ;
         count++ ;
      }
      sHx->Fill(sparArray[0]+x, sLikelihood);
      sDHx->Fill(sparArray[0]+x,sBeta[0][0]);
      if(count!=0) {
	  sNDHx->Fill(sparArray[0]+x, (sLikelihood-oldLikelihood)/10.0) ;
	  oldLikelihood = sLikelihood ;
	}	
     } // end of the x coordinate loop
    
     for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
     SetParameters();
     oldLikelihood=0.0 ;
     count =0 ;

     SetParameters();
     sHy->Reset(); sDHy-> Reset(); sNDHy-> Reset();
     for(y=-5000;y<5000;y+=10.0){
	*(sTrialparArray[1])=sparArray[1]+y;
	SumLikelihoods();
        if(count ==0){
            oldLikelihood = sLikelihood ;
            count++ ;
          }

	sHy->Fill(sparArray[1]+y, sLikelihood);
        sDHy->Fill(sparArray[1]+y,sBeta[1][0]);
        if(count!=0) {
	  sNDHy->Fill(sparArray[1]+y, (sLikelihood-oldLikelihood)/10.0) ;
                oldLikelihood = sLikelihood ;
          }

     } // end of the y coordinate loop

     for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
     SetParameters();
     oldLikelihood=0.0 ;
     count =0 ;

    sHz->Reset(); sDHz-> Reset(); sNDHz-> Reset();
    for(z=-5000;z<5000;z+=10.0){
	*(sTrialparArray[2])=sparArray[2]+z;
	SumLikelihoods();
        if(count ==0){
            oldLikelihood = sLikelihood ;
            count++ ;
          }

	sHz->Fill(sparArray[2]+z, sLikelihood);
        sDHz->Fill(sparArray[2]+z, sBeta[2][0]) ;
        if(count!=0) {
	  sNDHz->Fill(sparArray[2]+z, (sLikelihood-oldLikelihood)/10.0) ;
                oldLikelihood = sLikelihood ;
          }

     } // end of the z coordinate loop

     for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
     SetParameters();
     oldLikelihood=0.0 ;
     count =0 ;

     sHt->Reset(); sDHt-> Reset(); sNDHt-> Reset();
     for(t=-300;t<300;t+=1.0){
	*(sTrialparArray[3])=sparArray[3]+t;
	SumLikelihoods();
        if(count ==0){
            oldLikelihood = sLikelihood ;
            count++ ;
          }
	sHt->Fill(sparArray[3]+t, sLikelihood);
        sDHt->Fill(sparArray[3]+t,sBeta[3][0]);
        if(count!=0) {
                sNDHt->Fill(sparArray[3]+t, sLikelihood-oldLikelihood) ;
                oldLikelihood = sLikelihood ;
          }

    } // end of the time loop
    sFile->cd();
    sLLHFitPosition->Write(); sLLHoriginalPosition->Write(); TresFitPosition->Write(); TresOriginalPosition->Write(); sHx->Write(); sHy->Write(); sHz->Write(); sHt->Write(); sDHx->Write(); sDHy->Write(); sDHz->Write(); sDHt->Write(); sNDHx->Write();  sNDHy->Write(); sNDHz->Write(); sNDHt->Write(); AngleTres->Write(); pmtGEOLLH->Write();
    sFile->Close();
    } // end of the GTID loop for likelihood dump
  }


  void MultiPathFitter::DumpLikelihood(int sNumberOfEvents){
    //writes the likelihood surface and derivative only if asked for
    if(sEvent>sNumberOfEvents)return;
    std::list<MultiPathFitter *>::iterator itn;
    //loop to make histograms needed to fill
    for( itn = sList.begin(); itn != sList.end(); itn++ ){
      (*itn)->MakeParameterHistogram();
      break;
    }

    //how many different histograms there are
    if(sEvent==1){sNHist = vHists.size();}
    sFile = new TFile("dumpfit.root","recreate");
    for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
    SetParameters();
    SumLikelihoods();

    double oldLikelihood = 0.0 ;
    int count = 0 ;
    std::list<MultiPathFitter *>::iterator it;
    it = sList.begin() ;
    MultiPathFitter::SetTrialParameters();

    TVector3 PMTPosition ;
    oldLikelihood = 0.0 ;
    count = 0 ;
    for(it = sList.begin(); it != sList.end(); it++ ){
      sLikelihood = 0.0;
      (*it)->CalculatePMTLikelihood();
      double tres = (*it)->GetResidual() ;
      int PMTid = (*it)->GetPMTID() ;
      LLHFitPosition.Fill(PMTid, sLikelihood);
      TresFitPosition.Fill(PMTid, tres);
    }

    for(int i=0; i< 3; i++ ) *(sTrialparArray[i])=0.0 ;
    for(it = sList.begin(); it != sList.end(); it++ ){
      sLikelihood= 0.0;
      (*it)->CalculatePMTLikelihood();
      double tres = (*it)->GetResidual() ;
      int PMTid = (*it)->GetPMTID() ;
      TresOriginalPosition.Fill(PMTid, tres);
      LLHoriginalPosition.Fill(PMTid, sLikelihood);
    }
 
    for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
    for(int j = 0 ; j < nParameter;j++){
      oldLikelihood=0.0 ;
      count =0 ;
      vHists[j].Reset();
      //SetParameters();
      double lowBin = vHists[j].GetBinLowEdge(1);
      double highBin = vHists[j].GetBinLowEdge(vHists[j].GetNbinsX()+1);
      double increment = (vHists[j].GetBinLowEdge(1) - vHists[j].GetBinLowEdge(0));
      int NIncreament = (int)(highBin-lowBin)/increment ;
     
       for(int k = 0; k<NIncreament;k++){
	//std::cout<<" parameter "<< j+1<<" binID "<<k<<" parameter Value "<<lowBin+k*increment<<std::endl ;
        *(sTrialparArray[j])=lowBin+k*increment;
        SumLikelihoods();
        if(count ==0){
          oldLikelihood = sLikelihood ;
          count++ ;
        }
        unsigned int iLikelihood = vHists.size()-sNHist+j;
        unsigned int iDerivative =vHists.size()-sNHist+nParameter+j;
        unsigned int iNDerivative =vHists.size()-sNHist+(nParameter*2)+j;

        vHists[iLikelihood].Fill(lowBin+k*increment, sLikelihood);
        vHists[iDerivative].Fill(lowBin+k*increment, sBeta[j][0]);
        if(count!=0) {
            double w = sLikelihood-oldLikelihood;
            double nDerivative = w;
            if (w !=0.0 ){
              nDerivative = w/10.0;
            }
            vHists[iNDerivative].Fill(lowBin+k*increment, nDerivative) ;
            oldLikelihood = sLikelihood ;
          }
	//std::cout<<" bin number "<< k<< ", ";
        }
        for(int i=0; i< nParameter; i++ ) *(sTrialparArray[i])=sBestparArray[i] ;
      }
    // write file for last event
    if(sEvent==sNumberOfEvents){
      WriteHistograms();
      DeleteHistograms();
    }
  }
} /* namespace RAT */

