#include "StMyAnalysisMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "TLorentzVector.h"
// #include "StPicoDstMaker/StPicoV0.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoMtdPidTraits.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include <vector>
#include "phys_constants.h"
#include "StCuts.h"
#include "StBTofUtil/tofPathLength.hh"
//
// #include "TMinuit.h"
// #include "StCuts.h"
// #include "StPicoKFVertexFitter/StPicoKFVertexFitter.h"
// #include "PhysicalConstants.h"
// #include "phys_constants.h"
// #include "StBTofUtil/tofPathLength.hh"
// #include "StRoot/StRefMultCorr/StRefMultCorr.h"
//
//

//
ClassImp(StMyAnalysisMaker)
//-----------------------------------------------------------------------------
StMyAnalysisMaker::StMyAnalysisMaker(const char* name, StPicoDstMaker *picoMaker, const char* outName)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mOutName = outName;
}

//----------------------------------------------------------------------------- 
StMyAnalysisMaker::~StMyAnalysisMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Init() {
  DeclareHistograms();

  dcaG = new StDcaGeometry();
  dca = new StDcaGeometry();
  // mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Finish() {
  if(mOutName!="") {
    TFile *fout = new TFile(mOutName.Data(),"RECREATE");
    fout->cd();
    // WriteHistograms();
    tuplePair->Write();
    eventCheck->Write();
    trackCheckTof->Write();
    trackCheckNoTof->Write();
    deltaVZ->Write();
    fout->Write();
    fout->Close();
  }
  SafeDelete(dca);
  SafeDelete(dcaG);
	
  return kStOK;
}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::DeclareHistograms() {
  mD0Tuple = new TNtuple("mD0Tuple","","mass:pt:centrality");

  eventCheck = new TH1D("eventCheck","",10,-0.5,9.5);
  trackCheckTof = new TH1D("trackCheckTof","",10,-0.5,9.5);
  trackCheckNoTof = new TH1D("trackCheckNoTof","",10,-0.5,9.5);
   hPVZ = new TH2D("pvz",";tpc vz;vpd vz",200,-100,100,200,-100,100);
   hMult = new TH1D("hMult","",100,0,100);
  float multBin[6] = {0,7,12,16,22,100};
  float dzBin[2001];
  for(int i =0; i< 2001; i++)
    dzBin[i] = i*0.01-10.;
   deltaVZ = new TH2D("deltaVZ","",5,multBin,2000,dzBin);

   double hftPtBins[20] = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.,1.2,1.42,4,6,8,10};
   hEta = new TH1D("eta","",200,-1.5,1.5);
   hPhi = new TH1D("phi","",200,-4,4);
   hPt = new TH1D("pt","",19,hftPtBins);
   hPt->Sumw2();
   hDcaXY = new TH2D("dcaXY",";pt;dcaXY",200,0,5,1000,-0.5,0.5);
   hDcaZ = new TH2D("dcaZ",";pt;dcaZ",200,0,5,1000,-0.5,0.5);

   hTPCEta = new TH1D("tpc_eta","",200,-1.5,1.5);
   hTPCPhi = new TH1D("tpc_phi","",200,-4,4);
   hTPCPt = new TH1D("tpc_pt","",19,hftPtBins);
   hTPCPt->Sumw2();
   hTPCDcaXY = new TH2D("tpc_dcaXY",";pt;dcaXY",200,0,5,1000,-0.5,0.5);
   hTPCDcaZ = new TH2D("tpc_dcaZ",";pt;dcaZ",200,0,5,1000,-0.5,0.5);

   hHFTEta = new TH1D("hft_eta","",200,-1.5,1.5);
   hHFTPhi = new TH1D("hft_phi","",200,-4,4);
   hHFTPt = new TH1D("hft_pt","",19,hftPtBins);
   hHFTPt->Sumw2();
   hHFTDcaXY = new TH2D("hft_dcaXY",";pt;dcaXY",200,0,5,1000,-0.5,0.5);
   hHFTDcaZ = new TH2D("hft_dcaZ",";pt;dcaZ",200,0,5,1000,-0.5,0.5);

   hD0MassLike = new TH2D("d0MassLike","",100,0,10,500,1.6,2.1);
   hD0MassUnlike = new TH2D("d0MassUnlike","",100,0,10,500,1.6,2.1);
   hDecayLength = new TH1D("decayLength","",500,0,0.05);
   hKaonDca = new TH1D("kaonDca","",200,0,0.02);
   hPionDca = new TH1D("pionDca","",200,0,0.02);
   hDcaDaughters = new TH1D("dcaDaughters","",500,0,0.01);
   hPointing = new TH1D("pointing","",500,0,0.01);
   trackTest1 = new TH1D("trackTest1","",10,0,10);
   trackTest2 = new TH1D("trackTest2","",10,0,10);

   tuplePair = new TNtuple("tuplePair","","mass:pt:kaonDca:pionDca:dcaDaughters:pointingAngle:decayLenght:charge:kaonPt:pionPt");
   // -------------- USER VARIABLES -------------------------
   

   hVzVpdVz = new TH2F("hVzVpdVz","hVzVpdVz",200,-100,100,200,-100,100);
   ntpD0BackgroundSameSignPt01 = new TNtuple("ntpD0BackgroundSameSignPt01", "ntpD0BackgroundSameSignPt01", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt01 = new TNtuple("ntpD0BackgroundSideBandPt01", "ntpD0BackgroundSideBandPt01", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt12 = new TNtuple("ntpD0BackgroundSameSignPt12", "ntpD0BackgroundSameSignPt12", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt12 = new TNtuple("ntpD0BackgroundSideBandPt12", "ntpD0BackgroundSideBandPt12", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt23 = new TNtuple("ntpD0BackgroundSameSignPt23", "ntpD0BackgroundSameSignPt23", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt23 = new TNtuple("ntpD0BackgroundSideBandPt23", "ntpD0BackgroundSideBandPt23", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt35 = new TNtuple("ntpD0BackgroundSameSignPt35", "ntpD0BackgroundSameSignPt35", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt35 = new TNtuple("ntpD0BackgroundSideBandPt35", "ntpD0BackgroundSideBandPt35", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSameSignPt510 = new TNtuple("ntpD0BackgroundSameSignPt510", "ntpD0BackgroundSameSignPt510", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");
   ntpD0BackgroundSideBandPt510 = new TNtuple("ntpD0BackgroundSideBandPt510", "ntpD0BackgroundSideBandPt510", "m:pt:decayLength:dca12:dcaV0ToPv:cosTheta:cosThetaStar:ptKaon:dcaKaon:ptPion:dcaPion");


   // for v2 calcualtion
  TString flatten[5];
  flatten[0] = "v2";
  flatten[1] = "cosD";
  flatten[2] = "sinD";
  flatten[3] = "cosHadron";
  flatten[4] = "sinHadron";
  TString sb[8] = {"s1like","s3like","hSBlike","lSBlike","s1unlike","s3unlike","hSBunlike","lSBunlike"};
  // float xbin[7] = {0,1,2,3,4,5,10};
  const int xbinSize=5;
  float xbin[6] = {0,1,2,3,5,10};
  float binMass[2001];
  float binPhi[2001];
  candPt = new TProfile("candPt","",xbinSize,xbin);
  for(int i=0;i<2001;i++)
    binPhi[i] = 0.005*i-5;
  for(int i=0;i<2001;i++)
    binMass[i] = 0.01*i;
  massPt = new TH2D("massPt","",2000,binMass,xbinSize,xbin);
  massPtLike = new TH2D("massPtLike","",2000,binMass,xbinSize,xbin);
  massLike = new TH2D("massLike","",2000,binMass,xbinSize,xbin);
  massLike->Sumw2();
  massUnlike = new TH2D("massUnlike","",2000,binMass,xbinSize,xbin);
  massUnlike->Sumw2();
  float xWeight[6] = {0,7,12,16,22,100};
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        TString name = sb[i]+flatten[j]+Form("_%i",k);
        profV2[i][j][k] = new TProfile(name.Data(),"",xbinSize,xbin);
        profV2[i][j][k]->Sumw2();
      }
      TString weightName = sb[i]+Form("_%i_weigth",k);
      v2Weight[i][k] = new TH2D(weightName.Data(),"",5,xWeight,xbinSize,xbin);
      v2Weight[i][k]->Sumw2();

      TString namehPhi = "hadronPhi_"+sb[i]+Form("_%i",k);
      TString nameDPhi = "DPhi_"+sb[i]+Form("_%i",k);
      // hPhiHadron[i][k] = new TH2F(namehPhi.Data(),"",2000,binPhi,xbinSize,xbin);
      // hPhiD[i][k]= new TH2F(nameDPhi.Data(),"",2000,binPhi,xbinSize,xbin);
      // hPhiD[i][k]->Sumw2();
      // hPhiHadron[i][k]->Sumw2();
    }
  }
  // TH2D *testV2 = new TH2D("testV2","",10,0,10,200,-1,1);

  float ptbin1[12] = {0.225,0.375,0.525,0.675,0.825,0.975,1.12,1.27,1.42,1.58,1.73,1.88};
  float ptbin2[11];
  for(int i=0;i<11;i++)
    ptbin2[i] = 0.5*(ptbin1[i]+ptbin1[i+1]);
  for(int k=0;k<3;k++)
  {
    for(int i=0;i<5;i++)
    {
      hadronV2[i][k] = new TH1D(Form("hadron_%s_%i",flatten[i].Data(),k),"",5,xWeight);
      hadronV2[i][k]->Sumw2();
      hadronV2_sum[i][k] = new TH1D(Form("hadronsum_%s_%i",flatten[i].Data(),k),"",5,xWeight);
      hadronV2_sum[i][k]->Sumw2();
      for(int j=0;j<9;j++)
      {
        hadronV2_excl[i][j][k] = new TH1D(Form("hadron_%s_cent%i_%i",flatten[i].Data(),j,k),"",10,ptbin2);
        hadronV2_excl[i][j][k]->Sumw2();
        hadronV2_excl_sum[i][j][k] = new TH1D(Form("hadronsum_%s_cent%i_%i",flatten[i].Data(),j,k),"",10,ptbin2);
        hadronV2_excl_sum[i][j][k]->Sumw2();
      }
    }
  }

  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  // ifstream ifs("efficiency.txt");
  // for(int i=0; i<6; i++)
  //   for(int j=0; j<4; j++)
  //     ifs>>efficiency[j][i];
  for(int i=0;i<6;i++)//pt bin
  {
    for(int j=0;j<5;j++)//flatten
    {
      TString massName[2];
      massName[0] = Form("likeMass%i",i)+flatten[j];
      massName[1] = Form("unlikeMass%i",i)+flatten[j];
      for(int k=0;k<2;k++)
      {
        V2Mass[k][i][j] = new TProfile(massName[k].Data(),"",18,fitmean[i]-9*fitsigma[i],fitmean[i]+9*fitsigma[i]);
        V2Mass[k][i][j]->Sumw2();
      }
    }
  }

}

//-----------------------------------------------------------------------------
void StMyAnalysisMaker::WriteHistograms() {
  mD0Tuple->Write();
   eventCheck->Write();
  hMult->Write();
   hD0MassLike->Write();
   hD0MassUnlike->Write();
   hDecayLength->Write();
   hKaonDca->Write();
   hPionDca->Write();
   hDcaDaughters->Write();
   hPointing->Write();

   hPVZ->Write();
   deltaVZ->Write();
   hEta->Write();
   hTPCEta->Write();
   hHFTEta->Write();
   hPhi->Write();
   hTPCPhi->Write();
   hHFTPhi->Write();
   hPt->Write();
   hTPCPt->Write();
   hHFTPt->Write();
   hDcaZ->Write();
   hTPCDcaZ->Write();
   hHFTDcaZ->Write();
   hDcaXY->Write();
   hTPCDcaXY->Write();
   hHFTDcaXY->Write();

   tuplePair->Write();
   trackTest1->Write();
   trackTest2->Write();
   hVzVpdVz->Write();
   // testV2->Write();
   // ntpD0BackgroundSameSignPt01->Write();
   // ntpD0BackgroundSideBandPt01->Write();
   // ntpD0BackgroundSameSignPt12->Write();
   // ntpD0BackgroundSideBandPt12->Write();
   // ntpD0BackgroundSameSignPt23->Write();
   // ntpD0BackgroundSideBandPt23->Write();
   // ntpD0BackgroundSameSignPt35->Write();
   // ntpD0BackgroundSideBandPt35->Write();
   // ntpD0BackgroundSameSignPt510->Write();
   // ntpD0BackgroundSideBandPt510->Write();
   //
  //Saving for v2 calculation
  for(int i=0;i!=8;i++)
  {
    for(int k=0;k!=3;k++)
    {
      for(int j=0;j!=5;j++)
      {
        profV2[i][j][k]->Write();
      }
      v2Weight[i][k]->Write();
    }
  }
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<5;j++)
    {
      for(int k=0;k<2;k++)
        V2Mass[k][i][j]->Write();
    }
  }
  massLike->Write();
  candPt->Write();
  massUnlike->Write();
  for(int k=0;k<3;k++)
  {
    for(int i=0;i<5;i++)
    {
      hadronV2[i][k]->Write();
      hadronV2_sum[i][k]->Write();
    }
  }

}

//----------------------------------------------------------------------------- 
void StMyAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StMyAnalysisMaker::Make() {
  if(!mPicoDstMaker) {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }
//////////////////////////////////
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
   StThreeVectorF vtx = event->primaryVertex();
   float b = event->bField();
   double vpdvz = mPicoDst->event()->vzVpd();
   double vz = vtx.z();
   deltaVZ->Fill(event->grefMult(),vz-vpdvz);
  // cout<<"ref mult = "<<event->grefMult()<<endl;
  hMult->Fill(event->grefMult());
  eventCheck->Fill(0);
  // if(!(isGoodEvent()))
  // { 
  //   LOG_WARN << " Not Min Bias! Skip! " << endm;
  //   return kStWarn;
  // }
  // StThreeVectorF pVtx(-999.,-999.,-999.);
  // StPicoKFVertexFitter kfVertexFitter;
  // pVtx = kfVertexFitter.primaryVertexRefit(mPicoDst);

  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  // mGRefMultCorrUtil->init(mPicoDst->event()->runId());
  // mGRefMultCorrUtil->initEvent(mPicoDst->event()->grefMult(),pVtx.z(),mPicoDst->event()->ZDCx()) ;

  // float centrality = (float)mGRefMultCorrUtil->getCentralityBin9();
  bool isGoodTrigger = false;
   std::array<unsigned int, 5> const triggers = { 530003,    // vpdmb-5-p-nobsmd-hlt 
                                                  430903,    // vpdmb-5-p-nobsmd-hlt 
                                                  430904,    // vpdmb-5-p-nobsmd-hlt 
                                                  430906,    // vpdmb-5-p-nobsmd-hlt 
                                                  430907    // vpdmb-5-p-nobsmd-hlt 
                                                };
  for(auto trg: triggers)
  {
    if(event->isTrigger(trg)) isGoodTrigger = true;
  }
  // if(vz==0)
  //    exit(1);

   // if(!isGoodTrigger)
   // {
   //   // cout<<"not good trigger"<<endl;
   //   // exit(1);
   //   return kStWarn;
   // }
   //
   hPVZ->Fill(vz,vpdvz);
   hVzVpdVz->Fill(vz,vpdvz);

   // if( fabs(vz)>10 && fabs(vz-vpdvz)>3 ) 
   //   return kStWarn;

   // getHadronCorV2(0);
   
  eventCheck->Fill(1);
   // cout<<"number of tracks = "<<picoDst->numberOfTracks()<<endl;
   for(int i=0;i<mPicoDst->numberOfTracks();i++)
   {
     StPicoTrack const* trk = mPicoDst->track(i);
     double eta = trk->gMom(vtx, b).pseudoRapidity();
     double phi = trk->gMom(vtx, b).phi() ;
     double pt = trk->gPt();
     StPhysicalHelixD helix = trk->helix(b);
     double dcaZ = (trk->dcaPoint().z() - vtx.z());
     double dcaXY = helix.geometricSignedDistance(vtx.x(),vtx.y());
     double dca = (vtx-trk->dcaPoint()).mag();
     // TOF 
     int tofIndex = trk->bTofPidTraitsIndex();
     bool TofMatch = kFALSE;
     StPicoBTofPidTraits* tofPidTraits;
     if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex); //GNX
     if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
     if(!(trk->isHFTTrack())) continue;
     trackCheckTof->Fill(0);
     if(!TofMatch) continue; // Tof match for pile-up
     trackCheckTof->Fill(1);
     if(trk->nHitsFit() <= 20) continue;
     trackCheckTof->Fill(2);
     if(pt<0.6)  continue; // global pT cut
     trackCheckTof->Fill(3);
     if(dca>1.5)  continue; // dca
     trackCheckTof->Fill(4);
     if (isTpcPion(trk)) 
      trackCheckTof->Fill(5);

     bool tpcKaon = isTpcKaon(trk,&vtx);
     float kBeta = getTofBeta(trk,&vtx);
     bool tofAvailable = kBeta>0;
     bool tofKaon = tofAvailable && isTofKaon(trk,kBeta);
     bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
     if(goodKaon)
       trackCheckTof->Fill(6);
   }
   for(int i=0;i<mPicoDst->numberOfTracks();i++)
   {
     StPicoTrack const* trk = mPicoDst->track(i);
     double eta = trk->gMom(vtx, b).pseudoRapidity();
     double phi = trk->gMom(vtx, b).phi() ;
     double pt = trk->gPt();
     StPhysicalHelixD helix = trk->helix(b);
     double dcaZ = (trk->dcaPoint().z() - vtx.z());
     double dcaXY = helix.geometricSignedDistance(vtx.x(),vtx.y());
     double dca = (vtx-trk->dcaPoint()).mag();
     // TOF 
     int tofIndex = trk->bTofPidTraitsIndex();
     bool TofMatch = kFALSE;
     StPicoBTofPidTraits* tofPidTraits;
     if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex); //GNX
     if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
     if(!(trk->isHFTTrack())) continue;
     trackCheckNoTof->Fill(0);
     // if(!TofMatch) continue; // Tof match for pile-up
     // trackCheckNoTof->Fill(1);
     if(trk->nHitsFit() <= 20) continue;
     trackCheckNoTof->Fill(2);
     if(pt<0.6)  continue; // global pT cut
     trackCheckNoTof->Fill(3);
     if(dca>1.5)  continue; // dca
     trackCheckNoTof->Fill(4);
     if (isTpcPion(trk)) 
      trackCheckNoTof->Fill(5);

     bool tpcKaon = isTpcKaon(trk,&vtx);
     float kBeta = getTofBeta(trk,&vtx);
     bool tofAvailable = kBeta>0;
     bool tofKaon = tofAvailable && isTofKaon(trk,kBeta);
     bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
     if(goodKaon)
       trackCheckNoTof->Fill(6);
   }


   for(int i=0;i<mPicoDst->numberOfTracks();i++)
   {
     StPicoTrack const* itrk = mPicoDst->track(i);
     trackTest1->Fill(1);
     if(!isGoodTrack(itrk))  continue;
     trackTest1->Fill(2);
     if (!isTpcPion(itrk)) continue;
     trackTest1->Fill(3);
     if(!(itrk->isHft())) continue;
     trackTest1->Fill(4);
     for(int j=0;j<mPicoDst->numberOfTracks();j++)
     {
       StPicoTrack const* jtrk = mPicoDst->track(j);
       trackTest2->Fill(1);
       if(!isGoodTrack(jtrk))  continue; // for HFT track 
       trackTest2->Fill(2);
       // if(!(jtrk->isHFTTrack())) continue;
       // trackTest2->Fill(3);
       bool tpcKaon = isTpcKaon(jtrk,&vtx);
       float kBeta = getTofBeta(jtrk,&vtx);
       bool tofAvailable = kBeta>0;
       bool tofKaon = tofAvailable && isTofKaon(jtrk,kBeta);
       bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
       if(!goodKaon) continue;
       trackTest2->Fill(4);
       int charge = itrk->charge() * jtrk->charge();
       // Lukas method
       StThreeVectorF const p1Mom = itrk->gMom();
       StThreeVectorF const p2Mom = jtrk->gMom();

       StPhysicalHelixD const p1StraightLine(p1Mom, itrk->origin(), 0, itrk->charge());
       StPhysicalHelixD const p2StraightLine(p2Mom, jtrk->origin(), 0, jtrk->charge());

       // pair<double, double> const ss = (useStraightLine) ? p1StraightLine.pathLengths(p2StraightLine) : p1Helix.pathLengths(p2Helix);
       pair<double, double> const ss = p1StraightLine.pathLengths(p2StraightLine);
       StThreeVectorF const p1AtDcaToP2 = p1StraightLine.at(ss.first);
       StThreeVectorF const p2AtDcaToP1 = p2StraightLine.at(ss.second);
       double mDcaDaughters = (p1AtDcaToP2 - p2AtDcaToP1).mag();
       //////
       StKaonPion *kp = new StKaonPion(jtrk,itrk,j,i,vtx,b);
       hDecayLength->Fill(kp->decayLength());
       hKaonDca->Fill(kp->kaonDca());
       hPionDca->Fill(kp->pionDca());
       hDcaDaughters->Fill(kp->dcaDaughters());
       // cout<<"find pair"<<endl;
       hPointing->Fill(kp->pointingAngle() * kp->decayLength());

       StPicoTrack const* kaon = mPicoDst->track(kp->kaonIdx());
       StPicoTrack const* pion = mPicoDst->track(kp->pionIdx());
       if(kp->pt()>1 && kp->pt()<2)
         tuplePair->Fill(kp->m(),kp->pt(),kp->kaonDca(),kp->pionDca(),kp->dcaDaughters(),kp->pointingAngle(), kp->decayLength(),charge,kaon->gPt(),pion->gPt());
       if(kaon->charge()*pion->charge()>0 && kp->m()>1.813 && kp->m()<1.919)
       {
         // if(randomSample(kp->pt()))
         {
           if(kp->pt()>0 && kp->pt()<1)
             ntpD0BackgroundSameSignPt01->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>1 && kp->pt()<2)
             ntpD0BackgroundSameSignPt12->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>2 && kp->pt()<3)
             ntpD0BackgroundSameSignPt23->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>3 && kp->pt()<5)
             ntpD0BackgroundSameSignPt35->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>5 && kp->pt()<10)
             ntpD0BackgroundSameSignPt510->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
         }
       }
       if(kaon->charge()*pion->charge()<0 && ((kp->m()>1.761 && kp->m()<1.813) || (kp->m()>1.919 && kp->m()<1.971)))
       {
         // if(randomSample(kp->pt()))
         {
           if(kp->pt()>0 && kp->pt()<1)
             ntpD0BackgroundSideBandPt01->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>1 && kp->pt()<2)
             ntpD0BackgroundSideBandPt12->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>2 && kp->pt()<3)
             ntpD0BackgroundSideBandPt23->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>3 && kp->pt()<5)
             ntpD0BackgroundSideBandPt35->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
           if(kp->pt()>5 && kp->pt()<10)
             ntpD0BackgroundSideBandPt510->Fill(kp->m(), kp->pt(), kp->decayLength(), kp->dcaDaughters(), kp->decayLength()*sin(kp->pointingAngle()), cos(kp->pointingAngle()), kp->cosThetaStar(), kaon->gMom(vtx, b).perp(), kp->kaonDca(), pion->gMom(vtx, b).perp(), kp->pionDca());
         }
       }
 
       if((charge=isD0Pair(kp))!=0 )
       {
         cout<<"find good pair"<<endl;
         if(charge>0)
           hD0MassLike->Fill(kp->pt(),kp->m());
         if(charge==-1)
         {
           hD0MassUnlike->Fill(kp->pt(),kp->m());
           if(kp->m()>1.81&&kp->m()<1.91)
             candPt->Fill(kp->pt(),kp->pt(),1);
         }
         getCorV2(kp,1);
       }
       delete kp;
     }
   }


  return kStOK;
}




// bool StMyAnalysisMaker::isGoodEvent()
// {
//   StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
//      return (event->isHT()&&
//              fabs(event->primaryVertex().z()) < mycuts::vz &&
//            fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz);
//   //return event->triggerWord() & 0x1F;
// }

bool StMyAnalysisMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    // double ptot = trk->dcaGeometry().momentum().mag();
    double ptot = trk->gPtot();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }

  return tofKaon;
}
float StMyAnalysisMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{

  int index2tof = trk->bTofPidTraitsIndex();

  float beta = std::numeric_limits<float>::quiet_NaN();

  if(index2tof >= 0)
  {
    StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

    if(tofPid)
    {
      beta = tofPid->btofBeta();

      if (beta < 1e-4)
      {
        StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        // StPhysicalHelixD helix = trk->helix();
        StPhysicalHelixD helix = trk->helix(mPicoDstMaker->picoDst()->event()->bField());

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  return beta;
}
//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------

bool StMyAnalysisMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  // Tof match for pile-up tracks
  int tofIndex = trk->bTofPidTraitsIndex();
  bool TofMatch = kFALSE;
  StPicoBTofPidTraits* tofPidTraits;
  if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex); //GNX
  if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
  // Dca cut
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  StThreeVectorF vtx = event->primaryVertex();
  double dca = (vtx-trk->dcaPoint()).mag();
  // return dca<1.5 && TofMatch && trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  return dca<1.5 && TofMatch && trk->gPt() > 0.2 && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}

int StMyAnalysisMaker::isD0Pair(StKaonPion const* const kp) const
{

  StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  TLorentzVector d0Lorentz;
  d0Lorentz.SetPtEtaPhiM(kp->pt(),kp->eta(),kp->phi(),kp->m());
  if(fabs(d0Lorentz.Rapidity())>1.) return 0;
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0088 &&
      kp->pionDca() > 0.0112 && kp->kaonDca() > 0.0096 &&
      kp->dcaDaughters() < 0.0093 && kp->decayLength()>0.0106;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0075 &&
      kp->pionDca() > 0.0099 && kp->kaonDca() > 0.0087 &&
      kp->dcaDaughters() < 0.0134 && kp->decayLength()>0.0198;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0102 &&
      kp->pionDca() > 0.0069 && kp->kaonDca() > 0.0074 &&
      kp->dcaDaughters() < 0.0093 && kp->decayLength()>0.0232;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0066 &&
      kp->pionDca() > 0.0074 && kp->kaonDca() > 0.0070 &&
      kp->dcaDaughters() < 0.0093 && kp->decayLength()>0.0187;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0119 &&
      kp->pionDca() > 0.0057 && kp->kaonDca() > 0.0059 &&
      kp->dcaDaughters() < 0.0058 && kp->decayLength()>0.0046;  
  }
   // pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.005 &&
   //    kp->pionDca() > 0.01 && kp->kaonDca() > 0.01 &&
   //    kp->dcaDaughters() < 0.01 && kp->decayLength()>0.015;  
  
  int charge = kaon->charge() * pion->charge();
  if(charge>0)
    charge = kaon->charge()>0 ? 1:2;


  if(pairCuts)
    return charge;
  else
    return 0;
}

float StMyAnalysisMaker::sumCos(float phi,vector<float> &hadronsPhi) 
{
  float sumOfCos = 0;
  for(unsigned int i=0;i<hadronsPhi.size();++i)
  {
    sumOfCos += cos(2*(phi-hadronsPhi[i]));
  }
  return sumOfCos;
}

bool StMyAnalysisMaker::fixPhi(vector<float> &phi) 
{
  if(phi.size() == 0) return false;
  float sumPhi = 0;
  for(unsigned int i=0;i<phi.size();i++)
    sumPhi+=phi[i];
  float meanPhi = sumPhi/phi.size();
  for(unsigned int i=0;i<phi.size();i++)
    phi[i] = phi[i]-meanPhi;  
  return true;
}


bool StMyAnalysisMaker::getHadronCorV2(int idxGap)
{
  double etaGap[3] = {0,0.15,0.05};
  double mEtaGap = etaGap[idxGap];
  float hadronFill[7] = {0};
  const double reweight = 1;//mGRefMultCorrUtil->getWeight();
  // int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  int mult = event->grefMult();
  for(unsigned int i=0;i<mPicoDst->numberOfTracks();++i)
  {
    StPicoTrack const* hadron = mPicoDst->track(i);
    if(hadron->pMom().perp()<0.2) continue;
    if(!isGoodHadron(hadron)) continue;
    float etaHadron = hadron->pMom().pseudoRapidity();
    float phiHadron = hadron->pMom().phi();
    if(etaHadron<-0.5*mEtaGap)//backward sample 
    {
      hadronFill[0]++;
      hadronFill[1] += sin(2 * phiHadron);
      hadronFill[2] += cos(2 * phiHadron);
    }			
    if(etaHadron>0.5*mEtaGap)//forward sample
    {
      hadronFill[3]++;
      hadronFill[4] += sin(2 * phiHadron);
      hadronFill[5] += cos(2 * phiHadron);
    }			
  }
  hadronFill[6] = mult;
  hadronFill[7] = reweight;
  //mHadronTuple->Fill(hadronFill);
  if(hadronFill[0]==0 || hadronFill[3]==0)
    return false; 
  double temp = (hadronFill[1]*hadronFill[4]+hadronFill[2]*hadronFill[5]);
  hadronV2[0][idxGap]->Fill(mult,temp*reweight);
  hadronV2[1][idxGap]->Fill(mult,hadronFill[2]*reweight);
  hadronV2[2][idxGap]->Fill(mult,hadronFill[1]*reweight);
  hadronV2[3][idxGap]->Fill(mult,hadronFill[5]*reweight);
  hadronV2[4][idxGap]->Fill(mult,hadronFill[4]*reweight);
  hadronV2_sum[0][idxGap]->Fill(mult,hadronFill[0]*hadronFill[3]*reweight);
  hadronV2_sum[1][idxGap]->Fill(mult,hadronFill[0]*reweight);
  hadronV2_sum[2][idxGap]->Fill(mult,hadronFill[0]*reweight);
  hadronV2_sum[3][idxGap]->Fill(mult,hadronFill[3]*reweight);
  hadronV2_sum[4][idxGap]->Fill(mult,hadronFill[3]*reweight);
  //    StPicoTrack const* hadron = picoDst->track(i);
  //  hadronV2_excl[0][centrality]->Fill(hadron->pMom().perp(),temp*reweight);
  //  hadronV2_excl[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[2]*reweight);
  //  hadronV2_excl[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[1]*reweight);
  //  hadronV2_excl[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[5]*reweight);
  //  hadronV2_excl[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[4]*reweight);
  //  hadronV2_excl_sum[0][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*hadronFill[3]*reweight);
  //  hadronV2_excl_sum[1][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
  //  hadronV2_excl_sum[2][centrality]->Fill(hadron->pMom().perp(),hadronFill[0]*reweight);
  //  hadronV2_excl_sum[3][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
  //  hadronV2_excl_sum[4][centrality]->Fill(hadron->pMom().perp(),hadronFill[3]*reweight);
  return true;
}

//bool StMyAnalysisMaker::getCorV2(int idxCand,double etaGap, int &flagCor, double weight)
bool StMyAnalysisMaker::getCorV2(StKaonPion *kp,double weight)
{
  // int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  StPicoEvent *event = (StPicoEvent *)mPicoDst->event();
  int mult = event->grefMult();
  // TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();
  // StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idxCand);
  StPicoTrack const* kaon = mPicoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDst->track(kp->pionIdx());
  int charge = kaon->charge() * pion->charge();
  double dMass = kp->m();


  if(kp->pt()>10) return false;
  int ptIdx = 5;
  if(kp->pt()<5)
    ptIdx = static_cast<int>(kp->pt());
  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  // double fitmean[6] = {1.8620,1.8647,1.8617,1.86475,1.8608,1.8626};
  // double fitsigma[6] = {0.0190,0.0143,0.0143,0.0169282,0.0199567,0.0189131};
  double mean = fitmean[ptIdx];
  double sigma = fitsigma[ptIdx];
  bool fillSB[8];
  fillSB[0] =  (charge>0)&& (dMass>(mean-1*sigma)) &&  (dMass<(mean+1*sigma));
  fillSB[1] =  (charge>0)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
  fillSB[2] =  (charge>0) && (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
  fillSB[4] = (charge==-1)&& (((dMass>(mean+5.5*sigma)) &&  (dMass<(mean+7.5*sigma))) ||((dMass>(mean-7.5*sigma)) &&  (dMass<(mean-5.5*sigma))));
  fillSB[5] = (charge==-1)&& (dMass>(mean-3*sigma)) &&  (dMass<(mean+3*sigma));
  fillSB[6] = (charge==-1)&& (((dMass>(mean+4*sigma)) &&  (dMass<(mean+9*sigma))) ||((dMass>(mean-9*sigma)) &&  (dMass<(mean-4*sigma))));
  fillSB[3] = fillSB[1] || fillSB[2];
  fillSB[7] = fillSB[1] || fillSB[2] || fillSB[6];
  double etaGap[3] = {0,0.15,0.05};

  for(int k=0;k<3;k++)
  {
    double corFill[7] = {0};
    corFill[0] = 1 ;
    // corFill[1] = sin(2* kp->phi())/sqrt(hadronv2[centBin]);
    // corFill[2] = cos(2* kp->phi())/sqrt(hadronv2[centBin]);
    corFill[1] = sin(2* kp->phi());
    corFill[2] = cos(2* kp->phi());
    int chargeIdx = charge>0 ? 0:1;
    for(int j=0;j<8;j++)
    {
      if(fillSB[j])
      {
        profV2[j][1][k]->Fill(kp->pt(),corFill[1],weight);
        profV2[j][2][k]->Fill(kp->pt(),corFill[2],weight);
      }
    }
    for(unsigned int i=0; i<mPicoDst->numberOfTracks();i++)
    {
      StPicoTrack const* hadron = mPicoDst->track(i);
      if(hadron->pMom().perp()<0.2) continue;
      if(!isGoodHadron(hadron)) continue;
      if(i==kp->kaonIdx() || i==kp->pionIdx()) continue;
      float etaHadron = hadron->pMom().pseudoRapidity();
      float phiHadron = hadron->pMom().phi();
      if(!isEtaGap(kp->eta(),etaGap[k],etaHadron))  continue;
      corFill[3]++;
      corFill[4] += sin(2 * phiHadron);
      corFill[5] += cos(2 * phiHadron);
      for(int j=0;j<8;j++)
      {
        if(fillSB[j])
        {
          profV2[j][3][k]->Fill(kp->pt(),sin(2*phiHadron),weight);
          profV2[j][4][k]->Fill(kp->pt(),cos(2*phiHadron),weight);
          profV2[j][0][k]->Fill(kp->pt(),cos(2*(phiHadron-kp->phi())),weight);
          v2Weight[j][k]->Fill(mult,kp->pt(),weight);
        }
      }
    }// Loop over charged tracks
    // if(corFill[3]<=0) return false;
    // double cumulant = (corFill[1]*corFill[4]+corFill[2]*corFill[5])/(corFill[3]); 
    // cout<<"kp phi = "<<kp->phi()<<"\tcumulant = "<<cumulant<<endl;
    // cout<<"number of correlated hadrons = "<<corFill[3]<<endl;
    // for(int j=0;j<8;j++)
    // {
    //   if(fillSB[j])
    //   {
    //     // testV2->Fill(kp->pt(),cumulant);
    //   }
    // }
  }//Loop over different eta gap (k)
  return true;
}

bool StMyAnalysisMaker::isEtaGap(double dEta,double mGap,double hEta)
{
  if(mGap == 0) return true;
  double range =  2. - mGap*2;
  // if(dEta> (1.-2*mGap))
  //   return hEta<(dEta-mGap) && hEta>(dEta-mGap-range);
  // else if(dEta<(-1.+2*mGap))
  //   return hEta>(dEta+mGap) && hEta<(dEta+mGap+range);
  // else 
  //   return (hEta>(dEta+mGap) || hEta<(dEta-mGap));
  if(dEta>0)
    return hEta<-mGap;
  else
    return hEta>mGap;
}


//-----------------------------------------------------------------------------
bool StMyAnalysisMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  //return trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1.&&fabs(trk->nSigmaElectron())>3 && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
  int tofIndex = trk->bTofPidTraitsIndex();
  bool TofMatch = kFALSE;
  StPicoBTofPidTraits* tofPidTraits;
  if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex); //GNX
  if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
  return TofMatch && trk->pMom().perp() > mycuts::hadronPtMin &&trk->pMom().perp() < mycuts::hadronPtMax && trk->nHitsFit() >= 15 &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}


// int StMyAnalysisMaker::isD0Pair(StKaonPion const* const kp) const
// {
//
//   StPicoTrack const* kaon = mPicoDst->track(kp->kaonIdx());
//   StPicoTrack const* pion = mPicoDst->track(kp->pionIdx());
//   bool pairCuts = false;
//   if(kp->pt()<1)
//   {
//     pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
//       kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
//       kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;
//   }
//   else if(kp->pt()<2)
//   {
//     pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
//       kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
//       kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;
//   }
//   else if(kp->pt()<3)
//   {
//     pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
//       kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
//       kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;
//   }
//   else if(kp->pt()<5)
//   {
//     pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
//       kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
//       kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;
//   }
//   else
//   {
//     pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
//       kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
//       kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;
//   }
//
//   int charge = kaon->charge() * pion->charge();
//   if(charge>0)
//     charge = kaon->charge()>0 ? 1:2;
//
//
//   if(pairCuts)
//     return charge;
//   else
//     return 0;
// }
