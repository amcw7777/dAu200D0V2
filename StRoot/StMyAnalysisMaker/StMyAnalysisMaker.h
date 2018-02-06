#ifndef StMyAnalysisMaker_h
#define StMyAnalysisMaker_h

#include "StMaker.h"
#include "TNtuple.h"


//StKFVertexMaker includes
//#include "TObjArray.h"
#include "StEnumerations.h"
#include "StThreeVectorD.hh"
#include "StThreeVectorF.hh"
//
#include "TSpectrum.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer1D.h"
#include "TCanvas.h"
#include "TH1K.h"
#include "StiMaker/StKFVerticesCollection.h"
#include "StPhysicalHelixD.hh"
class StPrimaryVertex; 
class StEvent;
class StDcaGeometry; 
class KFParticle; 
class StKFVerticesCollection; 
class TSpectrum; 
class StAnneling;
//// 

class StPicoDst;
class StPicoDstMaker;
class StPicoTrack;
class TString;
class TH1F;
class TH2F;
class TH2D;
class TProfile;

class TMinuit;
class StKaonPion;
class StPicoPrescales;
class StRefMultCorr;


#include "TChain.h"
#include "StMaker.h"
#include "StAnaCuts.h"
#include "StThreeVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoD0Event;
class StKaonPion;
class StPicoDstMaker;
class StPicoTrack;
class TH1D;
class TH2D;


class StMyAnalysisMaker : public StMaker {
  public:
    StMyAnalysisMaker(const char *name, StPicoDstMaker *picoMaker, const char *outName);
    virtual ~StMyAnalysisMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();

    void    DeclareHistograms();
    void    WriteHistograms();
/*
*/
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    
    TString    mOutName;
    
    TNtuple*   mD0Tuple;
    TH1D*  mMult;
    StDcaGeometry *dca;
    StDcaGeometry *dcaG;
    int primaryVertexRefit(StThreeVectorD *, vector<int>& daughter);
    bool isGoodEvent();
    bool  isGoodTrack(StPicoTrack const*) const;
    bool  isTpcPion(StPicoTrack const*) const;
    bool  isTpcKaon(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    bool isTofKaon(StPicoTrack const* const, float beta) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const * pVtx) const;
    // int isD0Pair(StKaonPion const*) const;
    int isD0Pair(StKaonPion const* const ) const;

    StRefMultCorr* mGRefMultCorrUtil; 

    TH1D *eventCheck;
    TH1D *trackCheckTof;
    TH1D *trackCheckNoTof;
    TH2D *hPVZ;
    TH2D *deltaVZ;

    TH1D *hEta;
    TH1D *hMult;
    TH1D *hTPCEta;
    TH1D *hHFTEta;
    TH1D *hPhi;
    TH1D *hTPCPhi;
    TH1D *hHFTPhi;
    TH1D *hPt;
    TH1D *hTPCPt;
    TH1D *hHFTPt;
    TH2D *hDcaXY;
    TH2D *hTPCDcaXY;
    TH2D *hHFTDcaXY;
    TH2D *hDcaZ;
    TH2D *hTPCDcaZ;
    TH2D *hHFTDcaZ;

    TH2D *hD0MassLike;
    TH2D *hD0MassUnlike;
    TH1D *hDecayLength;
    TH1D *hKaonDca;
    TH1D *hPionDca;
    TH1D *hDcaDaughters;
    TH1D *hPointing;
    TH1D *trackTest1;
    TH1D *trackTest2;
    TNtuple *tuplePair;
  
    //d0 background ntuples
    TH2F*      hVzVpdVz;
    TNtuple* ntpD0BackgroundSameSignPt01;
    TNtuple* ntpD0BackgroundSideBandPt01;
    TNtuple* ntpD0BackgroundSameSignPt12;
    TNtuple* ntpD0BackgroundSideBandPt12;
    TNtuple* ntpD0BackgroundSameSignPt23;
    TNtuple* ntpD0BackgroundSideBandPt23;
    TNtuple* ntpD0BackgroundSameSignPt35;
    TNtuple* ntpD0BackgroundSideBandPt35;
    TNtuple* ntpD0BackgroundSameSignPt510;
    TNtuple* ntpD0BackgroundSideBandPt510; 
    //d0 v2 calculation
    bool  isGoodHadron(StPicoTrack const*) const;
    float sumCos(float phi, vector<float> &hadronsPhi);
    bool fixPhi(vector<float> &phi);
    bool getHadronCorV2(int );
    //bool getCorV2(int , double, int &);
    bool getCorV2(StKaonPion *, double);
    bool isEtaGap(double, double ,double);
    float getD0CorV2(int *sumPair, vector<const StKaonPion *> cand);

    TProfile *profV2[8][5][3];//i.S or B; j.flatten; k. differetn etaGap
    TH1D *hadronV2[5][3];
    TH1D *hadronV2_sum[5][3];
    TH1D *hadronV2_excl[5][9][3];
    TH1D *hadronV2_excl_sum[5][9][3];
    TH2D *testV2;
    TH2D *fitPhi[6];
    TH2D *massPt;
    TH2D *massPtLike;
    TH2D *massLike;
    TH2D *massLike2;
    TH2D *massUnlike;
    TH2D *v2Weight[8][3];
    TH2D *likeV2Mass[6][5];
    TH2D *likeV2Mass2[6][5];
    TH2D *unlikeV2Mass[6][5];
    TProfile *V2Mass[2][6][5];
    TProfile *candPt;


//
//	//StKFVertexMaker private
    ClassDef(StMyAnalysisMaker, 1)
};

#endif
