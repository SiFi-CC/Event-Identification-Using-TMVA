#ifndef __EIComp_H_
#define __EIComp_H_ 1
//#include "ComptonCone.hh"
#include "InputReader.hh"
#include "InputReaderGeant.hh"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TMath.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <TLorentzVector.h>

class EIComp /*: public TObject */{

public:
  EIComp();
  EIComp(TString path);
  ~EIComp();

  Bool_t Identify(void);
  Bool_t SetInputReader(void);
  
  Bool_t ReadConfig(TString path);
  Bool_t SaveToFile(TObject* ob);
  
  void Print(void);
  void Clear(void);

private:
  int fEventNumber;
  
  Int_t fStart;     ///< first event number
  Int_t fStop;      ///< last event number
   
  TString fInputName;
  Bool_t fFreshOutput;
  
  //Bool_t fVerbose;
  
  vector<PhysicVar>* fEnergyCluster3;
  vector<PhysicVar>* fEnergyCluster1;
  //Double_t fTotalEnergy0;
  //Double_t fEnergyLoss1;
  vector<PhysicVar>* fEnergyCluster2;
  //Double_t fTotalEnergy1;
  //Double_t fEnergyLoss2;
  
  Int_t fSize2;
  Int_t fSizeScat1;
  Int_t fSizeAbs1;
  Int_t fSizeScat2;
  Int_t fSizeAbs2;
  Int_t fSizeScat3;
  Int_t fSizeAbs3;
  Int_t fSizeScat4;
  Int_t fSizeAbs4;
  //vector<PhysicVar*>* fEnergyCluster3;
  //Double_t fTotalEnergy2;
  //TVector3* fPoint1;
  //PhysicVec* fPointReco1;
  PhysicVec* fPointRecoScat21;
  PhysicVec* fPointRecoAbs22;
  Double_t fTotalEnergy2;
  Double_t fEnergyRecoScat21;
  Double_t fEnergyRecoAbs22;
  
  Float_t fPosX_21;
  Float_t fPosX_22;
  
  Float_t fPosX_1;
  Float_t fPosY_1;
  Float_t fPosZ_1;
  Float_t fPosX_2;
  Float_t fPosY_2;
  Float_t fPosZ_2;
  Float_t fPosX_3;
  Float_t fPosY_3;
  Float_t fPosZ_3;
  Float_t fPosX_4;
  Float_t fPosY_4;
  Float_t fPosZ_4;
  Float_t fPosX_5;
  Float_t fPosY_5;
  Float_t fPosZ_5;
  
  Float_t fPosX_1unc;
  Float_t fPosY_1unc;
  Float_t fPosZ_1unc;
  Float_t fPosX_2unc;
  Float_t fPosY_2unc;
  Float_t fPosZ_2unc;
  Float_t fPosX_3unc;
  Float_t fPosY_3unc;
  Float_t fPosZ_3unc;
  Float_t fPosX_4unc;
  Float_t fPosY_4unc;
  Float_t fPosZ_4unc;
  Float_t fPosX_5unc;
  Float_t fPosY_5unc;
  Float_t fPosZ_5unc;
  

  TVector3* fPosScat;
  TVector3* fPosAbs;

  Float_t fEnergy_1;
  Float_t fEnergy_1unc;
  Float_t fEnergy_2;
  Float_t fEnergy_2unc;
  
  Float_t fEnDiff;
  
  
  Float_t fDiffPosX;
  Float_t fDiffEnergy;
  Float_t fRatioEnergy;
  
  
  Float_t fEnergy_3;
  Float_t fEnergy_3unc;
  Float_t fEnergyC;
  Float_t fEnergyCunc;
  Float_t fEnergySe;
  Float_t fEnergySeunc;
  Float_t fEnergyTh;
  Float_t fEnergyThunc;
  Float_t fEnergyFo;
  Float_t fEnergyFounc;
  
  Float_t fPriEnergy;
  
  
  Int_t fS;
  
  Float_t fTheta;
  
  TLorentzVector* fECII;
  
  TVector3* fPosReco_e;
  TVector3* fPosReco_p;
  PhysicVar* fRecoEnergy_e;
  PhysicVar* fRecoEnergy_p;
  //Double_t fRecoEnergy_e;
  //Double_t fRecoEnergy_p;
/////////////old version///////////////////////
  
//  TVector3* fPosScatReal;
//  TVector3* fPosAbsReal;

////////new version of file /////
  
  vector<TVector3>* fPosScatReal;
  vector<TVector3>* fPosAbsReal;
  
/////////////////////////  
  Float_t fPos_eX;
  Float_t fPos_eY;
  Float_t fPos_eZ;
  
  Float_t fPos_pX;
  Float_t fPos_pY;
  Float_t fPos_pZ;
  
  Float_t fRealEnergy_e;
  Float_t fRealEnergy_p;
  
  vector<PhysicVar>* fEnergy;
  vector<PhysicVec>* fPoint;
  
  
  PhysicVar fEnergyReco0;
  PhysicVar fEnergyReco1;
  PhysicVec fPointReco0;
  PhysicVec fPointReco1;
  
  PhysicVec fPosScatClus;
  PhysicVec fPosAbsClus;
  Double_t fEnergyRe0;
  Double_t fEnergyRe1;
  //TVector3* fPoint2;
  //PhysicVec* fPointReco2;
  //TVector3* fPoint3;
  //PhysicVec* fPointReco3;
  //Int_t fStart;
  //Int_t fStop;
  Bool_t fVerbose;
  Int_t fNentries;
  Int_t fNevS;
  TTree* fEvents;
  TTree* fTree_stat;
  TTree* fTreeM;
  TTree* fTreeS2;
  TTree* fTreeS3;
  TTree* fTreeS4;
  TTree* fTreeS5;
  TTree* fTreeB2;
  TTree* fTreeB3;
  TTree* fTreeB4;
  TTree* fTreeB5;
  TTree* fTreeBB2;
  TTree* fTreeBB3;
  TTree* fTreeBB4;
  TTree* fTreeBB5;
  TTree* fTreedB2;
  TTree* fTreedB3;
  TTree* fTreedB4;
  TTree* fTreedB5;
  TTree* fTree2;
  TTree* fTreeRB;
  TTree* fTreeRB2;
  TTree* fTreeRB3;
  TTree* fTreeRB4;
  TTree* fTreeRB5;
  TTree* fTreeBS2;
  TTree* fTreeBS3;
  TTree* fTreeBS4;
  TTree* fTreeBS5;
  TTree* fTreeSB2;
  TTree* fTreeSB3;
  TTree* fTreeSB4;
  TTree* fTreeSB5;
  TTree* fTreeSB;
  
  TFile* fOutputFile;
  TH1F* fHist;
  TH1F* fReco;
  TH1F* fRe0;
  TH1F* fRecoPos;
  TH1F* fRe0Pos;
  InputReaderGeant* fReader;
  

 // ClassDef(EI, 0)
};

#endif
