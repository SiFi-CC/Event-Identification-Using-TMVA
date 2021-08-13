 
#include <iostream> // Stream declarations
#include <vector>
#include <limits>
#include <string>

#include "TChain.h"
#include "TCut.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"

using namespace std;
/// This macro is used because of two below reasons:
/// 1- To have the second half statistics in the analysis phase. 
/// 2- To use Cross-Validation training, user should define a useful variable called "EventID"
/// to be used as a suitable spectator (a custom splitting variable) when dividing data sample into different k-folds.

void SecondHalfStatistics() {
    
    
    TString outfileName( "PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5SecondhalfModified2sigma.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
    
    TTree* fTree_stat = new TTree("TreeStat", "TreeStat");
    
    TTree* fTreeS2 = new TTree("TreeS2", "TreeS2");
    TTree* fTreeB2 = new TTree("TreeB2", "TreeB2");
    TTree* fTreeBB2 = new TTree("TreeBB2", "TreeBB2");
    TTree* fTreeRB2 = new TTree("TreeRB2", "TreeRB2");
    TTree* fTreeBS2 = new TTree("TreeBS2", "TreeBS2");
    TTree* fTreeSB2 = new TTree("TreeSB2", "TreeSB2");
     
    TTree* fTreeS3 = new TTree("TreeS3", "TreeS3");
    TTree* fTreeB3 = new TTree("TreeB3", "TreeB3");
    TTree* fTreeBB3 = new TTree("TreeBB3", "TreeBB3");
    TTree* fTreedB3 = new TTree("TreedB3", "TreedB3");
    TTree* fTreeRB3 = new TTree("TreeRB3", "TreeRB3");
    TTree* fTreeBS3 = new TTree("TreeBS3", "TreeBS3");
    TTree* fTreeSB3 = new TTree("TreeSB3", "TreeSB3");
    
    TTree* fTreeS4 = new TTree("TreeS4", "TreeS4");
    TTree* fTreeB4 = new TTree("TreeB4", "TreeB4");
    TTree* fTreeBB4 = new TTree("TreeBB4", "TreeBB4");
    TTree* fTreedB4 = new TTree("TreedB4", "TreedB4");
    TTree* fTreeRB4 = new TTree("TreeRB4", "TreeRB4");
    TTree* fTreeBS4 = new TTree("TreeBS4", "TreeBS4");
    TTree* fTreeSB4 = new TTree("TreeSB4", "TreeSB4");
    
    TTree* fTreeS5 = new TTree("TreeS5", "TreeS5");
    TTree* fTreeB5 = new TTree("TreeB5", "TreeB5");
    TTree* fTreeBB5 = new TTree("TreeBB5", "TreeBB5");
    TTree* fTreedB5 = new TTree("TreedB5", "TreedB5");
    TTree* fTreeRB5 = new TTree("TreeRB5", "TreeRB5");
    TTree* fTreeBS5 = new TTree("TreeBS5", "TreeBS5");
    TTree* fTreeSB5 = new TTree("TreeSB5", "TreeSB5");
    
    TTree* fTreeSB = new TTree("TreeSB", "TreeSB");
    
    Float_t Pos_eX, Pos_eY, Pos_eZ, Pos_pX, Pos_pY, Pos_pZ, RealEnergy_e, RealEnergy_p; 
   
    Float_t PosX_scat, PosY_scat, PosZ_scat,Energy_scat;
   
    Float_t PosX_abs, PosY_abs, PosZ_abs, Energy_abs;
   
    Float_t PosX_sec, PosY_sec, PosZ_sec;
    Float_t EnergyCluster_abs, AngularDistribution;
   
    Float_t PosX_scatunc, PosY_scatunc, PosZ_scatunc,Energy_scatunc;
    Float_t PosX_absunc, PosY_absunc, PosZ_absunc, Energy_absunc, EnDiff;
    Float_t PosX_secunc, PosY_secunc, PosZ_secunc;
   
    Float_t EnergyCluster_absunc,EnergySe, EnergySeunc;
   
    Float_t DiffEnergy, RatioEnergy, DiffPosition, PriEnergy, EnergySum, EnergySumunc;
    Int_t count = 0;
    Int_t EventNumber;
    Int_t Multiplicity;
    Float_t PosX_1, PosY_1, PosZ_1, Energy_1;
    Float_t PosX_2, PosY_2, PosZ_2, Energy_2;
    Float_t DEnergy, REnergy, DPosition, PEnergy, EnergyS;
    Int_t eventIDS2 = 0,eventIDB2 = 0,eventIDBB2 = 0,eventIDRB2 = 0,eventIDBS2 = 0;
    Int_t eventIDS3 = 0,eventIDB3 = 0,eventIDdB3 = 0,eventIDBB3 = 0,eventIDRB3 = 0,eventIDBS3 = 0;
    Int_t eventIDS4 = 0,eventIDB4 = 0,eventIDdB4 = 0,eventIDBB4 = 0,eventIDRB4 = 0,eventIDBS4 = 0;
    Int_t eventIDS5 = 0,eventIDB5 = 0,eventIDdB5 = 0,eventIDBB5 = 0,eventIDRB5 = 0,eventIDBS5 = 0;
    Int_t eventIDSB2 = 0, eventIDSB3 = 0, eventIDSB4 = 0, eventIDSB5 = 0, eventIDSB = 0 ;
    
    TLorentzVector *fECII = new TLorentzVector();
    
    fTreeS2->Branch("EventNumber", &EventNumber);
    fTreeS2->Branch("eventID", &eventIDS2);
    fTreeS2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS2->Branch("Pos_eX", &Pos_eX);
    fTreeS2->Branch("Pos_eY", &Pos_eY);
    fTreeS2->Branch("Pos_eZ", &Pos_eZ);
    fTreeS2->Branch("Pos_pX", &Pos_pX);
    fTreeS2->Branch("Pos_pY", &Pos_pY);
    fTreeS2->Branch("Pos_pZ", &Pos_pZ);
    fTreeS2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS2->Branch("PosX_Scat", &PosX_scat);
    fTreeS2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS2->Branch("PosY_Scat", &PosY_scat);
    fTreeS2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeS2->Branch("Pos_Scat", &fPosScat);
    fTreeS2->Branch("Energy_Scat", &Energy_scat);
    fTreeS2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS2->Branch("PosX_Abs", &PosX_abs);
    fTreeS2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS2->Branch("PosY_Abs", &PosY_abs);
    fTreeS2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS2->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeS2->Branch("Pos_Abs", &fPosAbs);
    fTreeS2->Branch("PosX_Sec", &PosX_sec);
    fTreeS2->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeS2->Branch("PosY_Sec", &PosY_sec);
    fTreeS2->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeS2->Branch("PosZ_Sec", &PosZ_sec);
    fTreeS2->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeS2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeS2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeS2->Branch("EnergySecond_abs", &EnergySe);
    fTreeS2->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeS2->Branch("Energy_Abs", &Energy_abs);
    fTreeS2->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeS2->Branch("EnDiff", &EnDiff);
    fTreeS2->Branch("EnergySum", &EnergySum);
    fTreeS2->Branch("EnergySumUnc", &EnergySumunc);
    fTreeS2->Branch("Multiplicity", &Multiplicity);
    fTreeS2->Branch("DiffPosition",  &DiffPosition );
    //fTreeS2->Branch("DiffEnergy", &DiffEnergy);
    fTreeS2->Branch("AngularDistribution", &AngularDistribution);
    fTreeS2->Branch("ECII", &fECII);
    fTreeS2->SetCircular(2000000); 
    
    fTreeB2->Branch("EventNumber", &EventNumber);
    fTreeB2->Branch("eventID", &eventIDB2);
    fTreeB2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB2->Branch("Pos_eX", &Pos_eX);
    fTreeB2->Branch("Pos_eY", &Pos_eY);
    fTreeB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeB2->Branch("Pos_pX", &Pos_pX);
    fTreeB2->Branch("Pos_pY", &Pos_pY);
    fTreeB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB2->Branch("PosX_Scat", &PosX_scat);
    fTreeB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB2->Branch("PosY_Scat", &PosY_scat);
    fTreeB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeB2->Branch("Pos_Scat", &fPosScat);
    fTreeB2->Branch("Energy_Scat", &Energy_scat);
    fTreeB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB2->Branch("PosX_Abs", &PosX_abs);
    fTreeB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB2->Branch("PosY_Abs", &PosY_abs);
    fTreeB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeB2->Branch("Pos_Abs", &fPosAbs);
    fTreeB2->Branch("PosX_Sec", &PosX_sec);
    fTreeB2->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeB2->Branch("PosY_Sec", &PosY_sec);
    fTreeB2->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeB2->Branch("PosZ_Sec", &PosZ_sec);
    fTreeB2->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeB2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeB2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeB2->Branch("EnergySecond_abs", &EnergySe);
    fTreeB2->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeB2->Branch("Energy_Abs", &Energy_abs);
    fTreeB2->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeB2->Branch("EnDiff", &EnDiff);
    fTreeB2->Branch("EnergySum", &EnergySum);
    fTreeB2->Branch("EnergySumUnc", &EnergySumunc);
    fTreeB2->Branch("Multiplicity", &Multiplicity);
    fTreeB2->Branch("DiffPosition",  &DiffPosition );
    //fTreeB2->Branch("DiffEnergy", &DiffEnergy);
    fTreeB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeB2->Branch("ECII", &fECII);
    fTreeB2->SetCircular(2000000);
    
    fTreeBB2->Branch("EventNumber", &EventNumber);
    fTreeBB2->Branch("eventID", &eventIDBB2);
    fTreeBB2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB2->Branch("Pos_eX", &Pos_eX);
    fTreeBB2->Branch("Pos_eY", &Pos_eY);
    fTreeBB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB2->Branch("Pos_pX", &Pos_pX);
    fTreeBB2->Branch("Pos_pY", &Pos_pY);
    fTreeBB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB2->Branch("PosX_Scat", &PosX_scat);
    fTreeBB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB2->Branch("PosY_Scat", &PosY_scat);
    fTreeBB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBB2->Branch("Pos_Scat", &fPosScat);
    fTreeBB2->Branch("Energy_Scat", &Energy_scat);
    fTreeBB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB2->Branch("PosX_Abs", &PosX_abs);
    fTreeBB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB2->Branch("PosY_Abs", &PosY_abs);
    fTreeBB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBB2->Branch("Pos_Abs", &fPosAbs);
    fTreeBB2->Branch("PosX_Sec", &PosX_sec);
    fTreeBB2->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBB2->Branch("PosY_Sec", &PosY_sec);
    fTreeBB2->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBB2->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBB2->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBB2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBB2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBB2->Branch("EnergySecond_abs", &EnergySe);
    fTreeBB2->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBB2->Branch("Energy_Abs", &Energy_abs);
    fTreeBB2->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBB2->Branch("EnDiff", &EnDiff);
    fTreeBB2->Branch("EnergySum", &EnergySum);
    fTreeBB2->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBB2->Branch("Multiplicity", &Multiplicity);
    fTreeBB2->Branch("DiffPosition",  &DiffPosition );
    //fTreeBB2->Branch("DiffEnergy", &DiffEnergy);
    fTreeBB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB2->Branch("ECII", &fECII);
    fTreeBB2->SetCircular(2000000); 
    
    fTreeRB2->Branch("EventNumber", &EventNumber);
    fTreeRB2->Branch("eventID", &eventIDRB2);
    fTreeRB2->Branch("PrimaryEnergy", &PriEnergy);
    //fTreeRB2->Branch("PERB", &PERB);
    fTreeRB2->Branch("Pos_eX", &Pos_eX);
    fTreeRB2->Branch("Pos_eY", &Pos_eY);
    fTreeRB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB2->Branch("Pos_pX", &Pos_pX);
    fTreeRB2->Branch("Pos_pY", &Pos_pY);
    fTreeRB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB2->Branch("PosX_Scat", &PosX_scat);
    fTreeRB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB2->Branch("PosY_Scat", &PosY_scat);
    fTreeRB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeRB2->Branch("Pos_Scat", &fPosScat);
    fTreeRB2->Branch("Energy_Scat", &Energy_scat);
    fTreeRB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB2->Branch("PosX_Abs", &PosX_abs);
    fTreeRB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB2->Branch("PosY_Abs", &PosY_abs);
    fTreeRB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeRB2->Branch("Pos_Abs", &fPosAbs);
    fTreeRB2->Branch("PosX_Sec", &PosX_sec);
    fTreeRB2->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeRB2->Branch("PosY_Sec", &PosY_sec);
    fTreeRB2->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeRB2->Branch("PosZ_Sec", &PosZ_sec);
    fTreeRB2->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeRB2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeRB2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeRB2->Branch("EnergySecond_abs", &EnergySe);
    fTreeRB2->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeRB2->Branch("Energy_Abs", &Energy_abs);
    fTreeRB2->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeRB2->Branch("EnDiff", &EnDiff);
    fTreeRB2->Branch("EnergySum", &EnergySum);
    fTreeRB2->Branch("EnergySumUnc", &EnergySumunc);
    fTreeRB2->Branch("Multiplicity", &Multiplicity);
    fTreeRB2->Branch("DiffPosition",  &DiffPosition );
    //fTreeRB2->Branch("DiffEnergy", &DiffEnergy);
    fTreeRB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB2->Branch("ECII", &fECII);
    fTreeRB2->SetCircular(2000000);
    
    fTreeBS2->Branch("EventNumber", &EventNumber);
    fTreeBS2->Branch("eventID", &eventIDBS2);
    fTreeBS2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS2->Branch("Pos_eX", &Pos_eX);
    fTreeBS2->Branch("Pos_eY", &Pos_eY);
    fTreeBS2->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS2->Branch("Pos_pX", &Pos_pX);
    fTreeBS2->Branch("Pos_pY", &Pos_pY);
    fTreeBS2->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS2->Branch("PosX_Scat", &PosX_scat);
    fTreeBS2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS2->Branch("PosY_Scat", &PosY_scat);
    fTreeBS2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBS2->Branch("Pos_Scat", &fPosScat);
    fTreeBS2->Branch("Energy_Scat", &Energy_scat);
    fTreeBS2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS2->Branch("PosX_Abs", &PosX_abs);
    fTreeBS2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS2->Branch("PosY_Abs", &PosY_abs);
    fTreeBS2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS2->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBS2->Branch("Pos_Abs", &fPosAbs);
    fTreeBS2->Branch("PosX_Sec", &PosX_sec);
    fTreeBS2->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBS2->Branch("PosY_Sec", &PosY_sec);
    fTreeBS2->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBS2->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBS2->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBS2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBS2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBS2->Branch("EnergySecond_abs", &EnergySe);
    fTreeBS2->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBS2->Branch("Energy_Abs", &Energy_abs);
    fTreeBS2->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBS2->Branch("EnDiff", &EnDiff);
    fTreeBS2->Branch("EnergySum", &EnergySum);
    fTreeBS2->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBS2->Branch("Multiplicity", &Multiplicity);
    fTreeBS2->Branch("DiffPosition",  &DiffPosition );
    //fTreeBS2->Branch("DiffEnergy", &DiffEnergy);
    fTreeBS2->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS2->Branch("ECII", &fECII);
    fTreeBS2->SetCircular(2000000);
    
    
    fTreeS3->Branch("EventNumber", &EventNumber);
    fTreeS3->Branch("eventID", &eventIDS3);
    fTreeS3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS3->Branch("Pos_eX", &Pos_eX);
    fTreeS3->Branch("Pos_eY", &Pos_eY);
    fTreeS3->Branch("Pos_eZ", &Pos_eZ);
    fTreeS3->Branch("Pos_pX", &Pos_pX);
    fTreeS3->Branch("Pos_pY", &Pos_pY);
    fTreeS3->Branch("Pos_pZ", &Pos_pZ);
    fTreeS3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS3->Branch("PosX_Scat", &PosX_scat);
    fTreeS3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS3->Branch("PosY_Scat", &PosY_scat);
    fTreeS3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeS3->Branch("Pos_Scat", &fPosScat);
    fTreeS3->Branch("Energy_Scat", &Energy_scat);
    fTreeS3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS3->Branch("PosX_Abs", &PosX_abs);
    fTreeS3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS3->Branch("PosY_Abs", &PosY_abs);
    fTreeS3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeS3->Branch("Pos_Abs", &fPosAbs);
    fTreeS3->Branch("PosX_Sec", &PosX_sec);
    fTreeS3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeS3->Branch("PosY_Sec", &PosY_sec);
    fTreeS3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeS3->Branch("PosZ_Sec", &PosZ_sec);
    fTreeS3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeS3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeS3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeS3->Branch("EnergySecond_abs", &EnergySe);
    fTreeS3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeS3->Branch("Energy_Abs", &Energy_abs);
    fTreeS3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeS3->Branch("EnDiff", &EnDiff);
    fTreeS3->Branch("EnergySum", &EnergySum);
    fTreeS3->Branch("EnergySumUnc", &EnergySumunc);
    fTreeS3->Branch("Multiplicity", &Multiplicity);
    fTreeS3->Branch("DiffPosition",  &DiffPosition );
    //fTreeS3->Branch("DiffEnergy", &DiffEnergy);
    fTreeS3->Branch("AngularDistribution", &AngularDistribution);
    fTreeS3->Branch("ECII", &fECII);
    fTreeS3->SetCircular(2000000); 
    
    fTreeB3->Branch("EventNumber", &EventNumber);
    fTreeB3->Branch("eventID", &eventIDB3);
    fTreeB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB3->Branch("Pos_eX", &Pos_eX);
    fTreeB3->Branch("Pos_eY", &Pos_eY);
    fTreeB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeB3->Branch("Pos_pX", &Pos_pX);
    fTreeB3->Branch("Pos_pY", &Pos_pY);
    fTreeB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB3->Branch("PosX_Scat", &PosX_scat);
    fTreeB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB3->Branch("PosY_Scat", &PosY_scat);
    fTreeB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeB3->Branch("Pos_Scat", &fPosScat);
    fTreeB3->Branch("Energy_Scat", &Energy_scat);
    fTreeB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB3->Branch("PosX_Abs", &PosX_abs);
    fTreeB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB3->Branch("PosY_Abs", &PosY_abs);
    fTreeB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeB3->Branch("Pos_Abs", &fPosAbs);
    fTreeB3->Branch("PosX_Sec", &PosX_sec);
    fTreeB3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeB3->Branch("PosY_Sec", &PosY_sec);
    fTreeB3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeB3->Branch("PosZ_Sec", &PosZ_sec);
    fTreeB3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeB3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeB3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeB3->Branch("EnergySecond_abs", &EnergySe);
    fTreeB3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeB3->Branch("Energy_Abs", &Energy_abs);
    fTreeB3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeB3->Branch("EnDiff", &EnDiff);
    fTreeB3->Branch("EnergySum", &EnergySum);
    fTreeB3->Branch("EnergySumUnc", &EnergySumunc);
    fTreeB3->Branch("Multiplicity", &Multiplicity);
    fTreeB3->Branch("DiffPosition",  &DiffPosition );
    //fTreeB3->Branch("DiffEnergy", &DiffEnergy);
    fTreeB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeB3->Branch("ECII", &fECII);
    fTreeB3->SetCircular(2000000);
    
    fTreeBB3->Branch("EventNumber", &EventNumber);
    fTreeBB3->Branch("eventID", &eventIDBB3);
    fTreeBB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB3->Branch("Pos_eX", &Pos_eX);
    fTreeBB3->Branch("Pos_eY", &Pos_eY);
    fTreeBB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB3->Branch("Pos_pX", &Pos_pX);
    fTreeBB3->Branch("Pos_pY", &Pos_pY);
    fTreeBB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB3->Branch("PosX_Scat", &PosX_scat);
    fTreeBB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB3->Branch("PosY_Scat", &PosY_scat);
    fTreeBB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBB3->Branch("Pos_Scat", &fPosScat);
    fTreeBB3->Branch("Energy_Scat", &Energy_scat);
    fTreeBB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB3->Branch("PosX_Abs", &PosX_abs);
    fTreeBB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB3->Branch("PosY_Abs", &PosY_abs);
    fTreeBB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBB3->Branch("Pos_Abs", &fPosAbs);
    fTreeBB3->Branch("PosX_Sec", &PosX_sec);
    fTreeBB3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBB3->Branch("PosY_Sec", &PosY_sec);
    fTreeBB3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBB3->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBB3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBB3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBB3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBB3->Branch("EnergySecond_abs", &EnergySe);
    fTreeBB3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBB3->Branch("Energy_Abs", &Energy_abs);
    fTreeBB3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBB3->Branch("EnDiff", &EnDiff);
    fTreeBB3->Branch("EnergySum", &EnergySum);
    fTreeBB3->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBB3->Branch("Multiplicity", &Multiplicity);
    fTreeBB3->Branch("DiffPosition",  &DiffPosition );
    //fTreeBB3->Branch("DiffEnergy", &DiffEnergy);
    fTreeBB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB3->Branch("ECII", &fECII);
    fTreeBB3->SetCircular(2000000);
    
    fTreedB3->Branch("EventNumber", &EventNumber);
    fTreedB3->Branch("eventID", &eventIDdB3);
    fTreedB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreedB3->Branch("Pos_eX", &Pos_eX);
    fTreedB3->Branch("Pos_eY", &Pos_eY);
    fTreedB3->Branch("Pos_eZ", &Pos_eZ);
    fTreedB3->Branch("Pos_pX", &Pos_pX);
    fTreedB3->Branch("Pos_pY", &Pos_pY);
    fTreedB3->Branch("Pos_pZ", &Pos_pZ);
    fTreedB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreedB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreedB3->Branch("PosX_Scat", &PosX_scat);
    fTreedB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreedB3->Branch("PosY_Scat", &PosY_scat);
    fTreedB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreedB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreedB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreedB3->Branch("Pos_Scat", &fPosScat);
    fTreedB3->Branch("Energy_Scat", &Energy_scat);
    fTreedB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreedB3->Branch("PosX_Abs", &PosX_abs);
    fTreedB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreedB3->Branch("PosY_Abs", &PosY_abs);
    fTreedB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreedB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreedB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreedB3->Branch("Pos_Abs", &fPosAbs);
    fTreedB3->Branch("PosX_Sec", &PosX_sec);
    fTreedB3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreedB3->Branch("PosY_Sec", &PosY_sec);
    fTreedB3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreedB3->Branch("PosZ_Sec", &PosZ_sec);
    fTreedB3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreedB3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreedB3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreedB3->Branch("EnergySecond_abs", &EnergySe);
    fTreedB3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreedB3->Branch("Energy_Abs", &Energy_abs);
    fTreedB3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreedB3->Branch("EnDiff", &EnDiff);
    fTreedB3->Branch("EnergySum", &EnergySum);
    fTreedB3->Branch("EnergySumUnc", &EnergySumunc);
    fTreedB3->Branch("Multiplicity", &Multiplicity);
    fTreedB3->Branch("DiffPosition",  &DiffPosition );
    //fTreedB3->Branch("DiffEnergy", &DiffEnergy);
    fTreedB3->Branch("AngularDistribution", &AngularDistribution);
    fTreedB3->Branch("ECII", &fECII);
    fTreedB3->SetCircular(2000000);
    
    fTreeRB3->Branch("EventNumber", &EventNumber);
    fTreeRB3->Branch("eventID", &eventIDRB3);
    fTreeRB3->Branch("PrimaryEnergy", &PriEnergy);
    //fTreeRB3->Branch("PERB", &PERB);
    fTreeRB3->Branch("Pos_eX", &Pos_eX);
    fTreeRB3->Branch("Pos_eY", &Pos_eY);
    fTreeRB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB3->Branch("Pos_pX", &Pos_pX);
    fTreeRB3->Branch("Pos_pY", &Pos_pY);
    fTreeRB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB3->Branch("PosX_Scat", &PosX_scat);
    fTreeRB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB3->Branch("PosY_Scat", &PosY_scat);
    fTreeRB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeRB3->Branch("Pos_Scat", &fPosScat);
    fTreeRB3->Branch("Energy_Scat", &Energy_scat);
    fTreeRB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB3->Branch("PosX_Abs", &PosX_abs);
    fTreeRB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB3->Branch("PosY_Abs", &PosY_abs);
    fTreeRB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeRB3->Branch("Pos_Abs", &fPosAbs);
    fTreeRB3->Branch("PosX_Sec", &PosX_sec);
    fTreeRB3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeRB3->Branch("PosY_Sec", &PosY_sec);
    fTreeRB3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeRB3->Branch("PosZ_Sec", &PosZ_sec);
    fTreeRB3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeRB3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeRB3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeRB3->Branch("EnergySecond_abs", &EnergySe);
    fTreeRB3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeRB3->Branch("Energy_Abs", &Energy_abs);
    fTreeRB3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeRB3->Branch("EnDiff", &EnDiff);
    fTreeRB3->Branch("EnergySum", &EnergySum);
    fTreeRB3->Branch("EnergySumUnc", &EnergySumunc);
    fTreeRB3->Branch("Multiplicity", &Multiplicity);
    fTreeRB3->Branch("DiffPosition",  &DiffPosition );
    //fTreeRB3->Branch("DiffEnergy", &DiffEnergy);
    fTreeRB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB3->Branch("ECII", &fECII);
    fTreeRB3->SetCircular(2000000);
    
    
    fTreeBS3->Branch("EventNumber", &EventNumber);
    fTreeBS3->Branch("eventID", &eventIDBS3);
    fTreeBS3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS3->Branch("Pos_eX", &Pos_eX);
    fTreeBS3->Branch("Pos_eY", &Pos_eY);
    fTreeBS3->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS3->Branch("Pos_pX", &Pos_pX);
    fTreeBS3->Branch("Pos_pY", &Pos_pY);
    fTreeBS3->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS3->Branch("PosX_Scat", &PosX_scat);
    fTreeBS3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS3->Branch("PosY_Scat", &PosY_scat);
    fTreeBS3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBS3->Branch("Pos_Scat", &fPosScat);
    fTreeBS3->Branch("Energy_Scat", &Energy_scat);
    fTreeBS3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS3->Branch("PosX_Abs", &PosX_abs);
    fTreeBS3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS3->Branch("PosY_Abs", &PosY_abs);
    fTreeBS3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBS3->Branch("Pos_Abs", &fPosAbs);
    fTreeBS3->Branch("PosX_Sec", &PosX_sec);
    fTreeBS3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBS3->Branch("PosY_Sec", &PosY_sec);
    fTreeBS3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBS3->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBS3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBS3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBS3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBS3->Branch("EnergySecond_abs", &EnergySe);
    fTreeBS3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBS3->Branch("Energy_Abs", &Energy_abs);
    fTreeBS3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBS3->Branch("EnDiff", &EnDiff);
    fTreeBS3->Branch("EnergySum", &EnergySum);
    fTreeBS3->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBS3->Branch("Multiplicity", &Multiplicity);
    fTreeBS3->Branch("DiffPosition",  &DiffPosition );
    //fTreeBS3->Branch("DiffEnergy", &DiffEnergy);
    fTreeBS3->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS3->Branch("ECII", &fECII);
    fTreeBS3->SetCircular(2000000);
    
    
    
    fTreeS4->Branch("EventNumber", &EventNumber);
    fTreeS4->Branch("eventID", &eventIDS4);
    fTreeS4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS4->Branch("Pos_eX", &Pos_eX);
    fTreeS4->Branch("Pos_eY", &Pos_eY);
    fTreeS4->Branch("Pos_eZ", &Pos_eZ);
    fTreeS4->Branch("Pos_pX", &Pos_pX);
    fTreeS4->Branch("Pos_pY", &Pos_pY);
    fTreeS4->Branch("Pos_pZ", &Pos_pZ);
    fTreeS4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS4->Branch("PosX_Scat", &PosX_scat);
    fTreeS4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS4->Branch("PosY_Scat", &PosY_scat);
    fTreeS4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeS4->Branch("Pos_Scat", &fPosScat);
    fTreeS4->Branch("Energy_Scat", &Energy_scat);
    fTreeS4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS4->Branch("PosX_Abs", &PosX_abs);
    fTreeS4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS4->Branch("PosY_Abs", &PosY_abs);
    fTreeS4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeS4->Branch("Pos_Abs", &fPosAbs);
    fTreeS4->Branch("PosX_Sec", &PosX_sec);
    fTreeS4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeS4->Branch("PosY_Sec", &PosY_sec);
    fTreeS4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeS4->Branch("PosZ_Sec", &PosZ_sec);
    fTreeS4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeS4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeS4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeS4->Branch("EnergySecond_abs", &EnergySe);
    fTreeS4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeS4->Branch("Energy_Abs", &Energy_abs);
    fTreeS4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeS4->Branch("EnDiff", &EnDiff);
    fTreeS4->Branch("EnergySum", &EnergySum);
    fTreeS4->Branch("EnergySumUnc", &EnergySumunc);
    fTreeS4->Branch("Multiplicity", &Multiplicity);
    fTreeS4->Branch("DiffPosition",  &DiffPosition );
    //fTreeS4->Branch("DiffEnergy", &DiffEnergy);
    fTreeS4->Branch("AngularDistribution", &AngularDistribution);
    fTreeS4->Branch("ECII", &fECII);
    fTreeS4->SetCircular(2000000); 
    
    fTreeB4->Branch("EventNumber", &EventNumber);
    fTreeB4->Branch("eventID", &eventIDB4);
    fTreeB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB4->Branch("Pos_eX", &Pos_eX);
    fTreeB4->Branch("Pos_eY", &Pos_eY);
    fTreeB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeB4->Branch("Pos_pX", &Pos_pX);
    fTreeB4->Branch("Pos_pY", &Pos_pY);
    fTreeB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB4->Branch("PosX_Scat", &PosX_scat);
    fTreeB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB4->Branch("PosY_Scat", &PosY_scat);
    fTreeB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeB4->Branch("Pos_Scat", &fPosScat);
    fTreeB4->Branch("Energy_Scat", &Energy_scat);
    fTreeB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB4->Branch("PosX_Abs", &PosX_abs);
    fTreeB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB4->Branch("PosY_Abs", &PosY_abs);
    fTreeB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeB4->Branch("Pos_Abs", &fPosAbs);
    fTreeB4->Branch("PosX_Sec", &PosX_sec);
    fTreeB4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeB4->Branch("PosY_Sec", &PosY_sec);
    fTreeB4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeB4->Branch("PosZ_Sec", &PosZ_sec);
    fTreeB4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeB4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeB4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeB4->Branch("EnergySecond_abs", &EnergySe);
    fTreeB4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeB4->Branch("Energy_Abs", &Energy_abs);
    fTreeB4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeB4->Branch("EnDiff", &EnDiff);
    fTreeB4->Branch("EnergySum", &EnergySum);
    fTreeB4->Branch("EnergySumUnc", &EnergySumunc);
    fTreeB4->Branch("Multiplicity", &Multiplicity);
    fTreeB4->Branch("DiffPosition",  &DiffPosition );
    //fTreeB4->Branch("DiffEnergy", &DiffEnergy);
    fTreeB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeB4->Branch("ECII", &fECII);
    fTreeB4->SetCircular(2000000);
    
    
    fTreeBB4->Branch("EventNumber", &EventNumber);
    fTreeBB4->Branch("eventID", &eventIDBB4);
    fTreeBB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB4->Branch("Pos_eX", &Pos_eX);
    fTreeBB4->Branch("Pos_eY", &Pos_eY);
    fTreeBB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB4->Branch("Pos_pX", &Pos_pX);
    fTreeBB4->Branch("Pos_pY", &Pos_pY);
    fTreeBB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB4->Branch("PosX_Scat", &PosX_scat);
    fTreeBB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB4->Branch("PosY_Scat", &PosY_scat);
    fTreeBB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBB4->Branch("Pos_Scat", &fPosScat);
    fTreeBB4->Branch("Energy_Scat", &Energy_scat);
    fTreeBB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB4->Branch("PosX_Abs", &PosX_abs);
    fTreeBB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB4->Branch("PosY_Abs", &PosY_abs);
    fTreeBB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBB4->Branch("Pos_Abs", &fPosAbs);
    fTreeBB4->Branch("PosX_Sec", &PosX_sec);
    fTreeBB4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBB4->Branch("PosY_Sec", &PosY_sec);
    fTreeBB4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBB4->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBB4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBB4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBB4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBB4->Branch("EnergySecond_abs", &EnergySe);
    fTreeBB4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBB4->Branch("Energy_Abs", &Energy_abs);
    fTreeBB4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBB4->Branch("EnDiff", &EnDiff);
    fTreeBB4->Branch("EnergySum", &EnergySum);
    fTreeBB4->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBB4->Branch("Multiplicity", &Multiplicity);
    fTreeBB4->Branch("DiffPosition",  &DiffPosition );
    //fTreeBB4->Branch("DiffEnergy", &DiffEnergy);
    fTreeBB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB4->Branch("ECII", &fECII);
    fTreeBB4->SetCircular(2000000);
    
    fTreedB4->Branch("EventNumber", &EventNumber);
    fTreedB4->Branch("eventID", &eventIDdB4);
    fTreedB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreedB4->Branch("Pos_eX", &Pos_eX);
    fTreedB4->Branch("Pos_eY", &Pos_eY);
    fTreedB4->Branch("Pos_eZ", &Pos_eZ);
    fTreedB4->Branch("Pos_pX", &Pos_pX);
    fTreedB4->Branch("Pos_pY", &Pos_pY);
    fTreedB4->Branch("Pos_pZ", &Pos_pZ);
    fTreedB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreedB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreedB4->Branch("PosX_Scat", &PosX_scat);
    fTreedB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreedB4->Branch("PosY_Scat", &PosY_scat);
    fTreedB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreedB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreedB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreedB4->Branch("Pos_Scat", &fPosScat);
    fTreedB4->Branch("Energy_Scat", &Energy_scat);
    fTreedB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreedB4->Branch("PosX_Abs", &PosX_abs);
    fTreedB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreedB4->Branch("PosY_Abs", &PosY_abs);
    fTreedB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreedB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreedB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreedB4->Branch("Pos_Abs", &fPosAbs);
    fTreedB4->Branch("PosX_Sec", &PosX_sec);
    fTreedB4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreedB4->Branch("PosY_Sec", &PosY_sec);
    fTreedB4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreedB4->Branch("PosZ_Sec", &PosZ_sec);
    fTreedB4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreedB4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreedB4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreedB4->Branch("EnergySecond_abs", &EnergySe);
    fTreedB4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreedB4->Branch("Energy_Abs", &Energy_abs);
    fTreedB4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreedB4->Branch("EnDiff", &EnDiff);
    fTreedB4->Branch("EnergySum", &EnergySum);
    fTreedB4->Branch("EnergySumUnc", &EnergySumunc);
    fTreedB4->Branch("Multiplicity", &Multiplicity);
    fTreedB4->Branch("DiffPosition",  &DiffPosition );
    //fTreedB4->Branch("DiffEnergy", &DiffEnergy);
    fTreedB4->Branch("AngularDistribution", &AngularDistribution);
    fTreedB4->Branch("ECII", &fECII);
    fTreedB4->SetCircular(2000000);
    
    fTreeRB4->Branch("EventNumber", &EventNumber);
    fTreeRB4->Branch("eventID", &eventIDRB4);
    fTreeRB4->Branch("PrimaryEnergy", &PriEnergy);
    //fTreeRB4->Branch("PERB", &PERB);
    fTreeRB4->Branch("Pos_eX", &Pos_eX);
    fTreeRB4->Branch("Pos_eY", &Pos_eY);
    fTreeRB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB4->Branch("Pos_pX", &Pos_pX);
    fTreeRB4->Branch("Pos_pY", &Pos_pY);
    fTreeRB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB4->Branch("PosX_Scat", &PosX_scat);
    fTreeRB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB4->Branch("PosY_Scat", &PosY_scat);
    fTreeRB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeRB4->Branch("Pos_Scat", &fPosScat);
    fTreeRB4->Branch("Energy_Scat", &Energy_scat);
    fTreeRB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB4->Branch("PosX_Abs", &PosX_abs);
    fTreeRB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB4->Branch("PosY_Abs", &PosY_abs);
    fTreeRB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeRB4->Branch("Pos_Abs", &fPosAbs);
    fTreeRB4->Branch("PosX_Sec", &PosX_sec);
    fTreeRB4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeRB4->Branch("PosY_Sec", &PosY_sec);
    fTreeRB4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeRB4->Branch("PosZ_Sec", &PosZ_sec);
    fTreeRB4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeRB4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeRB4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeRB4->Branch("EnergySecond_abs", &EnergySe);
    fTreeRB4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeRB4->Branch("Energy_Abs", &Energy_abs);
    fTreeRB4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeRB4->Branch("EnDiff", &EnDiff);
    fTreeRB4->Branch("EnergySum", &EnergySum);
    fTreeRB4->Branch("EnergySumUnc", &EnergySumunc);
    fTreeRB4->Branch("Multiplicity", &Multiplicity);
    fTreeRB4->Branch("DiffPosition",  &DiffPosition );
    //fTreeRB4->Branch("DiffEnergy", &DiffEnergy);
    fTreeRB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB4->Branch("ECII", &fECII);
    fTreeRB4->SetCircular(2000000);
    
    
    fTreeBS4->Branch("EventNumber", &EventNumber);
    fTreeBS4->Branch("eventID", &eventIDBS4);
    fTreeBS4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS4->Branch("Pos_eX", &Pos_eX);
    fTreeBS4->Branch("Pos_eY", &Pos_eY);
    fTreeBS4->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS4->Branch("Pos_pX", &Pos_pX);
    fTreeBS4->Branch("Pos_pY", &Pos_pY);
    fTreeBS4->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS4->Branch("PosX_Scat", &PosX_scat);
    fTreeBS4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS4->Branch("PosY_Scat", &PosY_scat);
    fTreeBS4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBS4->Branch("Pos_Scat", &fPosScat);
    fTreeBS4->Branch("Energy_Scat", &Energy_scat);
    fTreeBS4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS4->Branch("PosX_Abs", &PosX_abs);
    fTreeBS4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS4->Branch("PosY_Abs", &PosY_abs);
    fTreeBS4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBS4->Branch("Pos_Abs", &fPosAbs);
    fTreeBS4->Branch("PosX_Sec", &PosX_sec);
    fTreeBS4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBS4->Branch("PosY_Sec", &PosY_sec);
    fTreeBS4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBS4->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBS4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBS4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBS4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBS4->Branch("EnergySecond_abs", &EnergySe);
    fTreeBS4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBS4->Branch("Energy_Abs", &Energy_abs);
    fTreeBS4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBS4->Branch("EnDiff", &EnDiff);
    fTreeBS4->Branch("EnergySum", &EnergySum);
    fTreeBS4->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBS4->Branch("Multiplicity", &Multiplicity);
    fTreeBS4->Branch("DiffPosition",  &DiffPosition );
    //fTreeBS4->Branch("DiffEnergy", &DiffEnergy);
    fTreeBS4->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS4->Branch("ECII", &fECII);
    fTreeBS4->SetCircular(2000000);
    
    
    fTreeS5->Branch("EventNumber", &EventNumber);
    fTreeS5->Branch("eventID", &eventIDS5);
    fTreeS5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS5->Branch("Pos_eX", &Pos_eX);
    fTreeS5->Branch("Pos_eY", &Pos_eY);
    fTreeS5->Branch("Pos_eZ", &Pos_eZ);
    fTreeS5->Branch("Pos_pX", &Pos_pX);
    fTreeS5->Branch("Pos_pY", &Pos_pY);
    fTreeS5->Branch("Pos_pZ", &Pos_pZ);
    fTreeS5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS5->Branch("PosX_Scat", &PosX_scat);
    fTreeS5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS5->Branch("PosY_Scat", &PosY_scat);
    fTreeS5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeS5->Branch("Pos_Scat", &fPosScat);
    fTreeS5->Branch("Energy_Scat", &Energy_scat);
    fTreeS5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS5->Branch("PosX_Abs", &PosX_abs);
    fTreeS5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS5->Branch("PosY_Abs", &PosY_abs);
    fTreeS5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeS5->Branch("Pos_Abs", &fPosAbs);
    fTreeS5->Branch("PosX_Sec", &PosX_sec);
    fTreeS5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeS5->Branch("PosY_Sec", &PosY_sec);
    fTreeS5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeS5->Branch("PosZ_Sec", &PosZ_sec);
    fTreeS5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeS5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeS5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeS5->Branch("EnergySecond_abs", &EnergySe);
    fTreeS5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeS5->Branch("Energy_Abs", &Energy_abs);
    fTreeS5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeS5->Branch("EnDiff", &EnDiff);
    fTreeS5->Branch("EnergySum", &EnergySum);
    fTreeS5->Branch("EnergySumUnc", &EnergySumunc);
    fTreeS5->Branch("Multiplicity", &Multiplicity);
    fTreeS5->Branch("DiffPosition",  &DiffPosition );
    //fTreeS5->Branch("DiffEnergy", &DiffEnergy);
    fTreeS5->Branch("AngularDistribution", &AngularDistribution);
    fTreeS5->Branch("ECII", &fECII);
    fTreeS5->SetCircular(2000000); 
    
    fTreeB5->Branch("EventNumber", &EventNumber);
    fTreeB5->Branch("eventID", &eventIDB5);
    fTreeB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB5->Branch("Pos_eX", &Pos_eX);
    fTreeB5->Branch("Pos_eY", &Pos_eY);
    fTreeB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeB5->Branch("Pos_pX", &Pos_pX);
    fTreeB5->Branch("Pos_pY", &Pos_pY);
    fTreeB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB5->Branch("PosX_Scat", &PosX_scat);
    fTreeB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB5->Branch("PosY_Scat", &PosY_scat);
    fTreeB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeB5->Branch("Pos_Scat", &fPosScat);
    fTreeB5->Branch("Energy_Scat", &Energy_scat);
    fTreeB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB5->Branch("PosX_Abs", &PosX_abs);
    fTreeB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB5->Branch("PosY_Abs", &PosY_abs);
    fTreeB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeB5->Branch("Pos_Abs", &fPosAbs);
    fTreeB5->Branch("PosX_Sec", &PosX_sec);
    fTreeB5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeB5->Branch("PosY_Sec", &PosY_sec);
    fTreeB5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeB5->Branch("PosZ_Sec", &PosZ_sec);
    fTreeB5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeB5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeB5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeB5->Branch("EnergySecond_abs", &EnergySe);
    fTreeB5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeB5->Branch("Energy_Abs", &Energy_abs);
    fTreeB5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeB5->Branch("EnDiff", &EnDiff);
    fTreeB5->Branch("EnergySum", &EnergySum);
    fTreeB5->Branch("EnergySumUnc", &EnergySumunc);
    fTreeB5->Branch("Multiplicity", &Multiplicity);
    fTreeB5->Branch("DiffPosition",  &DiffPosition );
    //fTreeB5->Branch("DiffEnergy", &DiffEnergy);
    fTreeB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeB5->Branch("ECII", &fECII);
    fTreeB5->SetCircular(2000000);
    
    
    fTreeBB5->Branch("EventNumber", &EventNumber);
    fTreeBB5->Branch("eventID", &eventIDBB5);
    fTreeBB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB5->Branch("Pos_eX", &Pos_eX);
    fTreeBB5->Branch("Pos_eY", &Pos_eY);
    fTreeBB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB5->Branch("Pos_pX", &Pos_pX);
    fTreeBB5->Branch("Pos_pY", &Pos_pY);
    fTreeBB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB5->Branch("PosX_Scat", &PosX_scat);
    fTreeBB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB5->Branch("PosY_Scat", &PosY_scat);
    fTreeBB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBB5->Branch("Pos_Scat", &fPosScat);
    fTreeBB5->Branch("Energy_Scat", &Energy_scat);
    fTreeBB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB5->Branch("PosX_Abs", &PosX_abs);
    fTreeBB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB5->Branch("PosY_Abs", &PosY_abs);
    fTreeBB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBB5->Branch("Pos_Abs", &fPosAbs);
    fTreeBB5->Branch("PosX_Sec", &PosX_sec);
    fTreeBB5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBB5->Branch("PosY_Sec", &PosY_sec);
    fTreeBB5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBB5->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBB5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBB5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBB5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBB5->Branch("EnergySecond_abs", &EnergySe);
    fTreeBB5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBB5->Branch("Energy_Abs", &Energy_abs);
    fTreeBB5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBB5->Branch("EnDiff", &EnDiff);
    fTreeBB5->Branch("EnergySum", &EnergySum);
    fTreeBB5->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBB5->Branch("Multiplicity", &Multiplicity);
    fTreeBB5->Branch("DiffPosition",  &DiffPosition );
    //fTreeBB5->Branch("DiffEnergy", &DiffEnergy);
    fTreeBB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB5->Branch("ECII", &fECII);
    fTreeBB5->SetCircular(2000000);
    
    fTreedB5->Branch("EventNumber", &EventNumber);
    fTreedB5->Branch("eventID", &eventIDdB5);
    fTreedB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreedB5->Branch("Pos_eX", &Pos_eX);
    fTreedB5->Branch("Pos_eY", &Pos_eY);
    fTreedB5->Branch("Pos_eZ", &Pos_eZ);
    fTreedB5->Branch("Pos_pX", &Pos_pX);
    fTreedB5->Branch("Pos_pY", &Pos_pY);
    fTreedB5->Branch("Pos_pZ", &Pos_pZ);
    fTreedB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreedB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreedB5->Branch("PosX_Scat", &PosX_scat);
    fTreedB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreedB5->Branch("PosY_Scat", &PosY_scat);
    fTreedB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreedB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreedB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreedB5->Branch("Pos_Scat", &fPosScat);
    fTreedB5->Branch("Energy_Scat", &Energy_scat);
    fTreedB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreedB5->Branch("PosX_Abs", &PosX_abs);
    fTreedB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreedB5->Branch("PosY_Abs", &PosY_abs);
    fTreedB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreedB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreedB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreedB5->Branch("Pos_Abs", &fPosAbs);
    fTreedB5->Branch("PosX_Sec", &PosX_sec);
    fTreedB5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreedB5->Branch("PosY_Sec", &PosY_sec);
    fTreedB5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreedB5->Branch("PosZ_Sec", &PosZ_sec);
    fTreedB5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreedB5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreedB5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreedB5->Branch("EnergySecond_abs", &EnergySe);
    fTreedB5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreedB5->Branch("Energy_Abs", &Energy_abs);
    fTreedB5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreedB5->Branch("EnDiff", &EnDiff);
    fTreedB5->Branch("EnergySum", &EnergySum);
    fTreedB5->Branch("EnergySumUnc", &EnergySumunc);
    fTreedB5->Branch("Multiplicity", &Multiplicity);
    fTreedB5->Branch("DiffPosition",  &DiffPosition );
    //fTreedB5->Branch("DiffEnergy", &DiffEnergy);
    fTreedB5->Branch("AngularDistribution", &AngularDistribution);
    fTreedB5->Branch("ECII", &fECII);
    fTreedB5->SetCircular(2000000);
    
    fTreeRB5->Branch("EventNumber", &EventNumber);
    fTreeRB5->Branch("eventID", &eventIDRB5);
    fTreeRB5->Branch("PrimaryEnergy", &PriEnergy);
    //fTreeRB5->Branch("PERB", &PERB);
    fTreeRB5->Branch("Pos_eX", &Pos_eX);
    fTreeRB5->Branch("Pos_eY", &Pos_eY);
    fTreeRB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB5->Branch("Pos_pX", &Pos_pX);
    fTreeRB5->Branch("Pos_pY", &Pos_pY);
    fTreeRB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB5->Branch("PosX_Scat", &PosX_scat);
    fTreeRB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB5->Branch("PosY_Scat", &PosY_scat);
    fTreeRB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeRB5->Branch("Pos_Scat", &fPosScat);
    fTreeRB5->Branch("Energy_Scat", &Energy_scat);
    fTreeRB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB5->Branch("PosX_Abs", &PosX_abs);
    fTreeRB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB5->Branch("PosY_Abs", &PosY_abs);
    fTreeRB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeRB5->Branch("Pos_Abs", &fPosAbs);
    fTreeRB5->Branch("PosX_Sec", &PosX_sec);
    fTreeRB5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeRB5->Branch("PosY_Sec", &PosY_sec);
    fTreeRB5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeRB5->Branch("PosZ_Sec", &PosZ_sec);
    fTreeRB5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeRB5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeRB5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeRB5->Branch("EnergySecond_abs", &EnergySe);
    fTreeRB5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeRB5->Branch("Energy_Abs", &Energy_abs);
    fTreeRB5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeRB5->Branch("EnDiff", &EnDiff);
    fTreeRB5->Branch("EnergySum", &EnergySum);
    fTreeRB5->Branch("EnergySumUnc", &EnergySumunc);
    fTreeRB5->Branch("Multiplicity", &Multiplicity);
    fTreeRB5->Branch("DiffPosition",  &DiffPosition );
    //fTreeRB5->Branch("DiffEnergy", &DiffEnergy);
    fTreeRB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB5->Branch("ECII", &fECII);
    fTreeRB5->SetCircular(2000000);
    
    
    fTreeBS5->Branch("EventNumber", &EventNumber);
    fTreeBS5->Branch("eventID", &eventIDBS5);
    fTreeBS5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS5->Branch("Pos_eX", &Pos_eX);
    fTreeBS5->Branch("Pos_eY", &Pos_eY);
    fTreeBS5->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS5->Branch("Pos_pX", &Pos_pX);
    fTreeBS5->Branch("Pos_pY", &Pos_pY);
    fTreeBS5->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS5->Branch("PosX_Scat", &PosX_scat);
    fTreeBS5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS5->Branch("PosY_Scat", &PosY_scat);
    fTreeBS5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeBS5->Branch("Pos_Scat", &fPosScat);
    fTreeBS5->Branch("Energy_Scat", &Energy_scat);
    fTreeBS5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS5->Branch("PosX_Abs", &PosX_abs);
    fTreeBS5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS5->Branch("PosY_Abs", &PosY_abs);
    fTreeBS5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeBS5->Branch("Pos_Abs", &fPosAbs);
    fTreeBS5->Branch("PosX_Sec", &PosX_sec);
    fTreeBS5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeBS5->Branch("PosY_Sec", &PosY_sec);
    fTreeBS5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeBS5->Branch("PosZ_Sec", &PosZ_sec);
    fTreeBS5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeBS5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeBS5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeBS5->Branch("EnergySecond_abs", &EnergySe);
    fTreeBS5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeBS5->Branch("Energy_Abs", &Energy_abs);
    fTreeBS5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeBS5->Branch("EnDiff", &EnDiff);
    fTreeBS5->Branch("EnergySum", &EnergySum);
    fTreeBS5->Branch("EnergySumUnc", &EnergySumunc);
    fTreeBS5->Branch("Multiplicity", &Multiplicity);
    fTreeBS5->Branch("DiffPosition",  &DiffPosition );
    //fTreeBS5->Branch("DiffEnergy", &DiffEnergy);
    fTreeBS5->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS5->Branch("ECII", &fECII);
    fTreeBS5->SetCircular(2000000);
///////////////////////////////////////////////////////////////////////    
        
    
    fTreeSB2->Branch("EventNumber", &EventNumber);
    fTreeSB2->Branch("eventID", &eventIDSB2);
    fTreeSB2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeSB2->Branch("Pos_eX", &Pos_eX);
    fTreeSB2->Branch("Pos_eY", &Pos_eY);
    fTreeSB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeSB2->Branch("Pos_pX", &Pos_pX);
    fTreeSB2->Branch("Pos_pY", &Pos_pY);
    fTreeSB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeSB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeSB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeSB2->Branch("PosX_Scat", &PosX_scat);
    fTreeSB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeSB2->Branch("PosY_Scat", &PosY_scat);
    fTreeSB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeSB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeSB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeSB2->Branch("Pos_Scat", &PosScat);
    fTreeSB2->Branch("Energy_Scat", &Energy_scat);
    fTreeSB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeSB2->Branch("PosX_Abs", &PosX_abs);
    fTreeSB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeSB2->Branch("PosY_Abs", &PosY_abs);
    fTreeSB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeSB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeSB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeSB2->Branch("Pos_Abs", &PosAbs);
    fTreeSB2->Branch("PosX_Sec", &PosX_sec);
    fTreeSB2->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeSB2->Branch("PosY_Sec", &PosY_sec);
    fTreeSB2->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeSB2->Branch("PosZ_Sec", &PosZ_sec);
    fTreeSB2->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeSB2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeSB2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeSB2->Branch("EnergySecond_abs", &EnergySe);
    fTreeSB2->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeSB2->Branch("Energy_Abs", &Energy_abs);
    fTreeSB2->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeSB2->Branch("EnDiff", &EnDiff);
    fTreeSB2->Branch("EnergySum", &EnergySum);
    fTreeSB2->Branch("EnergySumUnc", &EnergySumunc);
    fTreeSB2->Branch("Multiplicity", &Multiplicity);
    fTreeSB2->Branch("DiffPosition", &DiffPosition);
    //fTreeSB2->Branch("DiffEnergy", &DiffEnergy);
    fTreeSB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeSB2->Branch("ECII", &fECII);
    fTreeSB2->SetCircular(2000000); 
    
    
    fTreeSB3->Branch("EventNumber", &EventNumber);
    fTreeSB3->Branch("eventID", &eventIDSB3);
    fTreeSB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeSB3->Branch("Pos_eX", &Pos_eX);
    fTreeSB3->Branch("Pos_eY", &Pos_eY);
    fTreeSB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeSB3->Branch("Pos_pX", &Pos_pX);
    fTreeSB3->Branch("Pos_pY", &Pos_pY);
    fTreeSB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeSB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeSB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeSB3->Branch("PosX_Scat", &PosX_scat);
    fTreeSB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeSB3->Branch("PosY_Scat", &PosY_scat);
    fTreeSB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeSB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeSB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeSB3->Branch("Pos_Scat", &PosScat);
    fTreeSB3->Branch("Energy_Scat", &Energy_scat);
    fTreeSB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeSB3->Branch("PosX_Abs", &PosX_abs);
    fTreeSB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeSB3->Branch("PosY_Abs", &PosY_abs);
    fTreeSB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeSB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeSB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeSB3->Branch("Pos_Abs", &PosAbs);
    fTreeSB3->Branch("PosX_Sec", &PosX_sec);
    fTreeSB3->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeSB3->Branch("PosY_Sec", &PosY_sec);
    fTreeSB3->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeSB3->Branch("PosZ_Sec", &PosZ_sec);
    fTreeSB3->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeSB3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeSB3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeSB3->Branch("EnergySecond_abs", &EnergySe);
    fTreeSB3->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeSB3->Branch("Energy_Abs", &Energy_abs);
    fTreeSB3->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeSB3->Branch("EnDiff", &EnDiff);
    fTreeSB3->Branch("EnergySum", &EnergySum);
    fTreeSB3->Branch("EnergySumUnc", &EnergySumunc);
    fTreeSB3->Branch("Multiplicity", &Multiplicity);
    fTreeSB3->Branch("DiffPosition", &DiffPosition);
    //fTreeSB3->Branch("DiffEnergy", &DiffEnergy);
    fTreeSB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeSB3->Branch("ECII", &fECII);
    fTreeSB3->SetCircular(2000000); 
    
    fTreeSB4->Branch("EventNumber", &EventNumber);
    fTreeSB4->Branch("eventID", &eventIDSB4);
    fTreeSB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeSB4->Branch("Pos_eX", &Pos_eX);
    fTreeSB4->Branch("Pos_eY", &Pos_eY);
    fTreeSB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeSB4->Branch("Pos_pX", &Pos_pX);
    fTreeSB4->Branch("Pos_pY", &Pos_pY);
    fTreeSB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeSB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeSB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeSB4->Branch("PosX_Scat", &PosX_scat);
    fTreeSB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeSB4->Branch("PosY_Scat", &PosY_scat);
    fTreeSB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeSB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeSB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeSB4->Branch("Pos_Scat", &PosScat);
    fTreeSB4->Branch("Energy_Scat", &Energy_scat);
    fTreeSB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeSB4->Branch("PosX_Abs", &PosX_abs);
    fTreeSB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeSB4->Branch("PosY_Abs", &PosY_abs);
    fTreeSB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeSB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeSB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeSB4->Branch("Pos_Abs", &PosAbs);
    fTreeSB4->Branch("PosX_Sec", &PosX_sec);
    fTreeSB4->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeSB4->Branch("PosY_Sec", &PosY_sec);
    fTreeSB4->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeSB4->Branch("PosZ_Sec", &PosZ_sec);
    fTreeSB4->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeSB4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeSB4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeSB4->Branch("EnergySecond_abs", &EnergySe);
    fTreeSB4->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeSB4->Branch("Energy_Abs", &Energy_abs);
    fTreeSB4->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeSB4->Branch("EnDiff", &EnDiff);
    fTreeSB4->Branch("EnergySum", &EnergySum);
    fTreeSB4->Branch("EnergySumUnc", &EnergySumunc);
    fTreeSB4->Branch("Multiplicity", &Multiplicity);
    fTreeSB4->Branch("DiffPosition", &DiffPosition);
    //fTreeSB4->Branch("DiffEnergy", &DiffEnergy);
    fTreeSB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeSB4->Branch("ECII", &fECII);
    fTreeSB4->SetCircular(2000000); 
    
    fTreeSB5->Branch("EventNumber", &EventNumber);
    fTreeSB5->Branch("eventID", &eventIDSB5);
    fTreeSB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeSB5->Branch("Pos_eX", &Pos_eX);
    fTreeSB5->Branch("Pos_eY", &Pos_eY);
    fTreeSB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeSB5->Branch("Pos_pX", &Pos_pX);
    fTreeSB5->Branch("Pos_pY", &Pos_pY);
    fTreeSB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeSB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeSB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeSB5->Branch("PosX_Scat", &PosX_scat);
    fTreeSB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeSB5->Branch("PosY_Scat", &PosY_scat);
    fTreeSB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeSB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeSB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeSB5->Branch("Pos_Scat", &PosScat);
    fTreeSB5->Branch("Energy_Scat", &Energy_scat);
    fTreeSB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeSB5->Branch("PosX_Abs", &PosX_abs);
    fTreeSB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeSB5->Branch("PosY_Abs", &PosY_abs);
    fTreeSB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeSB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeSB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeSB5->Branch("Pos_Abs", &PosAbs);
    fTreeSB5->Branch("PosX_Sec", &PosX_sec);
    fTreeSB5->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeSB5->Branch("PosY_Sec", &PosY_sec);
    fTreeSB5->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeSB5->Branch("PosZ_Sec", &PosZ_sec);
    fTreeSB5->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeSB5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeSB5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeSB5->Branch("EnergySecond_abs", &EnergySe);
    fTreeSB5->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeSB5->Branch("Energy_Abs", &Energy_abs);
    fTreeSB5->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeSB5->Branch("EnDiff", &EnDiff);
    fTreeSB5->Branch("EnergySum", &EnergySum);
    fTreeSB5->Branch("EnergySumUnc", &EnergySumunc);
    fTreeSB5->Branch("Multiplicity", &Multiplicity);
    fTreeSB5->Branch("DiffPosition", &DiffPosition);
    //fTreeSB5->Branch("DiffEnergy", &DiffEnergy);
    fTreeSB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeSB5->Branch("ECII", &fECII);
    fTreeSB5->SetCircular(2000000);
    
    fTreeSB->Branch("EventNumber", &EventNumber);
    fTreeSB->Branch("eventID", &eventIDSB5);
    fTreeSB->Branch("PrimaryEnergy", &PriEnergy);
    fTreeSB->Branch("Pos_eX", &Pos_eX);
    fTreeSB->Branch("Pos_eY", &Pos_eY);
    fTreeSB->Branch("Pos_eZ", &Pos_eZ);
    fTreeSB->Branch("Pos_pX", &Pos_pX);
    fTreeSB->Branch("Pos_pY", &Pos_pY);
    fTreeSB->Branch("Pos_pZ", &Pos_pZ);
    fTreeSB->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeSB->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeSB->Branch("PosX_Scat", &PosX_scat);
    fTreeSB->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeSB->Branch("PosY_Scat", &PosY_scat);
    fTreeSB->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeSB->Branch("PosZ_Scat", &PosZ_scat);
    fTreeSB->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    //fTreeSB5->Branch("Pos_Scat", &PosScat);
    fTreeSB->Branch("Energy_Scat", &Energy_scat);
    fTreeSB->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeSB->Branch("PosX_Abs", &PosX_abs);
    fTreeSB->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeSB->Branch("PosY_Abs", &PosY_abs);
    fTreeSB->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeSB->Branch("PosZ_Abs", &PosZ_abs);
    fTreeSB->Branch("PosZ_AbsUnc", &PosZ_absunc);
    //fTreeSB5->Branch("Pos_Abs", &PosAbs);
    fTreeSB->Branch("PosX_Sec", &PosX_sec);
    fTreeSB->Branch("PosX_SecUnc", &PosX_secunc);
    fTreeSB->Branch("PosY_Sec", &PosY_sec);
    fTreeSB->Branch("PosY_SecUnc", &PosY_secunc);
    fTreeSB->Branch("PosZ_Sec", &PosZ_sec);
    fTreeSB->Branch("PosZ_SecUnc", &PosZ_secunc);
    fTreeSB->Branch("EnergyCluster_abs", &EnergyCluster_abs);
    fTreeSB->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
    fTreeSB->Branch("EnergySecond_abs", &EnergySe);
    fTreeSB->Branch("EnergySecond_absUnc", &EnergySeunc);
    fTreeSB->Branch("Energy_Abs", &Energy_abs);
    fTreeSB->Branch("Energy_AbsUnc", &Energy_absunc);
    fTreeSB->Branch("EnDiff", &EnDiff);
    fTreeSB->Branch("EnergySum", &EnergySum);
    fTreeSB->Branch("EnergySumUnc", &EnergySumunc);
    fTreeSB->Branch("Multiplicity", &Multiplicity);
    fTreeSB->Branch("DiffPosition", &DiffPosition);
    //fTreeSB5->Branch("DiffEnergy", &DiffEnergy);
    fTreeSB->Branch("AngularDistribution", &AngularDistribution);
    fTreeSB->Branch("ECII", &fECII);
    fTreeSB->SetCircular(2000000); 
    
      // load the input file
    TFile *input(0);
    TString fname = "./PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5-full-sigma2.root";
    input = TFile::Open( fname );
    
    TTree *t = NULL;
    input->GetObject("TreeStat", t);
    
    TString filename= "PMMA180MeV0mmBPType3_EnoughStatistics.root";
    Int_t istart, istop, totalsimnev; 
    
    //t->SetBranchAddress("InputFilename", &filename);
    t->SetBranchAddress("StartEvent", &istart);
    t->SetBranchAddress("StopEvent", &istop);
    t->SetBranchAddress("TotalSimNev", &totalsimnev);
    
    //TString filename = input->GetName();
    Int_t istart1, istop1, totalsimnev1; 
    
    fTree_stat->Branch("InputFilename", &filename); 
    fTree_stat->Branch("StartEvent", &istart1); 
    fTree_stat->Branch("StopEvent", &istop1); 
    fTree_stat->Branch("TotalSimNev", &totalsimnev1); 
    
    const ULong64_t nEntries = t->GetEntries();
    //cout << nEntries << endl <<endl;
    for (ULong64_t e = 0; e < nEntries; ++e) {
        
        t->GetEntry(e);
        
        istart1 = totalsimnev/2;
        istop1 = istop;
        totalsimnev1 = totalsimnev/2;
        
        fTree_stat->Fill();
    }
    
    
    
    
    TTree* theTree = NULL;
    
    // loop through signal and all background trees
    for( int treeNumber = 0; treeNumber < 28; ++treeNumber ) {
        
        TStopwatch sw;
        sw.Start();
       
        if( treeNumber == 0 ){
            
            
            theTree = (TTree*)input->Get("TreeS2");
            std::cout << "--- Select signal sample" << std::endl;
            
             theTree->SetBranchAddress("EventNumber", &EventNumber);
             theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
             theTree->SetBranchAddress("Pos_eX", &Pos_eX);
             theTree->SetBranchAddress("Pos_eY", &Pos_eY);
             theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
             theTree->SetBranchAddress("Pos_pX", &Pos_pX);
             theTree->SetBranchAddress("Pos_pY", &Pos_pY);
            theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
            theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
            theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
            theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
            theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
            theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
            theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
            theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
             theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
             theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
            theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
            theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
            theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
            theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
            theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
            theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
            theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
            theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
            theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
            theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
            theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
            theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
            theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
            theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
            theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
            theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
            theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
            theTree->SetBranchAddress( "EnDiff", &EnDiff );
            theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
            theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
            theTree->SetBranchAddress( "EnergySum", &EnergySum );
            theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
            theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
            theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
            theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
             theTree->SetBranchAddress( "ECII", &fECII );
             theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
            
            std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
            Int_t nEvent = theTree->GetEntries();
            
            for (Long64_t ievt=4767; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDS2++;
               //AngularDistribution = 0;
               
               fTreeS2->Fill();
            }
            
        } else if( treeNumber == 1 ){
            
           
           theTree = (TTree*)input->Get("TreeB2");
           std::cout << "--- Select background sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           //Int_t nEvent = 100;
           for (Long64_t ievt=85684; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDB2++;

               //AngularDistribution = 0;
               fTreeB2->Fill();
           }
           
           
       } else if( treeNumber == 2 ){
            
           
           theTree = (TTree*)input->Get("TreeBB2");
           std::cout << "--- Select background sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=32117; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBB2++;
               //AngularDistribution = 0;

               fTreeBB2->Fill();
           }
           
           
       } else if( treeNumber == 3 ){
      
           theTree = (TTree*)input->Get("TreeRB2");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=50894; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDRB2++;
               //AngularDistribution = 0;
               fTreeRB2->Fill();
           }
           
           
           
       }else if( treeNumber == 4 ){
      
           theTree = (TTree*)input->Get("TreeBS2");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=2673; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBS2++;
 
               //AngularDistribution = 0;
               fTreeBS2->Fill();
           }
           
           
           
       } else if( treeNumber == 5 ){
           
           theTree = (TTree*)input->Get("TreeS3");
           std::cout << "--- Select signal sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=5768; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDS3++;
               
               fTreeS3->Fill();
           }
           
           
       } else if( treeNumber == 6 ){
           
           theTree = (TTree*)input->Get("TreeB3");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=75814; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDB3++;
               
               fTreeB3->Fill();
           }
           
           
       } else if( treeNumber == 7 ){
           
           theTree = (TTree*)input->Get("TreeBB3");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=13093; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBB3++;
               
               fTreeBB3->Fill();
           }
           
           
       } else if( treeNumber == 8 ){
           
           theTree = (TTree*)input->Get("TreedB3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=18861; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDdB3++;
               
               fTreedB3->Fill();
           }
           
           
       } else if( treeNumber == 9 ){
           
           theTree = (TTree*)input->Get("TreeRB3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=43184; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDRB3++;
               
               fTreeRB3->Fill();
           }
           
           
       } else if( treeNumber == 10 ){
           
           theTree = (TTree*)input->Get("TreeBS3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=676; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBS3++;
               
               fTreeBS3->Fill();
           }
           
           
       } else if( treeNumber == 11 ){
           
           theTree = (TTree*)input->Get("TreeS4");
           std::cout << "--- Select Real signal sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=2917; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               eventIDS4++;
               fTreeS4->Fill();
           }
           
       } else if( treeNumber == 12 ){
           
           theTree = (TTree*)input->Get("TreeB4");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=53540; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDB4++;
               fTreeB4->Fill();
           }
           
       } else if( treeNumber == 13 ){
           
           theTree = (TTree*)input->Get("TreeBB4");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=4481; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBB4++;
               
               fTreeBB4->Fill();
           }
           
       } else if( treeNumber == 14 ){
           
           theTree = (TTree*)input->Get("TreedB4");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=14796; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
                eventIDdB4++;
               fTreedB4->Fill();
           }
           
       } else if( treeNumber == 15 ){
           
           theTree = (TTree*)input->Get("TreeRB4");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
         
           for (Long64_t ievt=33918; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDRB4++;
               fTreeRB4->Fill();
           }
           
       } else if( treeNumber == 16 ){
           
           theTree = (TTree*)input->Get("TreeBS4");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=345; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
                eventIDBS4++;
               fTreeBS4->Fill();
           }
           
       } else if( treeNumber == 17 ){
           
           theTree = (TTree*)input->Get("TreeS5");
           std::cout << "--- Select Real signal sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=954; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              eventIDS5++;
               fTreeS5->Fill();
           }
           
       } else if( treeNumber == 18 ){
           
           theTree = (TTree*)input->Get("TreeB5");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=30514; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDB5++;
               fTreeB5->Fill();
           }
           
       } else if( treeNumber == 19 ){
           
           theTree = (TTree*)input->Get("TreeBB5");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=1273; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBB5++;
               fTreeBB5->Fill();
           }
           
       } else if( treeNumber == 20 ){
           
           theTree = (TTree*)input->Get("TreedB5");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=6681; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDdB5++;
               fTreedB5->Fill();
           }
           
       } else if( treeNumber == 21 ){
           
           theTree = (TTree*)input->Get("TreeRB5");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
         
           for (Long64_t ievt=22412; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
                eventIDRB5++;
               fTreeRB5->Fill();
           }
           
       } else if( treeNumber == 22 ){
           
           theTree = (TTree*)input->Get("TreeBS5");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "ECII", &fECII );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=148; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDBS5++;
               
               fTreeBS5->Fill();
           }
           
       } else if( treeNumber == 23 ){
           
           theTree = (TTree*)input->Get("TreeSB2");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=90451; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDSB2++;
               
               fTreeSB2->Fill();
           }
           
       } else if( treeNumber == 24 ){
           
           theTree = (TTree*)input->Get("TreeSB3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=81582; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDSB3++;
               
               fTreeSB3->Fill();
           }
           
       } else if( treeNumber == 25 ){
           
           theTree = (TTree*)input->Get("TreeSB4");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=56457; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDSB4++;
               
               fTreeSB4->Fill();
           }
           
       } else if( treeNumber == 26 ){
           
           theTree = (TTree*)input->Get("TreeSB5");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=31468; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDSB5++;
               
               fTreeSB5->Fill();
           }
           
       } else if( treeNumber == 27 ){
           
           theTree = (TTree*)input->Get("TreeSB");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
           theTree->SetBranchAddress("Pos_eX", &Pos_eX);
           theTree->SetBranchAddress("Pos_eY", &Pos_eY);
           theTree->SetBranchAddress("Pos_eZ", &Pos_eZ);
           theTree->SetBranchAddress("Pos_pX", &Pos_pX);
           theTree->SetBranchAddress("Pos_pY", &Pos_pY);
           theTree->SetBranchAddress("Pos_pZ", &Pos_pZ);
           theTree->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
           theTree->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
           theTree->SetBranchAddress( "PosX_Scat", &PosX_scat );
           theTree->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
           theTree->SetBranchAddress( "PosY_Scat", &PosY_scat );
           theTree->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
           theTree->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
           theTree->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
           theTree->SetBranchAddress( "Energy_Scat", &Energy_scat );
           theTree->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
           theTree->SetBranchAddress( "PosX_Abs", &PosX_abs );
           theTree->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
           theTree->SetBranchAddress( "PosY_Abs", &PosY_abs );
           theTree->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
           theTree->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
           theTree->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
           theTree->SetBranchAddress( "PosX_Sec", &PosX_sec );
           theTree->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
           theTree->SetBranchAddress( "PosY_Sec", &PosY_sec );
           theTree->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
           theTree->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
           theTree->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=259958; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               eventIDSB++;
               
               fTreeSB->Fill();
           }
           
       }
        
         // get elapsed time
        sw.Stop();
        std::cout << "--- End of event loop: "; sw.Print();
    }
    
    input->Close();
   // write output tree
/*   outputTree->SetDirectory(outputFile);
     outputTree->Write(); */
    outputFile->Write();
    outputFile->Close();
    std::cout << "--- Created root file: \"" << outfileName.Data() << "\" containing event classes trees" << std::endl;
           
            
}
    
