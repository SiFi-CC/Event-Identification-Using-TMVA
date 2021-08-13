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

#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"//required to load dataset
#include "TMVA/Reader.h"

#include "TRandom.h"
#include <algorithm>

#include <cstdlib>
#include <map>

#include "TTree.h"
#include "TString.h"

#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "TVector3.h"

#include "TMVA/TMVARegGui.h"

using namespace std;
using namespace TMVA;


// -----------------------------------------------------------------------------------------
// Regression of Predicted Compton events' energies
// -----------------------------------------------------------------------------------------
//
/// For ROOT version 6.18.04 and higher ///
/// The energy regression is done for the first half statistics.  
/// The deposited energy in the scatterer and the difference energy between the primary energy 
/// and the deposited energy in the scatterer are used as two variables in recovering total 
/// energy sum in the detector because of energy leaks in the energysum variable.   
/// This energy regression is done ONLY for true Compton events 
/// because their energy sum have a really good linear relationship with the primary energy
void Regression() {
    

    
     // Create a new root output file
   TString outfileName( "TMVARESReg.root" );
   
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
/*   
   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=N;D;P;G,D:AnalysisType=Regression" );
   */
   TMVA::Factory *factory = new TMVA::Factory( "TMVARegression", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
   
   TMVA::DataLoader *dataloader2=new TMVA::DataLoader("dataset2res");

   TMVA::DataLoader *dataloader3=new TMVA::DataLoader("dataset3res");
 
   TMVA::DataLoader *dataloader4=new TMVA::DataLoader("dataset4res");

   TMVA::DataLoader *dataloader5=new TMVA::DataLoader("dataset5res");

   //dataloader2->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader2->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   //dataloader2->AddVariable( "EnergySum","EnergySum", 'F', 0, 20 );
   dataloader2->AddVariable( "EnDiff","EnDiff", 'F', -20, 20 );      

   //dataloader3->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader3->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   //dataloader3->AddVariable( "EnergySum","EnergySum", 'F', 0, 20 );
   dataloader3->AddVariable( "EnDiff","EnDiff", 'F', -20, 20 );

   //dataloader4->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader4->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   //dataloader4->AddVariable( "EnergySum","EnergySum", 'F', 0, 20 );
   dataloader4->AddVariable( "EnDiff","EnDiff", 'F', -20, 20 );

   //dataloader5->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader5->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   //dataloader5->AddVariable( "EnergySum","EnergySum", 'F', 0, 20 );
   dataloader5->AddVariable( "EnDiff","EnDiff", 'F', -20, 20 );
   
   //dataloader2->AddTarget( "EnergySum" );
   dataloader2->AddTarget( "PrimaryEnergy" );

   //dataloader3->AddTarget( "EnergySum" ); 
   dataloader3->AddTarget( "PrimaryEnergy" );

   //dataloader4->AddTarget( "EnergySum" ); 
   dataloader4->AddTarget( "PrimaryEnergy" );

   //dataloader5->AddTarget( "EnergySum" );
   dataloader5->AddTarget( "PrimaryEnergy" );
 

// load the signal and background event samples from ROOT trees
   
   TString fname = "./PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5FirsthalfModified2sigma.root";
   
   TFile *input(0);
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); 
   }
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   
   std::cout << "--- TMVARegression           : Using input file: " << input->GetName() << std::endl;

   // Register the regression tree

   TTree *signal2      = (TTree*)input->Get("TreeS2");

   TTree *signal3      = (TTree*)input->Get("TreeS3");

   TTree *signal4      = (TTree*)input->Get("TreeS4");

   TTree *signal5      = (TTree*)input->Get("TreeS5");

   // global event weights per tree 
   Double_t regWeight  = 1.0;

// You can add an arbitrary number of regression trees
   dataloader2->AddRegressionTree( signal2, regWeight );

   dataloader3->AddRegressionTree( signal3, regWeight );

   dataloader4->AddRegressionTree( signal4, regWeight );

   dataloader5->AddRegressionTree( signal5, regWeight );
   
   // This would set individual event weights (the variables defined in the
   // expression need to exist in the original TTree)
/*   
   dataloader2->SetWeightExpression( "eventID", "Regression" );
   dataloader3->SetWeightExpression( "eventID", "Regression" );
   dataloader4->SetWeightExpression( "eventID", "Regression" );
   dataloader5->SetWeightExpression( "eventID", "Regression" );
  */ 
    // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = ""; 
   
   // tell the DataLoader to use all remaining events in the trees after training for testing:
   
   dataloader2->PrepareTrainingAndTestTree( mycut,
                                         "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=EqualNumEvents:!V" );
  
   dataloader3->PrepareTrainingAndTestTree( mycut,
                                         "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=EqualNumEvents:!V" );

   dataloader4->PrepareTrainingAndTestTree( mycut,
                                         "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=EqualNumEvents:!V" );
  
   dataloader5->PrepareTrainingAndTestTree( mycut,
                                         "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=EqualNumEvents:!V" );
  
///////////////////////////////////////////////////////////////// Gradient Boosted (BDTG) ///////////UseBaggedBoost/////////////////////////////////////

   factory->BookMethod( dataloader2, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=8000:MinNodeSize=0.1%:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=4" );

   factory->BookMethod( dataloader3, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=10000:MinNodeSize=0.2%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=4" );

   factory->BookMethod( dataloader4, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=8000:MinNodeSize=0.2%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=4" ); 

   factory->BookMethod( dataloader5, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=8000:MinNodeSize=0.1%:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=4" );

   // Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVARegression is done!" << std::endl;

   delete factory;
   
   delete dataloader2;
 
   delete dataloader3;

   delete dataloader4;

   delete dataloader5;

   // Launch the GUI for the root macros
   //if (!gROOT->IsBatch()) TMVA::TMVARegGui( outfileName );
}
// -----------------------------------------------------------------------------------------
// Regression Application
// -----------------------------------------------------------------------------------------
//
void RegressionApplied(){
   
   // Create a new root output file.
   TString outfileName( "PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5-FHRegression.root" );
   
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
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
    
    Float_t Pos_eX, Pos_eY, Pos_eZ, Pos_pX, Pos_pY, Pos_pZ, RealEnergy_e, RealEnergy_p; 
   
    Float_t PosX_scat, PosY_scat, PosZ_scat,Energy_scat,ReEnergy_scat;
   
    Float_t PosX_abs, PosY_abs, PosZ_abs, Energy_abs, ReEnergy_abs;
   
    Float_t PosX_sec, PosY_sec, PosZ_sec;
    Float_t EnergyCluster_abs, AngularDistribution;
   
    Float_t PosX_scatunc, PosY_scatunc, PosZ_scatunc,Energy_scatunc;
    Float_t PosX_absunc, PosY_absunc, PosZ_absunc, Energy_absunc, EnDiff;
    Float_t PosX_secunc, PosY_secunc, PosZ_secunc;
   
    Float_t EnergyCluster_absunc,EnergySe, EnergySeunc;
   
    Float_t DiffEnergy, RatioEnergy, DiffPosition, PriEnergy, EnergySum, EnergySumunc, PERB;
    Int_t count = 0;
    Int_t EventNumber;
    Int_t Multiplicity;
    Float_t PosX_1, PosY_1, PosZ_1, Energy_1;
    Float_t PosX_2, PosY_2, PosZ_2, Energy_2;
    Float_t DEnergy, REnergy, DPosition, PEnergy, EnergyS, ReES;
    Int_t eventID;

    
    fTreeS2->Branch("EventNumber", &EventNumber);
    fTreeS2->Branch("eventID", &eventID);
    fTreeS2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS2->Branch("ReEnergySum", &ReES);
    fTreeS2->Branch("Pos_eX", &Pos_eX);
    fTreeS2->Branch("Pos_eY", &Pos_eY);
    fTreeS2->Branch("Pos_eZ", &Pos_eZ);
    fTreeS2->Branch("Pos_pX", &Pos_pX);
    fTreeS2->Branch("Pos_pY", &Pos_pY);
    fTreeS2->Branch("Pos_pZ", &Pos_pZ);
    fTreeS2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS2->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeS2->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeS2->Branch("PosX_Scat", &PosX_scat);
    fTreeS2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS2->Branch("PosY_Scat", &PosY_scat);
    fTreeS2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeS2->Branch("Energy_Scat", &Energy_scat);
    fTreeS2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS2->Branch("PosX_Abs", &PosX_abs);
    fTreeS2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS2->Branch("PosY_Abs", &PosY_abs);
    fTreeS2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS2->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeS2->Branch("AngularDistribution", &AngularDistribution);
    fTreeS2->SetCircular(2000000); 
    
    fTreeB2->Branch("EventNumber", &EventNumber);
    fTreeB2->Branch("eventID", &eventID);
    fTreeB2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB2->Branch("ReEnergySum", &ReES);
    fTreeB2->Branch("Pos_eX", &Pos_eX);
    fTreeB2->Branch("Pos_eY", &Pos_eY);
    fTreeB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeB2->Branch("Pos_pX", &Pos_pX);
    fTreeB2->Branch("Pos_pY", &Pos_pY);
    fTreeB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB2->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeB2->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeB2->Branch("PosX_Scat", &PosX_scat);
    fTreeB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB2->Branch("PosY_Scat", &PosY_scat);
    fTreeB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeB2->Branch("Energy_Scat", &Energy_scat);
    fTreeB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB2->Branch("PosX_Abs", &PosX_abs);
    fTreeB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB2->Branch("PosY_Abs", &PosY_abs);
    fTreeB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeB2->SetCircular(2000000);
    
    fTreeBB2->Branch("EventNumber", &EventNumber);
    fTreeBB2->Branch("eventID", &eventID);
    fTreeBB2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB2->Branch("ReEnergySum", &ReES);
    fTreeBB2->Branch("Pos_eX", &Pos_eX);
    fTreeBB2->Branch("Pos_eY", &Pos_eY);
    fTreeBB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB2->Branch("Pos_pX", &Pos_pX);
    fTreeBB2->Branch("Pos_pY", &Pos_pY);
    fTreeBB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB2->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBB2->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBB2->Branch("PosX_Scat", &PosX_scat);
    fTreeBB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB2->Branch("PosY_Scat", &PosY_scat);
    fTreeBB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBB2->Branch("Energy_Scat", &Energy_scat);
    fTreeBB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB2->Branch("PosX_Abs", &PosX_abs);
    fTreeBB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB2->Branch("PosY_Abs", &PosY_abs);
    fTreeBB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB2->SetCircular(2000000); 
    
    fTreeRB2->Branch("EventNumber", &EventNumber);
    fTreeRB2->Branch("eventID", &eventID);
    fTreeRB2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeRB2->Branch("ReEnergySum", &ReES);
    fTreeRB2->Branch("Pos_eX", &Pos_eX);
    fTreeRB2->Branch("Pos_eY", &Pos_eY);
    fTreeRB2->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB2->Branch("Pos_pX", &Pos_pX);
    fTreeRB2->Branch("Pos_pY", &Pos_pY);
    fTreeRB2->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB2->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeRB2->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeRB2->Branch("PosX_Scat", &PosX_scat);
    fTreeRB2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB2->Branch("PosY_Scat", &PosY_scat);
    fTreeRB2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeRB2->Branch("Energy_Scat", &Energy_scat);
    fTreeRB2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB2->Branch("PosX_Abs", &PosX_abs);
    fTreeRB2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB2->Branch("PosY_Abs", &PosY_abs);
    fTreeRB2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB2->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeRB2->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB2->SetCircular(2000000);
    
    fTreeBS2->Branch("EventNumber", &EventNumber);
    fTreeBS2->Branch("eventID", &eventID);
    fTreeBS2->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS2->Branch("ReEnergySum", &ReES);
    fTreeBS2->Branch("Pos_eX", &Pos_eX);
    fTreeBS2->Branch("Pos_eY", &Pos_eY);
    fTreeBS2->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS2->Branch("Pos_pX", &Pos_pX);
    fTreeBS2->Branch("Pos_pY", &Pos_pY);
    fTreeBS2->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS2->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS2->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS2->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBS2->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBS2->Branch("PosX_Scat", &PosX_scat);
    fTreeBS2->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS2->Branch("PosY_Scat", &PosY_scat);
    fTreeBS2->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS2->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBS2->Branch("Energy_Scat", &Energy_scat);
    fTreeBS2->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS2->Branch("PosX_Abs", &PosX_abs);
    fTreeBS2->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS2->Branch("PosY_Abs", &PosY_abs);
    fTreeBS2->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS2->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS2->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBS2->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS2->SetCircular(2000000);
    
    
    fTreeS3->Branch("EventNumber", &EventNumber);
    fTreeS3->Branch("eventID", &eventID);
    fTreeS3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS3->Branch("ReEnergySum", &ReES);
    fTreeS3->Branch("Pos_eX", &Pos_eX);
    fTreeS3->Branch("Pos_eY", &Pos_eY);
    fTreeS3->Branch("Pos_eZ", &Pos_eZ);
    fTreeS3->Branch("Pos_pX", &Pos_pX);
    fTreeS3->Branch("Pos_pY", &Pos_pY);
    fTreeS3->Branch("Pos_pZ", &Pos_pZ);
    fTreeS3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS3->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeS3->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeS3->Branch("PosX_Scat", &PosX_scat);
    fTreeS3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS3->Branch("PosY_Scat", &PosY_scat);
    fTreeS3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeS3->Branch("Energy_Scat", &Energy_scat);
    fTreeS3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS3->Branch("PosX_Abs", &PosX_abs);
    fTreeS3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS3->Branch("PosY_Abs", &PosY_abs);
    fTreeS3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS3->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeS3->Branch("AngularDistribution", &AngularDistribution);
    fTreeS3->SetCircular(2000000); 
    
    fTreeB3->Branch("EventNumber", &EventNumber);
    fTreeB3->Branch("eventID", &eventID);
    fTreeB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB3->Branch("ReEnergySum", &ReES);
    fTreeB3->Branch("Pos_eX", &Pos_eX);
    fTreeB3->Branch("Pos_eY", &Pos_eY);
    fTreeB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeB3->Branch("Pos_pX", &Pos_pX);
    fTreeB3->Branch("Pos_pY", &Pos_pY);
    fTreeB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB3->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeB3->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeB3->Branch("PosX_Scat", &PosX_scat);
    fTreeB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB3->Branch("PosY_Scat", &PosY_scat);
    fTreeB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeB3->Branch("Energy_Scat", &Energy_scat);
    fTreeB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB3->Branch("PosX_Abs", &PosX_abs);
    fTreeB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB3->Branch("PosY_Abs", &PosY_abs);
    fTreeB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeB3->SetCircular(2000000);
    
    fTreeBB3->Branch("EventNumber", &EventNumber);
    fTreeBB3->Branch("eventID", &eventID);
    fTreeBB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB3->Branch("ReEnergySum", &ReES);
    fTreeBB3->Branch("Pos_eX", &Pos_eX);
    fTreeBB3->Branch("Pos_eY", &Pos_eY);
    fTreeBB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB3->Branch("Pos_pX", &Pos_pX);
    fTreeBB3->Branch("Pos_pY", &Pos_pY);
    fTreeBB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB3->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBB3->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBB3->Branch("PosX_Scat", &PosX_scat);
    fTreeBB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB3->Branch("PosY_Scat", &PosY_scat);
    fTreeBB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBB3->Branch("Energy_Scat", &Energy_scat);
    fTreeBB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB3->Branch("PosX_Abs", &PosX_abs);
    fTreeBB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB3->Branch("PosY_Abs", &PosY_abs);
    fTreeBB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB3->SetCircular(2000000);
    
    fTreedB3->Branch("EventNumber", &EventNumber);
    fTreedB3->Branch("eventID", &eventID);
    fTreedB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreedB3->Branch("ReEnergySum", &ReES);
    fTreedB3->Branch("Pos_eX", &Pos_eX);
    fTreedB3->Branch("Pos_eY", &Pos_eY);
    fTreedB3->Branch("Pos_eZ", &Pos_eZ);
    fTreedB3->Branch("Pos_pX", &Pos_pX);
    fTreedB3->Branch("Pos_pY", &Pos_pY);
    fTreedB3->Branch("Pos_pZ", &Pos_pZ);
    fTreedB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreedB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreedB3->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreedB3->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreedB3->Branch("PosX_Scat", &PosX_scat);
    fTreedB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreedB3->Branch("PosY_Scat", &PosY_scat);
    fTreedB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreedB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreedB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreedB3->Branch("Energy_Scat", &Energy_scat);
    fTreedB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreedB3->Branch("PosX_Abs", &PosX_abs);
    fTreedB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreedB3->Branch("PosY_Abs", &PosY_abs);
    fTreedB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreedB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreedB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreedB3->Branch("AngularDistribution", &AngularDistribution);
    fTreedB3->SetCircular(2000000);
    
    fTreeRB3->Branch("EventNumber", &EventNumber);
    fTreeRB3->Branch("eventID", &eventID);
    fTreeRB3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeRB3->Branch("ReEnergySum", &ReES);
    fTreeRB3->Branch("Pos_eX", &Pos_eX);
    fTreeRB3->Branch("Pos_eY", &Pos_eY);
    fTreeRB3->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB3->Branch("Pos_pX", &Pos_pX);
    fTreeRB3->Branch("Pos_pY", &Pos_pY);
    fTreeRB3->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB3->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeRB3->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeRB3->Branch("PosX_Scat", &PosX_scat);
    fTreeRB3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB3->Branch("PosY_Scat", &PosY_scat);
    fTreeRB3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeRB3->Branch("Energy_Scat", &Energy_scat);
    fTreeRB3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB3->Branch("PosX_Abs", &PosX_abs);
    fTreeRB3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB3->Branch("PosY_Abs", &PosY_abs);
    fTreeRB3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB3->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeRB3->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB3->SetCircular(2000000);
    
    
    fTreeBS3->Branch("EventNumber", &EventNumber);
    fTreeBS3->Branch("eventID", &eventID);
    fTreeBS3->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS3->Branch("ReEnergySum", &ReES);
    fTreeBS3->Branch("Pos_eX", &Pos_eX);
    fTreeBS3->Branch("Pos_eY", &Pos_eY);
    fTreeBS3->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS3->Branch("Pos_pX", &Pos_pX);
    fTreeBS3->Branch("Pos_pY", &Pos_pY);
    fTreeBS3->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS3->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS3->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS3->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBS3->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBS3->Branch("PosX_Scat", &PosX_scat);
    fTreeBS3->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS3->Branch("PosY_Scat", &PosY_scat);
    fTreeBS3->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS3->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBS3->Branch("Energy_Scat", &Energy_scat);
    fTreeBS3->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS3->Branch("PosX_Abs", &PosX_abs);
    fTreeBS3->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS3->Branch("PosY_Abs", &PosY_abs);
    fTreeBS3->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS3->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS3->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBS3->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS3->SetCircular(2000000);
    
    
    fTreeS4->Branch("EventNumber", &EventNumber);
    fTreeS4->Branch("eventID", &eventID);
    fTreeS4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS4->Branch("ReEnergySum", &ReES);
    fTreeS4->Branch("Pos_eX", &Pos_eX);
    fTreeS4->Branch("Pos_eY", &Pos_eY);
    fTreeS4->Branch("Pos_eZ", &Pos_eZ);
    fTreeS4->Branch("Pos_pX", &Pos_pX);
    fTreeS4->Branch("Pos_pY", &Pos_pY);
    fTreeS4->Branch("Pos_pZ", &Pos_pZ);
    fTreeS4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS4->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeS4->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeS4->Branch("PosX_Scat", &PosX_scat);
    fTreeS4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS4->Branch("PosY_Scat", &PosY_scat);
    fTreeS4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeS4->Branch("Energy_Scat", &Energy_scat);
    fTreeS4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS4->Branch("PosX_Abs", &PosX_abs);
    fTreeS4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS4->Branch("PosY_Abs", &PosY_abs);
    fTreeS4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS4->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeS4->Branch("AngularDistribution", &AngularDistribution);
    fTreeS4->SetCircular(2000000); 
    
    fTreeB4->Branch("EventNumber", &EventNumber);
    fTreeB4->Branch("eventID", &eventID);
    fTreeB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB4->Branch("ReEnergySum", &ReES);
    fTreeB4->Branch("Pos_eX", &Pos_eX);
    fTreeB4->Branch("Pos_eY", &Pos_eY);
    fTreeB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeB4->Branch("Pos_pX", &Pos_pX);
    fTreeB4->Branch("Pos_pY", &Pos_pY);
    fTreeB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB4->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeB4->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeB4->Branch("PosX_Scat", &PosX_scat);
    fTreeB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB4->Branch("PosY_Scat", &PosY_scat);
    fTreeB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeB4->Branch("Energy_Scat", &Energy_scat);
    fTreeB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB4->Branch("PosX_Abs", &PosX_abs);
    fTreeB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB4->Branch("PosY_Abs", &PosY_abs);
    fTreeB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeB4->SetCircular(2000000);
    
    
    fTreeBB4->Branch("EventNumber", &EventNumber);
    fTreeBB4->Branch("eventID", &eventID);
    fTreeBB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB4->Branch("ReEnergySum", &ReES);
    fTreeBB4->Branch("Pos_eX", &Pos_eX);
    fTreeBB4->Branch("Pos_eY", &Pos_eY);
    fTreeBB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB4->Branch("Pos_pX", &Pos_pX);
    fTreeBB4->Branch("Pos_pY", &Pos_pY);
    fTreeBB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB4->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBB4->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBB4->Branch("PosX_Scat", &PosX_scat);
    fTreeBB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB4->Branch("PosY_Scat", &PosY_scat);
    fTreeBB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBB4->Branch("Energy_Scat", &Energy_scat);
    fTreeBB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB4->Branch("PosX_Abs", &PosX_abs);
    fTreeBB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB4->Branch("PosY_Abs", &PosY_abs);
    fTreeBB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB4->SetCircular(2000000);
    
    fTreedB4->Branch("EventNumber", &EventNumber);
    fTreedB4->Branch("eventID", &eventID);
    fTreedB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreedB4->Branch("ReEnergySum", &ReES);
    fTreedB4->Branch("Pos_eX", &Pos_eX);
    fTreedB4->Branch("Pos_eY", &Pos_eY);
    fTreedB4->Branch("Pos_eZ", &Pos_eZ);
    fTreedB4->Branch("Pos_pX", &Pos_pX);
    fTreedB4->Branch("Pos_pY", &Pos_pY);
    fTreedB4->Branch("Pos_pZ", &Pos_pZ);
    fTreedB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreedB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreedB4->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreedB4->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreedB4->Branch("PosX_Scat", &PosX_scat);
    fTreedB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreedB4->Branch("PosY_Scat", &PosY_scat);
    fTreedB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreedB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreedB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreedB4->Branch("Energy_Scat", &Energy_scat);
    fTreedB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreedB4->Branch("PosX_Abs", &PosX_abs);
    fTreedB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreedB4->Branch("PosY_Abs", &PosY_abs);
    fTreedB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreedB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreedB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreedB4->Branch("AngularDistribution", &AngularDistribution);
    fTreedB4->SetCircular(2000000);
    
    fTreeRB4->Branch("EventNumber", &EventNumber);
    fTreeRB4->Branch("eventID", &eventID);
    fTreeRB4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeRB4->Branch("ReEnergySum", &ReES);
    fTreeRB4->Branch("Pos_eX", &Pos_eX);
    fTreeRB4->Branch("Pos_eY", &Pos_eY);
    fTreeRB4->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB4->Branch("Pos_pX", &Pos_pX);
    fTreeRB4->Branch("Pos_pY", &Pos_pY);
    fTreeRB4->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB4->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeRB4->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeRB4->Branch("PosX_Scat", &PosX_scat);
    fTreeRB4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB4->Branch("PosY_Scat", &PosY_scat);
    fTreeRB4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeRB4->Branch("Energy_Scat", &Energy_scat);
    fTreeRB4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB4->Branch("PosX_Abs", &PosX_abs);
    fTreeRB4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB4->Branch("PosY_Abs", &PosY_abs);
    fTreeRB4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB4->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeRB4->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB4->SetCircular(2000000);
    
    
    fTreeBS4->Branch("EventNumber", &EventNumber);
    fTreeBS4->Branch("eventID", &eventID);
    fTreeBS4->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS4->Branch("ReEnergySum", &ReES);
    fTreeBS4->Branch("Pos_eX", &Pos_eX);
    fTreeBS4->Branch("Pos_eY", &Pos_eY);
    fTreeBS4->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS4->Branch("Pos_pX", &Pos_pX);
    fTreeBS4->Branch("Pos_pY", &Pos_pY);
    fTreeBS4->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS4->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS4->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS4->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBS4->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBS4->Branch("PosX_Scat", &PosX_scat);
    fTreeBS4->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS4->Branch("PosY_Scat", &PosY_scat);
    fTreeBS4->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS4->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBS4->Branch("Energy_Scat", &Energy_scat);
    fTreeBS4->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS4->Branch("PosX_Abs", &PosX_abs);
    fTreeBS4->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS4->Branch("PosY_Abs", &PosY_abs);
    fTreeBS4->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS4->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS4->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBS4->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS4->SetCircular(2000000);
    
    
    fTreeS5->Branch("EventNumber", &EventNumber);
    fTreeS5->Branch("eventID", &eventID);
    fTreeS5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeS5->Branch("ReEnergySum", &ReES);
    fTreeS5->Branch("Pos_eX", &Pos_eX);
    fTreeS5->Branch("Pos_eY", &Pos_eY);
    fTreeS5->Branch("Pos_eZ", &Pos_eZ);
    fTreeS5->Branch("Pos_pX", &Pos_pX);
    fTreeS5->Branch("Pos_pY", &Pos_pY);
    fTreeS5->Branch("Pos_pZ", &Pos_pZ);
    fTreeS5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeS5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeS5->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeS5->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeS5->Branch("PosX_Scat", &PosX_scat);
    fTreeS5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeS5->Branch("PosY_Scat", &PosY_scat);
    fTreeS5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeS5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeS5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeS5->Branch("Energy_Scat", &Energy_scat);
    fTreeS5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeS5->Branch("PosX_Abs", &PosX_abs);
    fTreeS5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeS5->Branch("PosY_Abs", &PosY_abs);
    fTreeS5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeS5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeS5->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeS5->Branch("AngularDistribution", &AngularDistribution);
    fTreeS5->SetCircular(2000000); 
    
    fTreeB5->Branch("EventNumber", &EventNumber);
    fTreeB5->Branch("eventID", &eventID);
    fTreeB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeB5->Branch("ReEnergySum", &ReES);
    fTreeB5->Branch("Pos_eX", &Pos_eX);
    fTreeB5->Branch("Pos_eY", &Pos_eY);
    fTreeB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeB5->Branch("Pos_pX", &Pos_pX);
    fTreeB5->Branch("Pos_pY", &Pos_pY);
    fTreeB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeB5->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeB5->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeB5->Branch("PosX_Scat", &PosX_scat);
    fTreeB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeB5->Branch("PosY_Scat", &PosY_scat);
    fTreeB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeB5->Branch("Energy_Scat", &Energy_scat);
    fTreeB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeB5->Branch("PosX_Abs", &PosX_abs);
    fTreeB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeB5->Branch("PosY_Abs", &PosY_abs);
    fTreeB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeB5->SetCircular(2000000);
    
    
    fTreeBB5->Branch("EventNumber", &EventNumber);
    fTreeBB5->Branch("eventID", &eventID);
    fTreeBB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBB5->Branch("ReEnergySum", &ReES);
    fTreeBB5->Branch("Pos_eX", &Pos_eX);
    fTreeBB5->Branch("Pos_eY", &Pos_eY);
    fTreeBB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeBB5->Branch("Pos_pX", &Pos_pX);
    fTreeBB5->Branch("Pos_pY", &Pos_pY);
    fTreeBB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeBB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBB5->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBB5->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBB5->Branch("PosX_Scat", &PosX_scat);
    fTreeBB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBB5->Branch("PosY_Scat", &PosY_scat);
    fTreeBB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBB5->Branch("Energy_Scat", &Energy_scat);
    fTreeBB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBB5->Branch("PosX_Abs", &PosX_abs);
    fTreeBB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBB5->Branch("PosY_Abs", &PosY_abs);
    fTreeBB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeBB5->SetCircular(2000000);
    
    fTreedB5->Branch("EventNumber", &EventNumber);
    fTreedB5->Branch("eventID", &eventID);
    fTreedB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreedB5->Branch("ReEnergySum", &ReES);
    fTreedB5->Branch("Pos_eX", &Pos_eX);
    fTreedB5->Branch("Pos_eY", &Pos_eY);
    fTreedB5->Branch("Pos_eZ", &Pos_eZ);
    fTreedB5->Branch("Pos_pX", &Pos_pX);
    fTreedB5->Branch("Pos_pY", &Pos_pY);
    fTreedB5->Branch("Pos_pZ", &Pos_pZ);
    fTreedB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreedB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreedB5->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreedB5->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreedB5->Branch("PosX_Scat", &PosX_scat);
    fTreedB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreedB5->Branch("PosY_Scat", &PosY_scat);
    fTreedB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreedB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreedB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreedB5->Branch("Energy_Scat", &Energy_scat);
    fTreedB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreedB5->Branch("PosX_Abs", &PosX_abs);
    fTreedB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreedB5->Branch("PosY_Abs", &PosY_abs);
    fTreedB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreedB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreedB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreedB5->Branch("AngularDistribution", &AngularDistribution);
    fTreedB5->SetCircular(2000000);
    
    fTreeRB5->Branch("EventNumber", &EventNumber);
    fTreeRB5->Branch("eventID", &eventID);
    fTreeRB5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeRB5->Branch("ReEnergySum", &ReES);
    fTreeRB5->Branch("Pos_eX", &Pos_eX);
    fTreeRB5->Branch("Pos_eY", &Pos_eY);
    fTreeRB5->Branch("Pos_eZ", &Pos_eZ);
    fTreeRB5->Branch("Pos_pX", &Pos_pX);
    fTreeRB5->Branch("Pos_pY", &Pos_pY);
    fTreeRB5->Branch("Pos_pZ", &Pos_pZ);
    fTreeRB5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeRB5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeRB5->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeRB5->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeRB5->Branch("PosX_Scat", &PosX_scat);
    fTreeRB5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeRB5->Branch("PosY_Scat", &PosY_scat);
    fTreeRB5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeRB5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeRB5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeRB5->Branch("Energy_Scat", &Energy_scat);
    fTreeRB5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeRB5->Branch("PosX_Abs", &PosX_abs);
    fTreeRB5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeRB5->Branch("PosY_Abs", &PosY_abs);
    fTreeRB5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeRB5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeRB5->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeRB5->Branch("AngularDistribution", &AngularDistribution);
    fTreeRB5->SetCircular(2000000);
    
    
    fTreeBS5->Branch("EventNumber", &EventNumber);
    fTreeBS5->Branch("eventID", &eventID);
    fTreeBS5->Branch("PrimaryEnergy", &PriEnergy);
    fTreeBS5->Branch("ReEnergySum", &ReES);
    fTreeBS5->Branch("Pos_eX", &Pos_eX);
    fTreeBS5->Branch("Pos_eY", &Pos_eY);
    fTreeBS5->Branch("Pos_eZ", &Pos_eZ);
    fTreeBS5->Branch("Pos_pX", &Pos_pX);
    fTreeBS5->Branch("Pos_pY", &Pos_pY);
    fTreeBS5->Branch("Pos_pZ", &Pos_pZ);
    fTreeBS5->Branch("RealEnergy_e", &RealEnergy_e);
    fTreeBS5->Branch("RealEnergy_p", &RealEnergy_p);
    fTreeBS5->Branch("RecoEnergy_Scat", &ReEnergy_scat);
    fTreeBS5->Branch("RecoEnergy_Abs", &ReEnergy_abs);
    fTreeBS5->Branch("PosX_Scat", &PosX_scat);
    fTreeBS5->Branch("PosX_ScatUnc", &PosX_scatunc);
    fTreeBS5->Branch("PosY_Scat", &PosY_scat);
    fTreeBS5->Branch("PosY_ScatUnc", &PosY_scatunc);
    fTreeBS5->Branch("PosZ_Scat", &PosZ_scat);
    fTreeBS5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
    fTreeBS5->Branch("Energy_Scat", &Energy_scat);
    fTreeBS5->Branch("Energy_ScatUnc", &Energy_scatunc);
    fTreeBS5->Branch("PosX_Abs", &PosX_abs);
    fTreeBS5->Branch("PosX_AbsUnc", &PosX_absunc);
    fTreeBS5->Branch("PosY_Abs", &PosY_abs);
    fTreeBS5->Branch("PosY_AbsUnc", &PosY_absunc);
    fTreeBS5->Branch("PosZ_Abs", &PosZ_abs);
    fTreeBS5->Branch("PosZ_AbsUnc", &PosZ_absunc);
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
    fTreeBS5->Branch("AngularDistribution", &AngularDistribution);
    fTreeBS5->SetCircular(2000000);
    

// --- Create the Reader object

   TMVA::Reader *reader21 = new TMVA::Reader( "!Color:!Silent" );

   TMVA::Reader *reader31 = new TMVA::Reader( "!Color:!Silent" );

   TMVA::Reader *reader41 = new TMVA::Reader( "!Color:!Silent" );
  
   TMVA::Reader *reader51 = new TMVA::Reader( "!Color:!Silent" );

   //reader21->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader21->AddVariable( "Energy_Scat", &Energy_scat);
   //reader21->AddVariable( "EnergySum", &EnergySum);
   reader21->AddVariable( "EnDiff", &EnDiff);

   //reader31->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader31->AddVariable( "Energy_Scat", &Energy_scat);
   //reader31->AddVariable( "EnergySum", &EnergySum);
   reader31->AddVariable( "EnDiff", &EnDiff);

   //reader41->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader41->AddVariable( "Energy_Scat", &Energy_scat);
   //reader41->AddVariable( "EnergySum", &EnergySum);
   reader41->AddVariable( "EnDiff", &EnDiff );

   //reader51->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader51->AddVariable( "Energy_Scat", &Energy_scat);
   //reader51->AddVariable( "EnergySum", &EnergySum);
   reader51->AddVariable( "EnDiff", &EnDiff);

   // load the weight files for the readers
   TString method =  "BDT method";
   
   reader21->BookMVA( method, "dataset2res/weights/TMVARegression_BDTG.weights.xml" );


   reader31->BookMVA( method, "dataset3res/weights/TMVARegression_BDTG.weights.xml" );
   

   reader41->BookMVA( method, "dataset4res/weights/TMVARegression_BDTG.weights.xml" );


   reader51->BookMVA( method, "dataset5res/weights/TMVARegression_BDTG.weights.xml" );

    // load the input file
    TFile *input(0);
    TString fname = "./PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5FirsthalfModified2sigma.root";
    
    input = TFile::Open( fname );
    
    TTree* theTree = NULL;
    
    TStopwatch sw;
        sw.Start();
    // loop through signal and all background trees
    for( int treeNumber = 0; treeNumber < 27; ++treeNumber ) {
        
       
       
        if( treeNumber == 0 ){
            
            
            theTree = (TTree*)input->Get("TreeS2");
            std::cout << "--- Select signal sample" << std::endl;
            
             theTree->SetBranchAddress("EventNumber", &EventNumber);
             theTree->SetBranchAddress("eventID", &eventID);
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
           
             theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
            
            std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
            Int_t nEvent = theTree->GetEntries();
            
            for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader21->EvaluateRegression( method ))[0];
              
               
               fTreeS2->Fill();
            }
            
        } else if( treeNumber == 1 ){
            
           
           theTree = (TTree*)input->Get("TreeB2");
           std::cout << "--- Select background sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           //Int_t nEvent = 100;
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader21->EvaluateRegression( method ))[0];
               
               fTreeB2->Fill();
           }
           
           
       }else if( treeNumber == 2 ){
            
           
           theTree = (TTree*)input->Get("TreeBB2");
           std::cout << "--- Select background sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader21->EvaluateRegression( method ))[0];
               
               fTreeBB2->Fill();
           }
           
           
       }/* else if( treeNumber == 3 ){
      
           theTree = (TTree*)input->Get("TreeRB2");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
           //theTree->SetBranchAddress("PERB", &PERB);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader21->EvaluateRegression( method ))[0];
              
               
               fTreeRB2->Fill();
           }
           
           
           
       } else if( treeNumber == 4 ){
      
           theTree = (TTree*)input->Get("TreeBS2");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader21->EvaluateRegression( method ))[0];
               
               
               fTreeBS2->Fill();
           }
           
           
           
       }*/ else if( treeNumber == 5 ){
           
           theTree = (TTree*)input->Get("TreeS3");
           std::cout << "--- Select signal sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader31->EvaluateRegression( method ))[0];
               
               fTreeS3->Fill();
           }
           
           
       } else if( treeNumber == 6 ){
           
           theTree = (TTree*)input->Get("TreeB3");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader31->EvaluateRegression( method ))[0];
               
               
               fTreeB3->Fill();
           }
           
           
       } else if( treeNumber == 7 ){
           
           theTree = (TTree*)input->Get("TreeBB3");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader31->EvaluateRegression( method ))[0];
               
               
               fTreeBB3->Fill();
           }
           
           
       }/* else if( treeNumber == 8 ){
           
           theTree = (TTree*)input->Get("TreedB3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader31->EvaluateRegression( method ))[0];
               //ReEnergy_scat = (reader32->EvaluateRegression( method ))[0];
               //ReEnergy_abs = (reader33->EvaluateRegression( method ))[0];
               
               fTreedB3->Fill();
           }
           
           
       } else if( treeNumber == 9 ){
           
           theTree = (TTree*)input->Get("TreeRB3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
           //theTree->SetBranchAddress("PERB", &PERB);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader31->EvaluateRegression( method ))[0];
               //ReEnergy_scat = (reader32->EvaluateRegression( method ))[0];
               //ReEnergy_abs = (reader33->EvaluateRegression( method ))[0];
               
               fTreeRB3->Fill();
           }
           
           
       } else if( treeNumber == 10 ){
           
           theTree = (TTree*)input->Get("TreeBS3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader31->EvaluateRegression( method ))[0];
               //ReEnergy_scat = (reader32->EvaluateRegression( method ))[0];
               //ReEnergy_abs = (reader33->EvaluateRegression( method ))[0];
               
               fTreeBS3->Fill();
           }
           
           
       }*/ else if( treeNumber == 11 ){
           
           theTree = (TTree*)input->Get("TreeS4");
           std::cout << "--- Select Real signal sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
             
               ReES = (reader41->EvaluateRegression( method ))[0];
               
               
               fTreeS4->Fill();
           }
           
       } else if( treeNumber == 12 ){
           
           theTree = (TTree*)input->Get("TreeB4");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader41->EvaluateRegression( method ))[0];
               
               
               fTreeB4->Fill();
           }
           
       } else if( treeNumber == 13 ){
           
           theTree = (TTree*)input->Get("TreeBB4");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader41->EvaluateRegression( method ))[0];
               
               
               fTreeBB4->Fill();
           }
           
       }/* else if( treeNumber == 14 ){
           
           theTree = (TTree*)input->Get("TreedB4");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader41->EvaluateRegression( method ))[0];
              
               
               fTreedB4->Fill();
           }
           
       }else if( treeNumber == 15 ){
           
           theTree = (TTree*)input->Get("TreeRB4");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
           //theTree->SetBranchAddress("PERB", &PERB);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
         
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader41->EvaluateRegression( method ))[0];
               
               
               fTreeRB4->Fill();
           }
           
       }else if( treeNumber == 16 ){
           
           theTree = (TTree*)input->Get("TreeBS4");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader41->EvaluateRegression( method ))[0];
               
               
               fTreeBS4->Fill();
           }
           
       }*/ else if( treeNumber == 17 ){
           
           theTree = (TTree*)input->Get("TreeS5");
           std::cout << "--- Select Real signal sample" << std::endl;
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
          
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing signal events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader51->EvaluateRegression( method ))[0];
              
               
               fTreeS5->Fill();
           }
           
       } else if( treeNumber == 18 ){
           
           theTree = (TTree*)input->Get("TreeB5");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
              
               ReES = (reader51->EvaluateRegression( method ))[0];
               
               
               fTreeB5->Fill();
           }
           
       } else if( treeNumber == 19 ){
           
           theTree = (TTree*)input->Get("TreeBB5");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader51->EvaluateRegression( method ))[0];
               
               
               fTreeBB5->Fill();
           }
           
       }/* else if( treeNumber == 20 ){
           
           theTree = (TTree*)input->Get("TreedB5");
           std::cout << "--- Select background sample" << std::endl;
         
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               
               ReES = (reader51->EvaluateRegression( method ))[0];
              
               
               fTreedB5->Fill();
           }
           
       } else if( treeNumber == 21 ){
           
           theTree = (TTree*)input->Get("TreeRB5");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
         
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader51->EvaluateRegression( method ))[0];
               
               
               fTreeRB5->Fill();
           }
           
       } else if( treeNumber == 22 ){
           
           theTree = (TTree*)input->Get("TreeBS5");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           //theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader51->EvaluateRegression( method ))[0];
              
               
               fTreeBS5->Fill();
           }
           
       } else if( treeNumber == 23 ){
           
           theTree = (TTree*)input->Get("TreeSB2");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           //theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               
               ReES = (reader21->EvaluateRegression( method ))[0];
               //eventIDSB2++;
               
               fTreeSB2->Fill();
           }
           
       } else if( treeNumber == 24 ){
           
           theTree = (TTree*)input->Get("TreeSB3");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           //theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               ReES = (reader31->EvaluateRegression( method ))[0];
               //eventIDSB3++;
               
               fTreeSB3->Fill();
           }
           
       } else if( treeNumber == 25 ){
           
           theTree = (TTree*)input->Get("TreeSB4");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           //theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               ReES = (reader41->EvaluateRegression( method ))[0];
               //eventIDSB4++;
               
               fTreeSB4->Fill();
           }
           
       } else if( treeNumber == 26 ){
           
           theTree = (TTree*)input->Get("TreeSB5");
           std::cout << "--- Select background sample" << std::endl;
          
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
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
           theTree->SetBranchAddress("EnergySecond_abs", &EnergySe);
           theTree->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "DiffEnergy", &DiffEnergy );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           //theTree->SetBranchAddress( "ECII", &fECII );
           std::cout << "--- Processing: " << theTree->GetEntries() << " background events" << std::endl;
           
       
           Int_t nEvent = theTree->GetEntries();
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
               if (ievt%10000 == 0){
                   
               
                   std::cout << "--- ... Processing background events: " << ievt << std::endl;
               }
    
               theTree->GetEntry(ievt);
               ReES = (reader51->EvaluateRegression( method ))[0];
               //eventIDSB5++;
               
               fTreeSB5->Fill();
           }
           
       }
*/
        
    }
     // get elapsed time
        sw.Stop();
        std::cout << "--- End of event loop: "; sw.Print();
    
    input->Close();
   // write output tree
/*   outputTree->SetDirectory(outputFile);
     outputTree->Write(); */
    outputFile->Write();
    outputFile->Close();
    std::cout << "--- Created root file: \"" << outfileName.Data() << "\" containing event classes trees" << std::endl;
           
            
}

////////////////////////////////////////////////
/////////////////////////////////////////

void RESRegressionFinal(){
   // ----------------------------------------------------------------------------------------
   // Run all
   // ----------------------------------------------------------------------------------------
//    cout << "Start Test TMVAGAexample" << endl
//         << "========================" << endl
//         << endl;
//    TString createDataMacro = gROOT->GetTutorialDir() + "/tmva/createData.C";
//    gROOT->ProcessLine(TString::Format(".L %s",createDataMacro.Data()));
//    gROOT->ProcessLine("create_MultipleBackground(200)");
//    cout << endl;

   cout << "--- Predicted energy regression" << endl;
   Regression();
   cout << endl;
   cout << "========================" << endl;
   
   
   cout << "--- RegressionApplication" << endl;
   RegressionApplied();
   cout << endl;
   
 

}

int main( int argc, char** argv ) {
    
    RESRegressionFinal();
}
