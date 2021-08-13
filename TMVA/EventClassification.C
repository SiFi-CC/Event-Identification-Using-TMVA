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


////////////////////////////////////////////////////////////////////
/// To do Cross Validation training, please use ROOT version 6.18 or higher
// ----------------------------------------------------------------------------------------------
// Simple Training
// ----------------------------------------------------------------------------------------------
//
void Training(){
    
   std::string factoryOptions( "!V:!Silent:Transformations=N;I;D;P;G,D:AnalysisType=Classification" );
   TString fname = "./PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5FirsthalfModified2sigma.root";
   TFile *input(0);
   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   }
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   
//calling the signal and background events trees of each event class
   
   
   TTree *signal2      = (TTree*)input->Get("TreeS2");
   TTree *background2 = (TTree*)input->Get("TreeB2");

   TTree *signal3      = (TTree*)input->Get("TreeS3");
   TTree *background3 = (TTree*)input->Get("TreeB3");

   TTree *signal4      = (TTree*)input->Get("TreeS4");
   TTree *background4 = (TTree*)input->Get("TreeB4");
   
   TTree *signal5      = (TTree*)input->Get("TreeS5");
   TTree *background5 = (TTree*)input->Get("TreeB5");

   
// global event weights per tree 
   Double_t signalWeight      = 1.0;
   
   Double_t backgroundWeight = 1.0;
   
   
// Create a new root output file.
   TString outfileName( "TMVASignalBackground-s2s3s4s5-SimpleTrainAllVarExceptEabs.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
  
   
// ____________////////////////////////////////////////////////////////////
   TMVA::Factory *factory = new TMVA::Factory( "TMVASigBkg", outputFile, factoryOptions );
   
   TMVA::DataLoader *dataloader2=new TMVA::DataLoader("dataset-S28V");
   
   
   dataloader2->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader2->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   
   dataloader2->AddVariable( "PosX_Scat","PosX_Scat", 'F', 0, 450 );
   dataloader2->AddVariable( "PosY_Scat","PosY_Scat", 'F', -60, 60 );
   dataloader2->AddVariable( "PosZ_Scat","PosZ_Scat", 'F', -60, 60 );
   dataloader2->AddVariable( "PosX_Abs","PosX_Abs", 'F', 0, 450 );
   dataloader2->AddVariable( "PosY_Abs","PosY_Abs", 'F', -60, 60 );
   dataloader2->AddVariable( "PosZ_Abs","PosZ_Abs", 'F', -60, 60 );
   //dataloader2->AddVariable( "Energy_Abs","Energy_Abs", 'F', 0, 20 );
   
   
   dataloader2->AddSpectator("eventID");
   
   dataloader2->AddSignalTree    ( signal2,     signalWeight       );
   dataloader2->AddBackgroundTree( background2, backgroundWeight );
////////////////////////////////////////////////////////////////////////////////////////

   TMVA::DataLoader *dataloader3 = new TMVA::DataLoader("dataset-S39V");
   
   
   dataloader3->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader3->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   dataloader3->AddVariable( "AngularDistribution","AngularDistribution", 'F', -2, 8 );
   
   dataloader3->AddVariable( "PosX_Scat","PosX_Scat", 'F', 0, 450 );
   dataloader3->AddVariable( "PosY_Scat","PosY_Scat", 'F', -60, 60 );
   dataloader3->AddVariable( "PosZ_Scat","PosZ_Scat", 'F', -60, 60 );
   dataloader3->AddVariable( "PosX_Abs","PosX_Abs", 'F', 0, 450 );
   dataloader3->AddVariable( "PosY_Abs","PosY_Abs", 'F', -60, 60 );
   dataloader3->AddVariable( "PosZ_Abs","PosZ_Abs", 'F', -60, 60 );
   //dataloader3->AddVariable( "Energy_Abs","Energy_Abs", 'F', 0, 20 );
   
   dataloader3->AddSpectator("eventID");
   
   dataloader3->AddSignalTree    ( signal3,     signalWeight       );
   dataloader3->AddBackgroundTree( background3, backgroundWeight );

///////////////////////////////////////////////////////////////////////////////////////////
   
   TMVA::DataLoader *dataloader4=new TMVA::DataLoader("dataset-S49V");
   
   
   dataloader4->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader4->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   dataloader4->AddVariable( "AngularDistribution","AngularDistribution", 'F', -2, 8 );
   
   dataloader4->AddVariable( "PosX_Scat","PosX_Scat", 'F', 0, 450 );
   dataloader4->AddVariable( "PosY_Scat","PosY_Scat", 'F', -60, 60 );
   dataloader4->AddVariable( "PosZ_Scat","PosZ_Scat", 'F', -60, 60 );
   dataloader4->AddVariable( "PosX_Abs","PosX_Abs", 'F', 0, 450 );
   dataloader4->AddVariable( "PosY_Abs","PosY_Abs", 'F', -60, 60 );
   dataloader4->AddVariable( "PosZ_Abs","PosZ_Abs", 'F', -60, 60 );
   //dataloader4->AddVariable( "Energy_Abs","Energy_Abs", 'F', 0, 20 );
   
   dataloader4->AddSpectator("eventID");
   
   dataloader4->AddSignalTree    ( signal4,     signalWeight       );
   dataloader4->AddBackgroundTree( background4, backgroundWeight );
   
///////////////////////////////////////////////////////////////////////////////////////////////
   
   TMVA::DataLoader *dataloader5=new TMVA::DataLoader("dataset-S59V");
   
   dataloader5->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader5->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   dataloader5->AddVariable( "AngularDistribution","AngularDistribution", 'F', -2, 8 );
   
   dataloader5->AddVariable( "PosX_Scat","PosX_Scat", 'F', 0, 450 );
   dataloader5->AddVariable( "PosY_Scat","PosY_Scat", 'F', -60, 60 );
   dataloader5->AddVariable( "PosZ_Scat","PosZ_Scat", 'F', -60, 60 );
   dataloader5->AddVariable( "PosX_Abs","PosX_Abs", 'F', 0, 450 );
   dataloader5->AddVariable( "PosY_Abs","PosY_Abs", 'F', -60, 60 );
   dataloader5->AddVariable( "PosZ_Abs","PosZ_Abs", 'F', -60, 60 );
   //dataloader5->AddVariable( "Energy_Abs","Energy_Abs", 'F', 0, 20 );
   
   dataloader5->AddSpectator("eventID");
   
   dataloader5->AddSignalTree    ( signal5,     signalWeight       );
   dataloader5->AddBackgroundTree( background5, backgroundWeight );
   
///////////////////////////////////////////////////////////////////////////////////////////
 
   //dataloader->AddSignalTree    ( signal,     signalWeight, TMVA::Types::kTraining);
   //dataloader->AddBackgroundTree( background, backgroundWeight, TMVA::Types::kTraining);
   
   //dataloader->AddSignalTree    ( signal,     signalWeight, TMVA::Types::kTesting);
   //dataloader->AddBackgroundTree( background, backgroundWeight, TMVA::Types::kTesting);
   
   //     factory->SetBackgroundWeightExpression("weight");
   
   
   TCut mycuts = ""; 
   TCut mycutb = ""; 
// tell the factory to use all remaining events in the trees after training for testing:
   
   dataloader2->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   
   
   dataloader3->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   

   dataloader4->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
  
   
   dataloader5->PrepareTrainingAndTestTree( mycuts, mycutb,
                                        "nTrain_Signal=0:nTest_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
   

// Boosted Decision Trees //////////////////////////////////////////
   
   factory->BookMethod( dataloader2, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" ); 
   
   
   factory->BookMethod( dataloader3, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" ); 
   
   
   factory->BookMethod( dataloader4, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" ); 
   
   
   factory->BookMethod( dataloader5, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" ); 
   
// Gradient Boosted Decision Trees (BDTG) ///////////////////////////////////////////////////
/*   
   factory->BookMethod( dataloader2, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex" ); 
   
   
   factory->BookMethod( dataloader3, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex" ); 
   
   
   factory->BookMethod( dataloader4, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex" ); 
   
   
   factory->BookMethod( dataloader5, TMVA::Types::kBDT,  "BDTG",
                           "!H:!V:NTrees=1000:MinNodeSize=5%:BoostType=Grad:Shrinkage=0.1:UseBaggedGrad:BaggedSampleFraction=0.5:nCuts=30:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:SeparationType=GiniIndex" );    
*/

// MLP //////////////////////////////////////////
/*
   factory->BookMethod( dataloader2, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=MSE:VarTransform=N:NCycles=235:HiddenLayers=N,N:LearningRate=0.003:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator=F:ConvergenceTests=1"  ); 
   
   
   factory->BookMethod( dataloader3, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=MSE:VarTransform=N:NCycles=135:HiddenLayers=N,N:LearningRate=0.003:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator=F:ConvergenceTests=1" ); 
   
   
   factory->BookMethod( dataloader4, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=MSE:VarTransform=N:NCycles=115:HiddenLayers=N,N:LearningRate=0.005:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator=F:ConvergenceTests=1"  ); 
   
   
   factory->BookMethod( dataloader5, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=MSE:VarTransform=N:NCycles=85:HiddenLayers=N,N:LearningRate=0.02:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator=F:ConvergenceTests=1"  ); 
*/

// MLPBFGS /////////// Not regulated ///////////////////////////////
/*
   factory->BookMethod( dataloader2, TMVA::Types::kMLP, "MLPBFGS", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N,N:TestRate=10:TrainingMethod=BFGS:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:Tau=0.5" ); 
   
   
   factory->BookMethod( dataloader3, TMVA::Types::kMLP, "MLPBFGS", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N,N:TestRate=10:TrainingMethod=BFGS:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:Tau=0.5" ); 
   
   
   factory->BookMethod( dataloader4, TMVA::Types::kMLP, "MLPBFGS", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N,N:TestRate=10:TrainingMethod=BFGS:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:Tau=0.5" ); 
   
   
   factory->BookMethod( dataloader5, TMVA::Types::kMLP, "MLPBFGS", "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=200:HiddenLayers=N,N:TestRate=10:TrainingMethod=BFGS:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:Tau=0.5" ); 
  */ 
   
// k-NN //////////////////////////////////////////
/*
   factory->BookMethod( dataloader2, TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=80:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  ); 
   
   
   factory->BookMethod( dataloader3, TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=70:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  ); 
   
  
   factory->BookMethod( dataloader4, TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=70:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  ); 
   
   
   factory->BookMethod( dataloader5, TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=30:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  );  
   */
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();
   outputFile->Close();
   
   //if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
   
   delete factory;
   delete dataloader2;
   delete dataloader3;
   delete dataloader4;
   delete dataloader5;
   
} 
//
//---------------------------------------------------------------------------
// Cross-Validation (CV) training
//----------------------------------------------------------------------------
//int TMVACrossValidation(bool useRandomSplitting = false)
void TMVACrossValidation()
{
    
    // root -l -e 'TMVA::TMVAGui("TMVASB2.root")'
    
    
// This loads the library
   TMVA::Tools::Instance();
 
   // --------------------------------------------------------------------------
 
// Load the data into TTrees. 
   TString filename = "./PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5FirsthalfModified2sigma.root";
   TFile * input = TFile::Open( filename );
    
   TTree *signal2      = (TTree*)input->Get("TreeS2");
   TTree *background2 = (TTree*)input->Get("TreeB2");

   TTree *signal3      = (TTree*)input->Get("TreeS3");
   TTree *background3 = (TTree*)input->Get("TreeB3");

   TTree *signal4      = (TTree*)input->Get("TreeS4");
   TTree *background4 = (TTree*)input->Get("TreeB4");
   
   TTree *signal5      = (TTree*)input->Get("TreeS5");
   TTree *background5 = (TTree*)input->Get("TreeB5");

   Double_t signalWeight      = 1.0;
  
   Double_t backgroundWeight = 1.0;
   
// Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName("TMVASB2345CV.root");
    TFile *outputFile = TFile::Open(outfileName, "RECREATE");

// DataLoader definitions; We declare variables in the tree so that TMVA can
// find them.
    
//////////////////////////////////////////////////////////////
   
   TMVA::DataLoader *dataloader2=new TMVA::DataLoader("datasets2");
   
   
   dataloader2->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader2->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   
   
   dataloader2->AddSpectator("eventID");
   
   dataloader2->AddSignalTree    ( signal2,     signalWeight       );
   dataloader2->AddBackgroundTree( background2, backgroundWeight );
////////////////////////////////////////////////////////////////////////////////////////

   TMVA::DataLoader *dataloader3 = new TMVA::DataLoader("datasets3");
   
   
   dataloader3->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader3->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   dataloader3->AddVariable( "AngularDistribution","AngularDistribution", 'F', -2, 8 );
   
   
   dataloader3->AddSpectator("eventID");
   
   dataloader3->AddSignalTree    ( signal3,     signalWeight       );
   dataloader3->AddBackgroundTree( background3, backgroundWeight );

///////////////////////////////////////////////////////////////////////////////////////////
   
   TMVA::DataLoader *dataloader4=new TMVA::DataLoader("datasets4");
   
   dataloader4->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader4->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   dataloader4->AddVariable( "AngularDistribution","AngularDistribution", 'F', -2, 8 );
   
   dataloader4->AddSpectator("eventID");
   
   dataloader4->AddSignalTree    ( signal4,     signalWeight       );
   dataloader4->AddBackgroundTree( background4, backgroundWeight );
   
///////////////////////////////////////////////////////////////////////////////////////////
  
   TMVA::DataLoader *dataloader5=new TMVA::DataLoader("datasets5");
   
   dataloader5->AddVariable( "EnergyCluster_abs","EnergyCluster_abs", 'F', 0, 20 );
   dataloader5->AddVariable( "Energy_Scat","Energy_Scat", 'F', 0, 20 );
   dataloader5->AddVariable( "AngularDistribution","AngularDistribution", 'F', -2, 8 );
   
   
   dataloader5->AddSpectator("eventID");
   
   dataloader5->AddSignalTree    ( signal5,     signalWeight       );
   dataloader5->AddBackgroundTree( background5, backgroundWeight );
   
///////////////////////////////////////////////////////////////////////////////////////////

   
// The CV mechanism of TMVA splits up the training set into several folds.
// The test set is currently left unused. The `nTest_ClassName=1` assigns
// one event to the the test set for each class and puts the rest in the
// training set. A value of 0 is a special value and would split the
// datasets 50 / 50.
   
   dataloader2->PrepareTrainingAndTestTree("", "",
                                          "nTest_Signal=0"
                                          ":nTest_Background=0"
                                          ":SplitMode=Random"
                                          ":NormMode=NumEvents"
                                          ":!V");
   
  
   dataloader3->PrepareTrainingAndTestTree("", "",
                                          "nTest_Signal=0"
                                          ":nTest_Background=0"
                                          ":SplitMode=Random"
                                          ":NormMode=NumEvents"
                                          ":!V");


   dataloader4->PrepareTrainingAndTestTree("", "",
                                          "nTest_Signal=0"
                                          ":nTest_Background=0"
                                          ":SplitMode=Random"
                                          ":NormMode=NumEvents"
                                          ":!V");
   
  
   
   dataloader5->PrepareTrainingAndTestTree("", "",
                                          "nTest_Signal=0"
                                          ":nTest_Background=0"
                                          ":SplitMode=Random"
                                          ":NormMode=NumEvents"
                                          ":!V");
  
// ------------------------------------------------------------------------
// This sets up a CrossValidation class (which wraps a TMVA::Factory
// internally) for 10-folds cross validation.
   
// The split type can be "Random", "RandomStratified" or "Deterministic".
// For the last option, check the comment below. Random splitting randomises
// the order of events and distributes events as evenly as possible.
// RandomStratified applies the same logic but distributes events within a
// class as evenly as possible over the folds.
   
   
   bool useRandomSplitting = false;
   
   UInt_t numFolds = 10;
   TString analysisType = "Classification";
   
/////// For ROOT version 6.18 or higher, one can use the option below.
   
   TString splitType = (useRandomSplitting) ? "Random" : "Deterministic";
   
/////////////   For lower versions of ROOT 
   
   //TString splitType = "Deterministic";  
   
// One can also use a custom splitting function for producing the folds.
// The example uses a dataset spectator `eventID`.
   
// The idea here is that eventID should be an event number that is integral,
// random and independent of the data, generated only once. This last
// property ensures that if a calibration is changed the same event will
// still be assigned the same fold.
   
// If you want to run TMVACrossValidationApplication, make sure you have 
// run this training with Deterministic splitting type, i.e.
// with the option useRandomSPlitting = false
   
   
   //TString splitExpr = (useRandomSplitting) ? "" : "int(fabs([eventID]))%int([NumFolds])" ;
   
   TString splitExpr = "int(fabs([eventID]))%int([NumFolds])" ;
 
   TString cvOptions = Form("!V"
                            ":!Silent"
                            ":!ModelPersistence"
                            ":!FoldFileOutput"
                            ":OutputEnsembling=Avg"
                            ":AnalysisType=%s"
                            ":SplitType=%s"
                            ":NumFolds=%i"
                            ":SplitExpr=%s",
                            analysisType.Data(), splitType.Data(), numFolds,
                            splitExpr.Data());
   
/////////////   For lower versions of ROOT 
 /*
   TString cvOptions = Form("!V"
                            ":!Silent"
                            ":!ModelPersistence"
                            ":!FoldFileOutput"
                            ":Transformations=I;N;D;P;G,D"
                            ":OutputEnsembling=None"
                            ":AnalysisType=%s"
                            ":NumFolds=%i"
                            ":SplitExpr=%s",
                            analysisType.Data(), numFolds,
                            splitExpr.Data());
 */
/////////////////////////////////////////////////////////////

   TMVA::CrossValidation cv2{"TMVACrossValidation", dataloader2, outputFile, cvOptions};
   
   TMVA::CrossValidation cv3{"TMVACrossValidation", dataloader3, outputFile, cvOptions};
   
   TMVA::CrossValidation cv4{"TMVACrossValidation", dataloader4, outputFile, cvOptions};
   
   TMVA::CrossValidation cv5{"TMVACrossValidation", dataloader5, outputFile, cvOptions};
  
   // --------------------------------------------------------------------------
   
    //
   // Books a method to use for evaluation
   //
   
   cv2.BookMethod(TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=2000:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=200:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" );
   
 
   cv3.BookMethod(TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=2000:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=200:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" );
   
     
   cv4.BookMethod(TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=2000:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=200:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" );
   
   
   cv5.BookMethod(TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=2000:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:CreateMVAPdfs:NbinsMVAPdf=200:NsmoothMVAPdf=10:SeparationType=GiniIndex:nCuts=-1" );
   
// MLP ////////////////////////////////////////// 
/*
  cv2.BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=sigmoid:VarTransform=N:NCycles=235:HiddenLayers=N,N:LearningRate=0.003:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:ConvergenceTests=1"  ); 
   
   
  cv3.BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=sigmoid:VarTransform=N:NCycles=135:HiddenLayers=N,N:LearningRate=0.003:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:ConvergenceTests=1" ); 
   
   
  cv4.BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=sigmoid:VarTransform=N:NCycles=115:HiddenLayers=N,N:LearningRate=0.005:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:ConvergenceTests=1"  ); 
   
   
  cv5.BookMethod( TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=sigmoid:EstimatorType=sigmoid:VarTransform=N:NCycles=85:HiddenLayers=N,N:LearningRate=0.02:TestRate=10:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:UseRegulator:ConvergenceTests=1"  ); 
   
   */
// k-NN //////////////////////////////////////////
/*
   cv2.BookMethod( TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=80:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  ); 
   
   
   cv3.BookMethod( TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=70:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  ); 
   
  
   cv4.BookMethod( TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=70:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  ); 
   
   
   cv5.BookMethod( TMVA::Types::kKNN, "KNN",
                           "!H:nkNN=30:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=T:UseWeight=T:CreateMVAPdfs:NbinsMVAPdf=60:NsmoothMVAPdf=10:Trim"  );    

    */
 
// Train, test and evaluate the booked methods.
// Evaluates the booked methods once for each fold and aggregates the result
// in the specified output file.
   
   
   cv2.Evaluate();
  
   cv3.Evaluate();
   
   cv4.Evaluate();
  
   cv5.Evaluate();
 
 
   // --------------------------------------------------------------------------
   
   //
   // Save the output
   //
   outputFile->Close();
 
   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVACrossValidation is done!" << std::endl;
 
   // --------------------------------------------------------------------------
 
   //
   // Launch the GUI for the root macros
   //
   
   if (!gROOT->IsBatch()) {
      // Draw cv-specific graphs
       
      cout<< "Avg ROC for BDTS2 : " << cv2.GetResults()[0].GetROCAverage() << endl; 
      cout<< "Std ROC for BDTS2 : " << cv2.GetResults()[0].GetROCStandardDeviation() << endl; 
//       cout<< "Avg ROC for MLPS2 : " << cv2.GetResults()[1].GetROCAverage() << endl; 
//       cout<< "Std ROC for MLPS2 : " << cv2.GetResults()[1].GetROCStandardDeviation() << endl;
//       cout<< "Avg ROC for kNNS2 : " << cv2.GetResults()[2].GetROCAverage() << endl; 
//       cout<< "Std ROC for kNNS2 : " << cv2.GetResults()[2].GetROCStandardDeviation() << endl; 
//       cv2.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDT-S2");

      cout<< "Avg ROC for BDTS3 : " << cv3.GetResults()[0].GetROCAverage() << endl;
      cout<< "Std ROC for BDTS3 : " << cv3.GetResults()[0].GetROCStandardDeviation() << endl;
//       cout<< "Avg ROC for MLPS3 : " << cv3.GetResults()[1].GetROCAverage() << endl; 
//       cout<< "Std ROC for MLPS3 : " << cv3.GetResults()[1].GetROCStandardDeviation() << endl;
//       cout<< "Avg ROC for kNNS3 : " << cv3.GetResults()[2].GetROCAverage() << endl; 
//       cout<< "Std ROC for kNNS3 : " << cv3.GetResults()[2].GetROCStandardDeviation() << endl; 
//       cv3.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDT-S3");
      
      cout<< "Avg ROC for BDTS4 : " << cv4.GetResults()[0].GetROCAverage() << endl;
      cout<< "Std ROC for BDTS4 : " << cv4.GetResults()[0].GetROCStandardDeviation() << endl;
//       cout<< "Avg ROC for MLPS4 : " << cv4.GetResults()[1].GetROCAverage() << endl; 
//       cout<< "Std ROC for MLPS4 : " << cv4.GetResults()[1].GetROCStandardDeviation() << endl;
//       cout<< "Avg ROC for kNNS4 : " << cv4.GetResults()[2].GetROCAverage() << endl; 
//       cout<< "Std ROC for kNNS4 : " << cv4.GetResults()[2].GetROCStandardDeviation() << endl; 
//       cv4.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDT-S4");
       
      cout<< "Avg ROC for BDTS5 : " << cv5.GetResults()[0].GetROCAverage() << endl;
      cout<< "Std ROC for BDTS5 : " << cv5.GetResults()[0].GetROCStandardDeviation() << endl; 
//       cout<< "Avg ROC for MLPS5 : " << cv5.GetResults()[1].GetROCAverage() << endl; 
//       cout<< "Std ROC for MLPS5 : " << cv5.GetResults()[1].GetROCStandardDeviation() << endl;
//       cout<< "Avg ROC for kNNS5 : " << cv5.GetResults()[2].GetROCAverage() << endl; 
//       cout<< "Std ROC for kNNS5 : " << cv5.GetResults()[2].GetROCStandardDeviation() << endl; 
//       cv5.GetResults()[0].DrawAvgROCCurve(kTRUE, "Avg ROC for BDT-S5");
      
 
      // You can also use the classical gui
      //TMVA::TMVAGui(outfileName);
   }
   
}
//
//-----------------------------------------------------------------------------------------------
// Application
// ----------------------------------------------------------------------------------------------
//
/// create a summary tree with all signal and background events and for each event class. 
void ApplicationCreateCombinedTree(){
   // Create a new root output file.
   //TString outfileName( "tmva__applied-s2s3s4s5CV.root" ); // CV training
   TString outfileName( "tmva__applied-s2s3s4s5ST-9var.root" ); // Simple Training
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
   TTree* outputTree2 = new TTree("SigBkg2","SigBkg tree");
   TTree* outputTree3 = new TTree("SigBkg3","SigBkg tree");
   TTree* outputTree4 = new TTree("SigBkg4","SigBkg tree");
   TTree* outputTree5 = new TTree("SigBkg5","SigBkg tree");
   
   Float_t Pos_eX, Pos_eY, Pos_eZ, Pos_pX, Pos_pY, Pos_pZ, RealEnergy_e, RealEnergy_p; 
   
   Float_t PosX_scat, PosY_scat, PosZ_scat,Energy_scat, ReEnergy_scat;
   
   Float_t PosX_abs, PosY_abs, PosZ_abs, Energy_abs, ReEnergy_abs;
   
   Float_t PosX_sec, PosY_sec, PosZ_sec;
   Float_t EnergyCluster_abs, AngularDistribution;
   
   Float_t PosX_scatunc, PosY_scatunc, PosZ_scatunc,Energy_scatunc;
   Float_t PosX_absunc, PosY_absunc, PosZ_absunc, Energy_absunc, EnDiff;
   Float_t PosX_secunc, PosY_secunc, PosZ_secunc;
   
   Float_t EnergyCluster_absunc,EnergySe, EnergySeunc;
   
   Int_t Multiplicity;
   
   Float_t DiffEnergy, DiffPosition, PriEnergy, EnergySum, EnergySumunc;
   
   Int_t count2 = 0, count3 = 0, count4 = 0, count5 = 0;
   Int_t EventNumber;
   Float_t PosX_1, PosY_1, PosZ_1, Energy_1;
   Float_t PosX_2, PosY_2, PosZ_2, Energy_2;
   Float_t DEnergy, REnergy, DPosition, PEnergy, EnergyS;
  
   Int_t   classID = 0;
   Float_t weight = 1.f;
   Float_t classifier21, classifier22, classifier23, classifier31, classifier32, classifier33, classifier34, classifier41, classifier42, classifier43, classifier44, classifier51, classifier52, classifier53, classifier54;
   Float_t regval1, regval2, regval3;
   
   Float_t ReES,ReESunc, sigma;
   
   Int_t eventID;
   
   Float_t prob21, prob22, prob23, prob31, prob32, prob33, prob34, prob41, prob42, prob43, prob44, prob51, prob52, prob53, prob54;  
   Float_t pSig21, pSig22, pSig23, pSig31, pSig32, pSig33, pSig34, pSig41, pSig42, pSig43, pSig44, pSig51, pSig52, pSig53, pSig54;
   Float_t rar21, rar22, rar23, rar31, rar32, rar33, rar34, rar41, rar42, rar43, rar44, rar51, rar52, rar53, rar54;
   
   TLorentzVector *fECII = new TLorentzVector();
/// Real info. are being saved to find a potentially useful discriminator
/// during the study. Of course, they won't be used in the analysis phase at all.
///Only clusters' positions and energies from the detector response are applicable in whole study. 
   outputTree2->Branch("EventNumber", &EventNumber);
   outputTree2->Branch("PrimaryEnergy", &PriEnergy);
   outputTree2->Branch("eventID", &eventID);
   outputTree2->Branch("ECII", &fECII);
   outputTree2->Branch("Pos_eX", &Pos_eX);
   outputTree2->Branch("Pos_eY", &Pos_eY);
   outputTree2->Branch("Pos_eZ", &Pos_eZ);
   outputTree2->Branch("Pos_pX", &Pos_pX);
   outputTree2->Branch("Pos_pY", &Pos_pY);
   outputTree2->Branch("Pos_pZ", &Pos_pZ);
   outputTree2->Branch("RealEnergy_e", &RealEnergy_e);
   outputTree2->Branch("RealEnergy_p", &RealEnergy_p);
   outputTree2->Branch("PosX_Scat", &PosX_scat);
   outputTree2->Branch("PosX_ScatUnc", &PosX_scatunc);
   outputTree2->Branch("PosY_Scat", &PosY_scat);
   outputTree2->Branch("PosY_ScatUnc", &PosY_scatunc);
   outputTree2->Branch("PosZ_Scat", &PosZ_scat);
   outputTree2->Branch("PosZ_ScatUnc", &PosZ_scatunc);
   outputTree2->Branch("Energy_Scat", &Energy_scat);
   outputTree2->Branch("Energy_ScatUnc", &Energy_scatunc);
   outputTree2->Branch("PosX_Abs", &PosX_abs);
   outputTree2->Branch("PosX_AbsUnc", &PosX_absunc);
   outputTree2->Branch("PosY_Abs", &PosY_abs);
   outputTree2->Branch("PosY_AbsUnc", &PosY_absunc);
   outputTree2->Branch("PosZ_Abs", &PosZ_abs);
   outputTree2->Branch("PosZ_AbsUnc", &PosZ_absunc);
   outputTree2->Branch("EnergyCluster_abs", &EnergyCluster_abs);
   outputTree2->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc);
   outputTree2->Branch("Energy_Abs", &Energy_abs);
   outputTree2->Branch("Energy_AbsUnc", &Energy_absunc);
   outputTree2->Branch("EnDiff", &EnDiff);
   outputTree2->Branch("EnergySecond_abs", &EnergySe);
   outputTree2->Branch("EnergySecond_absUnc", &EnergySeunc);
   outputTree2->Branch("EnergySum", &EnergySum);
   outputTree2->Branch("EnergySumUnc", &EnergySumunc);
   outputTree2->Branch("PosX_Sec", &PosX_sec);
   outputTree2->Branch("PosX_SecUnc", &PosX_secunc);
   outputTree2->Branch("PosY_Sec", &PosY_sec);
   outputTree2->Branch("PosY_SecUnc", &PosY_secunc);
   outputTree2->Branch("PosZ_Sec", &PosZ_sec);
   outputTree2->Branch("PosZ_SecUnc", &PosZ_secunc);
   outputTree2->Branch("DiffPosition", &DiffPosition);
   outputTree2->Branch("AngularDistribution", &AngularDistribution);
   outputTree2->Branch("classID", &classID, "classID/I");
   outputTree2->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree2->Branch("cls21", &classifier21, "cls21/F");
   outputTree2->Branch("rarity21", &rar21, "rarity21/F");
   
   
   outputTree3->Branch("EventNumber", &EventNumber);
   outputTree3->Branch("PrimaryEnergy", &PriEnergy);
   outputTree3->Branch("eventID", &eventID);
   outputTree3->Branch("ECII", &fECII);
   outputTree3->Branch("Pos_eX", &Pos_eX);
   outputTree3->Branch("Pos_eY", &Pos_eY);
   outputTree3->Branch("Pos_eZ", &Pos_eZ);
   outputTree3->Branch("Pos_pX", &Pos_pX);
   outputTree3->Branch("Pos_pY", &Pos_pY);
   outputTree3->Branch("Pos_pZ", &Pos_pZ);
   outputTree3->Branch("RealEnergy_e", &RealEnergy_e);
   outputTree3->Branch("RealEnergy_p", &RealEnergy_p);
   outputTree3->Branch("PosX_Scat", &PosX_scat);
   outputTree3->Branch("PosX_ScatUnc", &PosX_scatunc);
   outputTree3->Branch("PosY_Scat", &PosY_scat);
   outputTree3->Branch("PosY_ScatUnc", &PosY_scatunc);
   outputTree3->Branch("PosZ_Scat", &PosZ_scat);
   outputTree3->Branch("PosZ_ScatUnc", &PosZ_scatunc);
   outputTree3->Branch("Energy_Scat", &Energy_scat);
   outputTree3->Branch("Energy_ScatUnc", &Energy_scatunc);
   outputTree3->Branch("PosX_Abs", &PosX_abs);
   outputTree3->Branch("PosX_AbsUnc", &PosX_absunc);
   outputTree3->Branch("PosY_Abs", &PosY_abs);
   outputTree3->Branch("PosY_AbsUnc", &PosY_absunc);
   outputTree3->Branch("PosZ_Abs", &PosZ_abs);
   outputTree3->Branch("PosZ_AbsUnc", &PosZ_absunc);
   outputTree3->Branch("EnergyCluster_abs", &EnergyCluster_abs);
   outputTree3->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc );
   outputTree3->Branch("Energy_Abs", &Energy_abs);
   outputTree3->Branch("Energy_AbsUnc", &Energy_absunc);
   outputTree3->Branch("EnDiff", &EnDiff);
   outputTree3->Branch("EnergySecond_abs", &EnergySe);
   outputTree3->Branch("EnergySecond_absUnc", &EnergySeunc);
   outputTree3->Branch("EnergySum", &EnergySum);
   outputTree3->Branch("EnergySumUnc", &EnergySumunc);
   outputTree3->Branch("PosX_Sec", &PosX_sec);
   outputTree3->Branch("PosX_SecUnc", &PosX_secunc);
   outputTree3->Branch("PosY_Sec", &PosY_sec);
   outputTree3->Branch("PosY_SecUnc", &PosY_secunc);
   outputTree3->Branch("PosZ_Sec", &PosZ_sec);
   outputTree3->Branch("PosZ_SecUnc", &PosZ_secunc);
   outputTree3->Branch("DiffPosition", &DiffPosition);
   outputTree3->Branch("AngularDistribution", &AngularDistribution);
   outputTree3->Branch("classID", &classID, "classID/I");
   outputTree3->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree3->Branch("cls31", &classifier31, "cls31/F");
   outputTree3->Branch("rarity31", &rar31, "rarity31/F");
   
   
   outputTree4->Branch("EventNumber", &EventNumber);
   outputTree4->Branch("PrimaryEnergy", &PriEnergy);
   outputTree4->Branch("eventID", &eventID);
   outputTree4->Branch("ECII", &fECII);
   outputTree4->Branch("Pos_eX", &Pos_eX);
   outputTree4->Branch("Pos_eY", &Pos_eY);
   outputTree4->Branch("Pos_eZ", &Pos_eZ);
   outputTree4->Branch("Pos_pX", &Pos_pX);
   outputTree4->Branch("Pos_pY", &Pos_pY);
   outputTree4->Branch("Pos_pZ", &Pos_pZ);
   outputTree4->Branch("RealEnergy_e", &RealEnergy_e);
   outputTree4->Branch("RealEnergy_p", &RealEnergy_p);
   outputTree4->Branch("PosX_Scat", &PosX_scat);
   outputTree4->Branch("PosX_ScatUnc", &PosX_scatunc);
   outputTree4->Branch("PosY_Scat", &PosY_scat);
   outputTree4->Branch("PosY_ScatUnc", &PosY_scatunc);
   outputTree4->Branch("PosZ_Scat", &PosZ_scat);
   outputTree4->Branch("PosZ_ScatUnc", &PosZ_scatunc);
   outputTree4->Branch("Energy_Scat", &Energy_scat);
   outputTree4->Branch("Energy_ScatUnc", &Energy_scatunc);
   outputTree4->Branch("PosX_Abs", &PosX_abs);
   outputTree4->Branch("PosX_AbsUnc", &PosX_absunc);
   outputTree4->Branch("PosY_Abs", &PosY_abs);
   outputTree4->Branch("PosY_AbsUnc", &PosY_absunc);
   outputTree4->Branch("PosZ_Abs", &PosZ_abs);
   outputTree4->Branch("PosZ_AbsUnc", &PosZ_absunc);
   outputTree4->Branch("EnergyCluster_abs", &EnergyCluster_abs);
   outputTree4->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc );
   outputTree4->Branch("Energy_Abs", &Energy_abs);
   outputTree4->Branch("Energy_AbsUnc", &Energy_absunc);
   outputTree4->Branch("EnDiff", &EnDiff);
   outputTree4->Branch("EnergySecond_abs", &EnergySe);
   outputTree4->Branch("EnergySecond_absUnc", &EnergySeunc);
   outputTree4->Branch("EnergySum", &EnergySum);
   outputTree4->Branch("EnergySumUnc", &EnergySumunc);
   outputTree4->Branch("PosX_Sec", &PosX_sec);
   outputTree4->Branch("PosX_SecUnc", &PosX_secunc);
   outputTree4->Branch("PosY_Sec", &PosY_sec);
   outputTree4->Branch("PosY_SecUnc", &PosY_secunc);
   outputTree4->Branch("PosZ_Sec", &PosZ_sec);
   outputTree4->Branch("PosZ_SecUnc", &PosZ_secunc);
   outputTree4->Branch("DiffPosition", &DiffPosition);
   outputTree4->Branch("AngularDistribution", &AngularDistribution);
   outputTree4->Branch("classID", &classID, "classID/I");
   outputTree4->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree4->Branch("cls41", &classifier41, "cls41/F");
   outputTree4->Branch("rarity41", &rar41, "rarity41/F");
   
   
   outputTree5->Branch("EventNumber", &EventNumber);
   outputTree5->Branch("PrimaryEnergy", &PriEnergy);
   outputTree5->Branch("eventID", &eventID);
   outputTree5->Branch("ECII", &fECII);
   outputTree5->Branch("Pos_eX", &Pos_eX);
   outputTree5->Branch("Pos_eY", &Pos_eY);
   outputTree5->Branch("Pos_eZ", &Pos_eZ);
   outputTree5->Branch("Pos_pX", &Pos_pX);
   outputTree5->Branch("Pos_pY", &Pos_pY);
   outputTree5->Branch("Pos_pZ", &Pos_pZ);
   outputTree5->Branch("RealEnergy_e", &RealEnergy_e);
   outputTree5->Branch("RealEnergy_p", &RealEnergy_p);
   outputTree5->Branch("PosX_Scat", &PosX_scat);
   outputTree5->Branch("PosX_ScatUnc", &PosX_scatunc);
   outputTree5->Branch("PosY_Scat", &PosY_scat);
   outputTree5->Branch("PosY_ScatUnc", &PosY_scatunc);
   outputTree5->Branch("PosZ_Scat", &PosZ_scat);
   outputTree5->Branch("PosZ_ScatUnc", &PosZ_scatunc);
   outputTree5->Branch("Energy_Scat", &Energy_scat);
   outputTree5->Branch("Energy_ScatUnc", &Energy_scatunc);
   outputTree5->Branch("PosX_Abs", &PosX_abs);
   outputTree5->Branch("PosX_AbsUnc", &PosX_absunc);
   outputTree5->Branch("PosY_Abs", &PosY_abs);
   outputTree5->Branch("PosY_AbsUnc", &PosY_absunc);
   outputTree5->Branch("PosZ_Abs", &PosZ_abs);
   outputTree5->Branch("PosZ_AbsUnc", &PosZ_absunc);
   outputTree5->Branch("EnergyCluster_abs", &EnergyCluster_abs);
   outputTree5->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc );
   outputTree5->Branch("Energy_Abs", &Energy_abs);
   outputTree5->Branch("Energy_AbsUnc", &Energy_absunc);
   outputTree5->Branch("EnDiff", &EnDiff);
   outputTree5->Branch("EnergySecond_abs", &EnergySe);
   outputTree5->Branch("EnergySecond_absUnc", &EnergySeunc);
   outputTree5->Branch("EnergySum", &EnergySum);
   outputTree5->Branch("EnergySumUnc", &EnergySumunc);
   outputTree5->Branch("PosX_Sec", &PosX_sec);
   outputTree5->Branch("PosX_SecUnc", &PosX_secunc);
   outputTree5->Branch("PosY_Sec", &PosY_sec);
   outputTree5->Branch("PosY_SecUnc", &PosY_secunc);
   outputTree5->Branch("PosZ_Sec", &PosZ_sec);
   outputTree5->Branch("PosZ_SecUnc", &PosZ_secunc);
   outputTree5->Branch("DiffPosition", &DiffPosition);
   outputTree5->Branch("AngularDistribution", &AngularDistribution);
   outputTree5->Branch("classID", &classID, "classID/I");
   outputTree5->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree5->Branch("cls51", &classifier51, "cls51/F");
   outputTree5->Branch("rarity51", &rar51, "rarity51/F");
   
  
   // create four readers for the four different signal/background classifications,
   
   TMVA::Reader *reader21 = new TMVA::Reader( "!Color:!Silent" );
   
   TMVA::Reader *reader31 = new TMVA::Reader( "!Color:!Silent" );
   
   TMVA::Reader *reader41 = new TMVA::Reader( "!Color:!Silent" );
  
   TMVA::Reader *reader51 = new TMVA::Reader( "!Color:!Silent" );
 

   
   reader21->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader21->AddVariable( "Energy_Scat", &Energy_scat);
   
   reader21->AddVariable( "PosX_Scat", &PosX_scat );
   reader21->AddVariable( "PosY_Scat", &PosY_scat );
   reader21->AddVariable( "PosZ_Scat", &PosZ_scat );
   reader21->AddVariable( "PosX_Abs", &PosX_abs );
   reader21->AddVariable( "PosY_Abs", &PosY_abs );
   reader21->AddVariable( "PosZ_Abs", &PosZ_abs );
   
   reader21->AddSpectator("eventID", &eventID);
   

   
   reader31->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader31->AddVariable( "Energy_Scat", &Energy_scat);
   reader31->AddVariable( "AngularDistribution", &AngularDistribution);
   
   reader31->AddVariable( "PosX_Scat", &PosX_scat );
   reader31->AddVariable( "PosY_Scat", &PosY_scat );
   reader31->AddVariable( "PosZ_Scat", &PosZ_scat );
   reader31->AddVariable( "PosX_Abs", &PosX_abs );
   reader31->AddVariable( "PosY_Abs", &PosY_abs );
   reader31->AddVariable( "PosZ_Abs", &PosZ_abs );
   
   reader31->AddSpectator("eventID", &eventID);
   

   
   reader41->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader41->AddVariable( "Energy_Scat", &Energy_scat);
   reader41->AddVariable( "AngularDistribution", &AngularDistribution);
   
   reader41->AddVariable( "PosX_Scat", &PosX_scat );
   reader41->AddVariable( "PosY_Scat", &PosY_scat );
   reader41->AddVariable( "PosZ_Scat", &PosZ_scat );
   reader41->AddVariable( "PosX_Abs", &PosX_abs );
   reader41->AddVariable( "PosY_Abs", &PosY_abs );
   reader41->AddVariable( "PosZ_Abs", &PosZ_abs );
   
   reader41->AddSpectator("eventID", &eventID);
   
   
   reader51->AddVariable( "EnergyCluster_abs", &EnergyCluster_abs );
   reader51->AddVariable( "Energy_Scat", &Energy_scat);
   reader51->AddVariable( "AngularDistribution", &AngularDistribution);
  
   reader51->AddVariable( "PosX_Scat", &PosX_scat );
   reader51->AddVariable( "PosY_Scat", &PosY_scat );
   reader51->AddVariable( "PosZ_Scat", &PosZ_scat );
   reader51->AddVariable( "PosX_Abs", &PosX_abs );
   reader51->AddVariable( "PosY_Abs", &PosY_abs );
   reader51->AddVariable( "PosZ_Abs", &PosZ_abs );
   
   reader51->AddSpectator("eventID", &eventID);
   
  
   // load the weight files for the readers
   TString method =  "BDT method";
   
   /////////////////////// Simple training-testing ////////////////////
   
   reader21->BookMVA( method, "dataset-S28V/weights/TMVASigBkg_BDT.weights.xml" );
   
   reader31->BookMVA( method, "dataset-S39V/weights/TMVASigBkg_BDT.weights.xml" );
   
   reader41->BookMVA( method, "dataset-S49V/weights/TMVASigBkg_BDT.weights.xml" );
   
   reader51->BookMVA( method, "dataset-S59V/weights/TMVASigBkg_BDT.weights.xml" );
   
  
////////////////////////////////////////////// Cross-Validation /////////////////////////////////////
  /* 
   reader21->BookMVA( method, "datasets2/weights/TMVACrossValidation_BDT.weights.xml" );
   
   reader31->BookMVA( method, "datasets3/weights/TMVACrossValidation_BDT.weights.xml" );
  
   reader41->BookMVA( method, "datasets4/weights/TMVACrossValidation_BDT.weights.xml" );
   
   reader51->BookMVA( method, "datasets5/weights/TMVACrossValidation_BDT.weights.xml" );
  
 */  
///////////// For CV in TMVA, The GetProbability variable is not available!!! Use GetRarity instead////////////////////////
   
// load the input file
   TFile *input(0);
   TString fname = "./PMMA180MeV0mmBPType3_EnoughStatistics_EI-s2s3s4s5SecondhalfModified2sigma.root";

   input = TFile::Open( fname );
   
   TTree* theTree2 = NULL;
   TTree* theTree3 = NULL;
   TTree* theTree4 = NULL;
   TTree* theTree5 = NULL;
   
   TStopwatch sw2;
   sw2.Start();
   
   // loop through signal and all background trees
   for( int treeNumber = 0; treeNumber < 2; ++treeNumber ) {
       

       if( treeNumber == 0 ){
      
           theTree2 = (TTree*)input->Get("TreeS2");
           std::cout << "--- Select signal sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
          
           weight = 1;
           classID = 0;
       } else if( treeNumber == 1 ){
           
           theTree2 = (TTree*)input->Get("TreeB2");
           std::cout << "--- Select background sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 1;
       }
       
       theTree2->SetBranchAddress("EventNumber", &EventNumber);
       theTree2->SetBranchAddress("PrimaryEnergy", &PriEnergy);
       theTree2->SetBranchAddress("eventID", &eventID);
       theTree2->SetBranchAddress("Pos_eX", &Pos_eX);
       theTree2->SetBranchAddress("Pos_eY", &Pos_eY);
       theTree2->SetBranchAddress("Pos_eZ", &Pos_eZ);
       theTree2->SetBranchAddress("Pos_pX", &Pos_pX);
       theTree2->SetBranchAddress("Pos_pY", &Pos_pY);
       theTree2->SetBranchAddress("Pos_pZ", &Pos_pZ);
       theTree2->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
       theTree2->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
       theTree2->SetBranchAddress( "PosX_Scat", &PosX_scat );
       theTree2->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
       theTree2->SetBranchAddress( "PosY_Scat", &PosY_scat );
       theTree2->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
       theTree2->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
       theTree2->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
       theTree2->SetBranchAddress( "Energy_Scat", &Energy_scat );
       theTree2->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
       theTree2->SetBranchAddress( "PosX_Abs", &PosX_abs );
       theTree2->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
       theTree2->SetBranchAddress( "PosY_Abs", &PosY_abs );
       theTree2->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
       theTree2->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
       theTree2->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
       theTree2->SetBranchAddress( "PosX_Sec", &PosX_sec );
       theTree2->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
       theTree2->SetBranchAddress( "PosY_Sec", &PosY_sec );
       theTree2->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
       theTree2->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
       theTree2->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
       theTree2->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
       theTree2->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
       theTree2->SetBranchAddress( "Energy_Abs", &Energy_abs );
       theTree2->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
       theTree2->SetBranchAddress("EnDiff", &EnDiff);
       theTree2->SetBranchAddress("EnergySecond_abs", &EnergySe);
       theTree2->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
       theTree2->SetBranchAddress( "EnergySum", &EnergySum );
       theTree2->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
       theTree2->SetBranchAddress( "Multiplicity", &Multiplicity );
       theTree2->SetBranchAddress( "DiffPosition", &DiffPosition );
       theTree2->SetBranchAddress( "ECII", &fECII );
       theTree2->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           
       std::cout << "--- Processing: " << theTree2->GetEntries() << " events" << std::endl;
           
      
           
       Int_t nEvent2 = theTree2->GetEntries();
          
       for (Long64_t ievt=0; ievt<nEvent2; ievt++) {
           
           
           if (ievt%10000 == 0){
 
               std::cout << "--- ... Processing events: " << ievt << std::endl;
           }
    
           theTree2->GetEntry(ievt);
               
           // get the classifier's response for each of the signal/background classifications
           classifier21 = reader21->EvaluateMVA( method );
           
           rar21 = reader21->GetRarity(method);
           
           fECII->SetT(rar21);
       
           count2++;
           outputTree2->Fill();
       
       }
   }

   sw2.Stop();
      std::cout << "--- End of event loop - event class 2: "; sw2.Print(); 
      
   TStopwatch sw3;
   sw3.Start();
      
   for( int treeNumber = 2; treeNumber < 4; ++treeNumber ) {
       

       if( treeNumber == 2 ){

           theTree3 = (TTree*)input->Get("TreeS3");
           std::cout << "--- Select signal sample" << std::endl;
          //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 2;
           
       } else if( treeNumber == 3 ){

           theTree3 = (TTree*)input->Get("TreeB3");
           std::cout << "--- Select background sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 3;
           
       }
       
       theTree3->SetBranchAddress("EventNumber", &EventNumber);
       theTree3->SetBranchAddress("PrimaryEnergy", &PriEnergy);
       theTree3->SetBranchAddress("eventID", &eventID);
       theTree3->SetBranchAddress("Pos_eX", &Pos_eX);
       theTree3->SetBranchAddress("Pos_eY", &Pos_eY);
       theTree3->SetBranchAddress("Pos_eZ", &Pos_eZ);
       theTree3->SetBranchAddress("Pos_pX", &Pos_pX);
       theTree3->SetBranchAddress("Pos_pY", &Pos_pY);
       theTree3->SetBranchAddress("Pos_pZ", &Pos_pZ);
       theTree3->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
       theTree3->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
       theTree3->SetBranchAddress( "PosX_Scat", &PosX_scat );
       theTree3->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
       theTree3->SetBranchAddress( "PosY_Scat", &PosY_scat );
       theTree3->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
       theTree3->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
       theTree3->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
       theTree3->SetBranchAddress( "Energy_Scat", &Energy_scat );
       theTree3->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
       theTree3->SetBranchAddress( "PosX_Abs", &PosX_abs );
       theTree3->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
       theTree3->SetBranchAddress( "PosY_Abs", &PosY_abs );
       theTree3->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
       theTree3->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
       theTree3->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
       theTree3->SetBranchAddress( "PosX_Sec", &PosX_sec );
       theTree3->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
       theTree3->SetBranchAddress( "PosY_Sec", &PosY_sec );
       theTree3->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
       theTree3->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
       theTree3->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
       theTree3->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
       theTree3->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
       theTree3->SetBranchAddress( "Energy_Abs", &Energy_abs );
       theTree3->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
       theTree3->SetBranchAddress("EnDiff", &EnDiff);
       theTree3->SetBranchAddress("EnergySecond_abs", &EnergySe);
       theTree3->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
       theTree3->SetBranchAddress( "EnergySum", &EnergySum );
       theTree3->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
       theTree3->SetBranchAddress( "Multiplicity", &Multiplicity );
       theTree3->SetBranchAddress( "DiffPosition", &DiffPosition );
       theTree3->SetBranchAddress( "ECII", &fECII );
       theTree3->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           
       std::cout << "--- Processing: " << theTree3->GetEntries() << " events" << std::endl;
 
       Int_t nEvent3 = theTree3->GetEntries();
          
       for (Long64_t ievt=0; ievt<nEvent3; ievt++) {
           
           
           if (ievt%10000 == 0){
    
               std::cout << "--- ... Processing events: " << ievt << std::endl;
           }
    
           theTree3->GetEntry(ievt);
           // get the classifier's response for each of the signal/background classifications
           classifier31 = reader31->EvaluateMVA( method );
           
           rar31 = reader31->GetRarity(method);
           
           fECII->SetT(rar31);
           
           count3++;
           outputTree3->Fill();
       }
      
   }
   sw3.Stop();
   std::cout << "--- End of event loop - event class 3: "; sw3.Print();   
  
   TStopwatch sw4;
   sw4.Start();
      
   for( int treeNumber = 4; treeNumber < 6; ++treeNumber ) {
       
       
     
       if( treeNumber == 4 ){
           
           theTree4 = (TTree*)input->Get("TreeS4");
           std::cout << "--- Select signal sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 4;
           
       } else if( treeNumber == 5 ){
           
           theTree4 = (TTree*)input->Get("TreeB4");
           std::cout << "--- Select background sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 5;
           
       }
       
       theTree4->SetBranchAddress("EventNumber", &EventNumber);
       theTree4->SetBranchAddress("PrimaryEnergy", &PriEnergy);
       theTree4->SetBranchAddress("eventID", &eventID);
       theTree4->SetBranchAddress("Pos_eX", &Pos_eX);
       theTree4->SetBranchAddress("Pos_eY", &Pos_eY);
       theTree4->SetBranchAddress("Pos_eZ", &Pos_eZ);
       theTree4->SetBranchAddress("Pos_pX", &Pos_pX);
       theTree4->SetBranchAddress("Pos_pY", &Pos_pY);
       theTree4->SetBranchAddress("Pos_pZ", &Pos_pZ);
       theTree4->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
       theTree4->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
       theTree4->SetBranchAddress( "PosX_Scat", &PosX_scat );
       theTree4->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
       theTree4->SetBranchAddress( "PosY_Scat", &PosY_scat );
       theTree4->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
       theTree4->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
       theTree4->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
       theTree4->SetBranchAddress( "Energy_Scat", &Energy_scat );
       theTree4->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
       theTree4->SetBranchAddress( "PosX_Abs", &PosX_abs );
       theTree4->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
       theTree4->SetBranchAddress( "PosY_Abs", &PosY_abs );
       theTree4->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
       theTree4->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
       theTree4->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
       theTree4->SetBranchAddress( "PosX_Sec", &PosX_sec );
       theTree4->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
       theTree4->SetBranchAddress( "PosY_Sec", &PosY_sec );
       theTree4->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
       theTree4->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
       theTree4->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
       theTree4->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
       theTree4->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
       theTree4->SetBranchAddress( "Energy_Abs", &Energy_abs );
       theTree4->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
       theTree4->SetBranchAddress("EnDiff", &EnDiff);
       theTree4->SetBranchAddress("EnergySecond_abs", &EnergySe);
       theTree4->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
       theTree4->SetBranchAddress( "EnergySum", &EnergySum );
       theTree4->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
       theTree4->SetBranchAddress( "Multiplicity", &Multiplicity );
       theTree4->SetBranchAddress( "DiffPosition", &DiffPosition );
       theTree4->SetBranchAddress( "ECII", &fECII );
       theTree4->SetBranchAddress( "AngularDistribution", &AngularDistribution );
      
       std::cout << "--- Processing: " << theTree4->GetEntries() << " events" << std::endl;
      
      
       Int_t nEvent4 = theTree4->GetEntries();
           
       for (Long64_t ievt=0; ievt<nEvent4; ievt++) {
           

           if (ievt%10000 == 0){

              std::cout << "--- ... Processing events: " << ievt << std::endl;
           }
    
           theTree4->GetEntry(ievt);
            
           classifier41 = reader41->EvaluateMVA( method );
           
           rar41 = reader41->GetRarity(method);
           
           fECII->SetT(rar41);
           
           count4++;
           outputTree4->Fill();
       }
   }
 
// get elapsed time
   sw4.Stop();
   std::cout << "--- End of event loop - event class 4: "; sw4.Print();
       
   TStopwatch sw5;
   sw5.Start();
      
   for( int treeNumber = 6; treeNumber < 8; ++treeNumber ) {
       
       if( treeNumber == 6 ){

           theTree5 = (TTree*)input->Get("TreeS5");
           std::cout << "--- Select signal sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 6;
           
       } else if( treeNumber == 7 ){
           
           theTree5 = (TTree*)input->Get("TreeB5");
           std::cout << "--- Select background sample" << std::endl;
           //theTree->SetBranchAddress( "weight", &weight );
           weight = 1;
           classID = 7;
           
       }

       theTree5->SetBranchAddress("EventNumber", &EventNumber);
       theTree5->SetBranchAddress("PrimaryEnergy", &PriEnergy);
       theTree5->SetBranchAddress("eventID", &eventID);
       theTree5->SetBranchAddress("Pos_eX", &Pos_eX);
       theTree5->SetBranchAddress("Pos_eY", &Pos_eY);
       theTree5->SetBranchAddress("Pos_eZ", &Pos_eZ);
       theTree5->SetBranchAddress("Pos_pX", &Pos_pX);
       theTree5->SetBranchAddress("Pos_pY", &Pos_pY);
       theTree5->SetBranchAddress("Pos_pZ", &Pos_pZ);
       theTree5->SetBranchAddress("RealEnergy_e", &RealEnergy_e);
       theTree5->SetBranchAddress("RealEnergy_p", &RealEnergy_p);
       theTree5->SetBranchAddress( "PosX_Scat", &PosX_scat );
       theTree5->SetBranchAddress("PosX_ScatUnc", &PosX_scatunc);
       theTree5->SetBranchAddress( "PosY_Scat", &PosY_scat );
       theTree5->SetBranchAddress("PosY_ScatUnc", &PosY_scatunc);
       theTree5->SetBranchAddress( "PosZ_Scat", &PosZ_scat );
       theTree5->SetBranchAddress("PosZ_ScatUnc", &PosZ_scatunc);
       theTree5->SetBranchAddress( "Energy_Scat", &Energy_scat );
       theTree5->SetBranchAddress("Energy_ScatUnc", &Energy_scatunc);
       theTree5->SetBranchAddress( "PosX_Abs", &PosX_abs );
       theTree5->SetBranchAddress("PosX_AbsUnc", &PosX_absunc);
       theTree5->SetBranchAddress( "PosY_Abs", &PosY_abs );
       theTree5->SetBranchAddress("PosY_AbsUnc", &PosY_absunc);
       theTree5->SetBranchAddress( "PosZ_Abs", &PosZ_abs );
       theTree5->SetBranchAddress("PosZ_AbsUnc", &PosZ_absunc);
       theTree5->SetBranchAddress( "PosX_Sec", &PosX_sec );
       theTree5->SetBranchAddress("PosX_SecUnc", &PosX_secunc);
       theTree5->SetBranchAddress( "PosY_Sec", &PosY_sec );
       theTree5->SetBranchAddress("PosY_SecUnc", &PosY_secunc);
       theTree5->SetBranchAddress( "PosZ_Sec", &PosZ_sec );
       theTree5->SetBranchAddress("PosZ_SecUnc", &PosZ_secunc);
       theTree5->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
       theTree5->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
       theTree5->SetBranchAddress( "Energy_Abs", &Energy_abs );
       theTree5->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
       theTree5->SetBranchAddress("EnDiff", &EnDiff);
       theTree5->SetBranchAddress("EnergySecond_abs", &EnergySe);
       theTree5->SetBranchAddress("EnergySecond_absUnc", &EnergySeunc);
       theTree5->SetBranchAddress( "EnergySum", &EnergySum );
       theTree5->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
       theTree5->SetBranchAddress( "Multiplicity", &Multiplicity );
       theTree5->SetBranchAddress( "DiffPosition", &DiffPosition );
       theTree5->SetBranchAddress( "ECII", &fECII );
       theTree5->SetBranchAddress( "AngularDistribution", &AngularDistribution );
       
       std::cout << "--- Processing: " << theTree5->GetEntries() << " events" << std::endl;

       Int_t nEvent5 = theTree5->GetEntries();
           
       for (Long64_t ievt=0; ievt<nEvent5; ievt++) {
           
           if (ievt%10000 == 0){
               
              std::cout << "--- ... Processing events: " << ievt << std::endl;
           }
    
           theTree5->GetEntry(ievt);
        
           classifier51 = reader51->EvaluateMVA( method );
           
           rar51 = reader51->GetRarity(method);
           
           fECII->SetT(rar51);
          
           count5++;
           outputTree5->Fill();
       }
   }

   // get elapsed time
   sw5.Stop();
   std::cout << "--- End of event loop - event class 5: "; sw5.Print();
   
   std::cout << "Total No. of event classes : " << " \t " << " S = 2 : " << count2 << " \t " << " S = 3 : " << count3 << " \t " << " S = 4 : " << count4 << " \t " << " S = 5 : " << count5 << endl; 
   
   input->Close();
   // write output tree
/*   outputTree->SetDirectory(outputFile);
     outputTree->Write(); */
   outputFile->Write();
   outputFile->Close();
   std::cout << "--- Created root file: \"" << outfileName.Data() << "\" containing the MVA output histograms" << std::endl;
   delete reader21;
   
   delete reader31;
  
   delete reader41;
  
   delete reader51;
  
   
   std::cout << "==> Application of readers is done! combined tree created" << std::endl << std::endl;
}
////////////////////////////////////////////////////////////////////////////////////////////

void EventClassification(){
   // ----------------------------------------------------------------------------------------
   // Run all
   // ----------------------------------------------------------------------------------------
    
   cout << "========================" << endl;
   cout << "--- Training" << endl;
   Training();
   cout << endl;
   
/*

   cout << "========================" << endl;
   cout << "--- CrossValidation" << endl;
   TMVACrossValidation();
   cout << endl;
   
   */

   cout << "========================" << endl;
   cout << "--- Application & create combined tree" << endl;
   ApplicationCreateCombinedTree();
   cout << endl;



}
//////////////////////////////////////////

int main( int argc, char** argv ) {
    
   EventClassification();
}
