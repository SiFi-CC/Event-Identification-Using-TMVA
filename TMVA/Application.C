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

#include <TLorentzVector.h>

#include "TMVA/TMVARegGui.h"

using namespace std;
using namespace TMVA;

// -----------------------------------------------------------------------------------------
// Genetic Algorithm Fitness definition (event class 2 clusters)
// -----------------------------------------------------------------------------------------
//
class MyFitness2 : public IFitterTarget {
public:
   // constructor
   MyFitness2( TChain* _chain ) : IFitterTarget() {
       
       chain2 = _chain;
       hSignal2 = new TH1F("hsignal2","hsignal2",100,-1,1);
       hFP2 = new TH1F("hfp2","hfp2",100,-1,1);
       hTP2 = new TH1F("htp2","htp2",100,-1,1);
       
       TString cutsAndWeightSignal  = "(classID==0)";
       nSignal2 = chain2->Draw("Entry$/Entries$>>hsignal2",cutsAndWeightSignal,"goff");
       
       weightsSignal2 = hSignal2->Integral();
   }

// the output of this function will be minimized
   Double_t EstimatorFunction( std::vector<Double_t> & factors2 ){
       
       
       TString cutsAndWeightTruePositive  = Form("((classID==0) && cls21>%f)",factors2.at(0));
       
       TString cutsAndWeightFalsePositive = Form("((classID > 0) && (classID < 2) && cls21>%f )",factors2.at(0));
       
      // Entry$/Entries$ just draws something reasonable. Could in principle anything
       Float_t nTP = chain2->Draw("Entry$/Entries$>>htp2",cutsAndWeightTruePositive,"goff");
       Float_t nFP = chain2->Draw("Entry$/Entries$>>hfp2",cutsAndWeightFalsePositive,"goff");
      
       weightsTruePositive2 = hTP2->Integral();
       weightsFalsePositive2 = hFP2->Integral();
       
       efficiency2 = 0;
       if( weightsSignal2 > 0 )  efficiency2 = weightsTruePositive2/weightsSignal2;
       
       purity2 = 0;
       if( weightsTruePositive2+weightsFalsePositive2 > 0 )  purity2 = weightsTruePositive2/(weightsTruePositive2+weightsFalsePositive2);
       
       Float_t effTimesPur = efficiency2*purity2;
       
       Float_t toMinimize = std::numeric_limits<float>::max(); // set to the highest existing number
       
       if( effTimesPur > 0 ) {
           // if larger than 0, take 1/x. This is the value to minimize
           toMinimize = 1./(effTimesPur); // we want to minimize 1/efficiency*purity
       }
      // Print();
      return toMinimize;
   }
   
   void Print(){
       
       std::cout << std::endl;
       std::cout << "======================" << std::endl
       << "Efficiency2 : " << efficiency2 << std::endl
       << "Purity2     : " << purity2 << std::endl << std::endl
       << "Signal weights2        : " << weightsSignal2 << std::endl
       << "True positive weights2 : " << weightsTruePositive2 << std::endl
       << "False positive weights2: " << weightsFalsePositive2 << std::endl;
       
       
  
   }
   
   Float_t nSignal2;
   Float_t efficiency2;
   Float_t purity2;
   Float_t weightsTruePositive2;
   Float_t weightsFalsePositive2;
   Float_t weightsSignal2;
  
private:
    
    TChain* chain2;
    TH1F* hSignal2;
    TH1F* hFP2;
    TH1F* hTP2;
    
};
// -----------------------------------------------------------------------------------------
// Genetic Algorithm Fitness definition (event class 3 clusters)
// -----------------------------------------------------------------------------------------
//

class MyFitness3 : public IFitterTarget {
public:
   // constructor
   MyFitness3( TChain* _chain ) : IFitterTarget() {
       
       chain3 = _chain;
       hSignal3 = new TH1F("hsignal3","hsignal3",100,-1,1);
       hFP3 = new TH1F("hfp3","hfp3",100,-1,1);
       hTP3 = new TH1F("htp3","htp3",100,-1,1);
       
       TString cutsAndWeightSignal  = "(classID==2)";
       nSignal3 = chain3->Draw("Entry$/Entries$>>hsignal3",cutsAndWeightSignal,"goff");
      
       weightsSignal3 = hSignal3->Integral();
   }

   // the output of this function will be minimized
   Double_t EstimatorFunction( std::vector<Double_t> & factors3 ){
      
       TString cutsAndWeightTruePositive  = Form("((classID==2) && cls31>%f )",factors3.at(0));
       
       TString cutsAndWeightFalsePositive = Form("((classID > 2) && (classID < 4) && cls31>%f )",factors3.at(0));
       
      
       Float_t nTP = chain3->Draw("Entry$/Entries$>>htp3",cutsAndWeightTruePositive,"goff");
       Float_t nFP = chain3->Draw("Entry$/Entries$>>hfp3",cutsAndWeightFalsePositive,"goff");
       
       weightsTruePositive3 = hTP3->Integral();
       weightsFalsePositive3 = hFP3->Integral();
       
       efficiency3 = 0;
       if( weightsSignal3 > 0 )  efficiency3 = weightsTruePositive3/weightsSignal3;
       
       purity3 = 0;
       if( weightsTruePositive3+weightsFalsePositive3 > 0 )  purity3 = weightsTruePositive3/(weightsTruePositive3+weightsFalsePositive3);
       
       Float_t effTimesPur = efficiency3*purity3;
       
       Float_t toMinimize = std::numeric_limits<float>::max(); // set to the highest existing number
       
       if( effTimesPur > 0 ) {
           
           toMinimize = 1./(effTimesPur); 
       }
      
      return toMinimize;
   }
   
   void Print(){
       
       std::cout << std::endl;
       std::cout << "======================" << std::endl
       << "Efficiency3 : " << efficiency3 << std::endl
       << "Purity3     : " << purity3 << std::endl << std::endl
       << "Signal weights3        : " << weightsSignal3 << std::endl
       << "True positive weights3 : " << weightsTruePositive3 << std::endl
       << "False positive weights3: " << weightsFalsePositive3 << std::endl;
       
       
  
   }
   
   Float_t nSignal3;
   Float_t efficiency3;
   Float_t purity3;
   Float_t weightsTruePositive3;
   Float_t weightsFalsePositive3;
   Float_t weightsSignal3;
  
private:
    
    TChain* chain3;
    TH1F* hSignal3;
    TH1F* hFP3;
    TH1F* hTP3;
    
};

// -----------------------------------------------------------------------------------------
// Genetic Algorithm Fitness definition (event class 4 clusters)
// -----------------------------------------------------------------------------------------
//

class MyFitness4 : public IFitterTarget {
public:
   // constructor
  
   MyFitness4( TChain* _chain ) : IFitterTarget() {
       
       chain4 = _chain;
       hSignal4 = new TH1F("hsignal4","hsignal4",100,-1,1);
       hFP4 = new TH1F("hfp4","hfp4",100,-1,1);
       hTP4 = new TH1F("htp4","htp4",100,-1,1);
     
       TString cutsAndWeightSignal  = "(classID==4)";
       nSignal4 = chain4->Draw("Entry$/Entries$>>hsignal4",cutsAndWeightSignal,"goff");
      
       weightsSignal4 = hSignal4->Integral();
   }

   // the output of this function will be minimized
 
   Double_t EstimatorFunction( std::vector<Double_t> & factors4 ){
       
       
       TString cutsAndWeightTruePositive  = Form("((classID==4) && cls41>%f )",factors4.at(0));
       
       TString cutsAndWeightFalsePositive = Form("((classID > 4) && (classID < 6) && cls41>%f )",factors4.at(0));
      
       Float_t nTP = chain4->Draw("Entry$/Entries$>>htp4",cutsAndWeightTruePositive,"goff");
       Float_t nFP = chain4->Draw("Entry$/Entries$>>hfp4",cutsAndWeightFalsePositive,"goff");
      
       weightsTruePositive4 = hTP4->Integral();
       weightsFalsePositive4 = hFP4->Integral();
       
       efficiency4 = 0;
       if( weightsSignal4 > 0 )  efficiency4 = weightsTruePositive4/weightsSignal4;
       
       purity4 = 0;
       if( weightsTruePositive4+weightsFalsePositive4 > 0 )  purity4 = weightsTruePositive4/(weightsTruePositive4+weightsFalsePositive4);
       
       Float_t effTimesPur = efficiency4*purity4;
       
       Float_t toMinimize = std::numeric_limits<float>::max(); 
       
       if( effTimesPur > 0 ) {
           
           toMinimize = 1./(effTimesPur); 
       }
      
      return toMinimize;
   }
   void Print(){
       
       std::cout << "======================" << std::endl
       << "Efficiency4 : " << efficiency4 << std::endl
       << "Purity4     : " << purity4 << std::endl << std::endl
       << "Signal weights4        : " << weightsSignal4 << std::endl
       << "True positive weights4 : " << weightsTruePositive4 << std::endl
       << "False positive weights4: " << weightsFalsePositive4 << std::endl;
       
   }
   
   
   Float_t nSignal4;
   Float_t efficiency4;
   Float_t purity4;
   Float_t weightsTruePositive4;
   Float_t weightsFalsePositive4;
   Float_t weightsSignal4;

private:
    
    
    TChain* chain4;
    TH1F* hSignal4;
    TH1F* hFP4;
    TH1F* hTP4;
};
// -----------------------------------------------------------------------------------------
// Genetic Algorithm Fitness definition (event class 5 clusters)
// -----------------------------------------------------------------------------------------
//
class MyFitness5 : public IFitterTarget {
public:
   // constructor
  
   MyFitness5( TChain* _chain ) : IFitterTarget() {
       
       chain5 = _chain;
       hSignal5 = new TH1F("hsignal5","hsignal5",100,-1,1);
       hFP5 = new TH1F("hfp5","hfp5",100,-1,1);
       hTP5 = new TH1F("htp5","htp5",100,-1,1);
      
       TString cutsAndWeightSignal  = "(classID==6)";
       nSignal5 = chain5->Draw("Entry$/Entries$>>hsignal5",cutsAndWeightSignal,"goff");
       
       weightsSignal5 = hSignal5->Integral();
   }

 
   Double_t EstimatorFunction( std::vector<Double_t> & factors5 ){
       
       
       TString cutsAndWeightTruePositive  = Form("((classID==6) && cls51>%f)",factors5.at(0));
       
       TString cutsAndWeightFalsePositive = Form("((classID > 6) && (classID < 8) && cls51>%f)",factors5.at(0));
       
       Float_t nTP = chain5->Draw("Entry$/Entries$>>htp5",cutsAndWeightTruePositive,"goff");
       Float_t nFP = chain5->Draw("Entry$/Entries$>>hfp5",cutsAndWeightFalsePositive,"goff");
       
       weightsTruePositive5 = hTP5->Integral();
       weightsFalsePositive5 = hFP5->Integral();
       
       efficiency5 = 0;
       if( weightsSignal5 > 0 )  efficiency5 = weightsTruePositive5/weightsSignal5;
       
       purity5 = 0;
       if( weightsTruePositive5+weightsFalsePositive5 > 0 )  purity5 = weightsTruePositive5/(weightsTruePositive5+weightsFalsePositive5);
       
       Float_t effTimesPur = efficiency5*purity5;
       
       Float_t toMinimize = std::numeric_limits<float>::max(); 
       
       if( effTimesPur > 0 ) {
           
           toMinimize = 1./(effTimesPur);
       }
      
      return toMinimize;
   }
   void Print(){
       
       std::cout << "======================" << std::endl
       << "Efficiency5 : " << efficiency5 << std::endl
       << "Purity5     : " << purity5 << std::endl << std::endl
       << "Signal weights5        : " << weightsSignal5 << std::endl
       << "True positive weights5 : " << weightsTruePositive5 << std::endl
       << "False positive weights5: " << weightsFalsePositive5 << std::endl;
       
   }
   
   
   Float_t nSignal5;
   Float_t efficiency5;
   Float_t purity5;
   Float_t weightsTruePositive5;
   Float_t weightsFalsePositive5;
   Float_t weightsSignal5;

private:
    
    
    TChain* chain5;
    TH1F* hSignal5;
    TH1F* hFP5;
    TH1F* hTP5;
};

// ----------------------------------------------------------------------------------------------
// Call of Genetic algorithm and then the function applies the optimized cuts on each event class
// ----------------------------------------------------------------------------------------------
//// One can either use the optimized cuts offered by the training to get the maximum significance
//// or pass the dataset to Genetic Fitter to obtain the optimal cuts  
void MaxiSignificanceCutApp(){
    
    std::vector<Double_t> cut;
    // define all the parameters by their minimum and maximum value
    
    
    vector<Interval*> ranges2;
    ranges2.push_back( new Interval(-1,1) );
    
    vector<Interval*> ranges3;
    ranges3.push_back( new Interval(-1,1) ); 
    
    vector<Interval*> ranges4;
    ranges4.push_back( new Interval(-1,1) );
    
    vector<Interval*> ranges5;
    ranges5.push_back( new Interval(-1,1) );
   
    std::cout << "Classifier ranges (defined by the user)" << std::endl;
    for( std::vector<Interval*>::iterator it = ranges2.begin(); it != ranges2.end(); it++ ){
        
        std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
    }

    TChain* chain2 = new TChain("SigBkg2");
    //chain2->Add("tmva__applied-s2s3s4s5CV.root");
    chain2->Add("tmva__applied-s2s3s4s5ST-9var.root");
    
    IFitterTarget* myFitness2 = new MyFitness2( chain2 );
    // mind: big population sizes will help in searching the domain space of the solution
    // but you have to weight this out to the number of generations
    // the extreme case of 1 generation and population size n is equal to
    // a Monte Carlo calculation with n entries
    const TString name2( "multipleBackgroundGA" );
    const TString opts2( "PopSize=100:Steps=30:Cycles=3" );
    GeneticFitter mg2( *myFitness2, name2, ranges2, opts2);
    // mg.SetParameters( 4, 30, 200, 10,5, 0.95, 0.001 );
    std::vector<Double_t> result2;
    Double_t estimator2 = mg2.Run(result2);
    dynamic_cast<MyFitness2*>(myFitness2)->Print();
    std::cout << std::endl;
    int n2 = 0;
    for( std::vector<Double_t>::iterator it = result2.begin(); it<result2.end(); it++ ){
        
        std::cout << "  cutValue[" << n2 << "] = " << (*it) << ";"<< std::endl;
        cut.push_back(*it);
        n2++;
    }
    
    

/////////////////////////////////////////////////////////////////////////  

    std::cout << "Classifier ranges (defined by the user)" << std::endl;
    for( std::vector<Interval*>::iterator it = ranges3.begin(); it != ranges3.end(); it++ ){
        
        std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
    }

    TChain* chain3 = new TChain("SigBkg3");
    //chain3->Add("tmva__applied-s2s3s4s5CV.root");
    chain3->Add("tmva__applied-s2s3s4s5ST-9var.root");
    
    IFitterTarget* myFitness3 = new MyFitness3( chain3 );
    
    const TString name3( "multipleBackgroundGA" );
    const TString opts3( "PopSize=100:Steps=30:Cycles=3:SaveBestCycle=10:SaveBestGen=1" );
    GeneticFitter mg3( *myFitness3, name3, ranges3, opts3);
    // mg.SetParameters( 4, 30, 200, 10,5, 0.95, 0.001 );
    std::vector<Double_t> result3;
    Double_t estimator3 = mg3.Run(result3);
    dynamic_cast<MyFitness3*>(myFitness3)->Print();
    std::cout << std::endl;
    int n3 = 0;
    for( std::vector<Double_t>::iterator it = result3.begin(); it<result3.end(); it++ ){
        
        std::cout << "  cutValue[" << n3 << "] = " << (*it) << ";"<< std::endl;
        cut.push_back(*it);
        n3++;
    }
    
    
//////////////////////////////////////////////////////////////////////////////

    std::cout << "Classifier ranges (defined by the user)" << std::endl;
    for( std::vector<Interval*>::iterator it = ranges4.begin(); it != ranges4.end(); it++ ){
        
        std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
    }

    TChain* chain4 = new TChain("SigBkg4");
    //chain4->Add("tmva__applied-s2s3s4s5CV.root");
    chain4->Add("tmva__applied-s2s3s4s5ST-9var.root");
    
    IFitterTarget* myFitness4 = new MyFitness4( chain4 );
 
    const TString name4( "multipleBackgroundGA" );
    const TString opts4( "PopSize=100:Steps=30:Cycles=3" );
    GeneticFitter mg4( *myFitness4, name4, ranges4, opts4);
    // mg.SetParameters( 4, 30, 200, 10,5, 0.95, 0.001 );
    std::vector<Double_t> result4;
    Double_t estimator4 = mg4.Run(result4);
    dynamic_cast<MyFitness4*>(myFitness4)->Print();
    std::cout << std::endl;
    int n4 = 0;
    for( std::vector<Double_t>::iterator it = result4.begin(); it<result4.end(); it++ ){
        
        std::cout << "  cutValue[" << n4 << "] = " << (*it) << ";"<< std::endl;
        cut.push_back(*it);
        n4++;
    }

//////////////////////////////////////////////////////////////////////////////

    std::cout << "Classifier ranges (defined by the user)" << std::endl;
    for( std::vector<Interval*>::iterator it = ranges5.begin(); it != ranges5.end(); it++ ){
        
        std::cout << " range: " << (*it)->GetMin() << "   " << (*it)->GetMax() << std::endl;
    }

    TChain* chain5 = new TChain("SigBkg5");
    //chain5->Add("tmva__applied-s2s3s4s5CV.root");
    chain5->Add("tmva__applied-s2s3s4s5ST-9var.root");
    
    IFitterTarget* myFitness5 = new MyFitness5( chain5 );
   
    const TString name5( "multipleBackgroundGA" );
    const TString opts5( "PopSize=100:Steps=30:Cycles=3" );
    GeneticFitter mg5( *myFitness5, name5, ranges5, opts5);
    // mg.SetParameters( 4, 30, 200, 10,5, 0.95, 0.001 );
    std::vector<Double_t> result5;
    Double_t estimator5 = mg5.Run(result5);
    dynamic_cast<MyFitness5*>(myFitness5)->Print();
    std::cout << std::endl;
    int n5 = 0;
    for( std::vector<Double_t>::iterator it = result5.begin(); it<result5.end(); it++ ){
        
        std::cout << "  cutValue[" << n5 << "] = " << (*it) << ";"<< std::endl;
        cut.push_back(*it);
        n5++;
    }
       
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

    //TString outfileName( "tmvapp-s2s3s4s5-RecoEnCVGA.root" );
    //TString outfileName( "tmvapp-s2s3s4s5-STRegInAna.root" );
    TString outfileName( "tmvapp-s2s3s4s5-ST-9var.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
    TTree* outputTree2 = new TTree("SigBkg2","SigBkg tree");
    TTree* outputTree3 = new TTree("SigBkg3","SigBkg tree");
    TTree* outputTree4 = new TTree("SigBkg4","SigBkg tree");
    TTree* outputTree5 = new TTree("SigBkg5","SigBkg tree");

    TTree* outputTree = new TTree("Tree","Tree");
    
    Float_t Pos_eX, Pos_eY, Pos_eZ, Pos_pX, Pos_pY, Pos_pZ, RealEnergy_e, RealEnergy_p; 
    
    Float_t PosX_scat, PosY_scat, PosZ_scat,Energy_scat, ReEnergy_scat;
    Float_t PosX_abs, PosY_abs, PosZ_abs, Energy_abs, ReEnergy_abs, EnDiff;
    Float_t PosX_sec, PosY_sec, PosZ_sec;
    Float_t PosX_secunc, PosY_secunc, PosZ_secunc;
    
    Float_t PosX_scatunc, PosY_scatunc, PosZ_scatunc,Energy_scatunc;
    Float_t PosX_absunc, PosY_absunc, PosZ_absunc, Energy_absunc;
    Float_t EnergyCluster_absunc,EnergySe, EnergySeunc;
    
    Int_t Multiplicity;
    Float_t DiffEnergy, RatioEnergy, DiffPosition, EnergySum, PriEnergy, EnergySumunc;
    Float_t EnergyCluster_abs, AngularDistribution;
    Int_t count = 0,counts = 0,countb = 0,countdup = 0,countrb = 0,countbs = 0;
    Int_t countrb2 = 0, countrb3 = 0, countrb4 = 0, countrb5 = 0;
    Int_t countbs2 = 0, countbs3 = 0, countbs4 = 0, countbs5 = 0;
    Int_t EventNumber;
    
    Int_t count20 = 0, count21 = 0, count22 = 0, count23 = 0;
    Int_t count30 = 0, count31 = 0, count32 = 0, count33 = 0, count34 = 0;
    Int_t count40 = 0, count41 = 0, count42 = 0, count43 = 0, count44 = 0;
    Int_t count50 = 0, count51 = 0, count52 = 0, count53 = 0, count54 = 0;
    
    Int_t countRB3 = 0, countRB4 = 0, countRB33 = 0, countRB44 = 0,countRB55 = 0, countRB5 = 0;
    Int_t countBS3 = 0, countBS4 = 0, countBS33 = 0, countBS44 = 0, countBS55 = 0, countBS5 = 0;
    Int_t count311 = 0, count411 = 0, count511 = 0;
   
    Int_t countp3 = 0, countp31 = 0, countp32 = 0, countunique3 = 0;
    Int_t countp4 = 0, countp4D1 = 0, countp4D2 = 0, countp41 = 0, countp42 = 0, countp43 = 0, countp44 = 0, countp45 = 0, countp46 = 0, countp47 = 0, countunique4 = 0;
    Int_t countp5 = 0, countp51 = 0, countp52 = 0, countp53 = 0, countp54 = 0, countunique5 = 0;
    
    Int_t count2non = 0, count3non = 0, count4non = 0, count5non = 0;
    Int_t count2bs = 0, count3bs = 0, count4bs = 0, count5bs = 0;
    Int_t count2bb = 0, count3bb = 0, count4bb = 0, count5bb = 0;
    
    Int_t   classID;
    Float_t weight = 1.f;
    Float_t classifier21, classifier22, classifier23, classifier31, classifier32, classifier33, classifier34, classifier41, classifier42, classifier43, classifier44, classifier51, classifier52, classifier53, classifier54;
    
    Float_t ReES,CPriEn, ReESunc, sigma;
    
    Int_t eventID;
    
    Float_t prob21, prob22, prob23, prob31, prob32, prob33, prob34, prob41, prob42, prob43, prob44, prob51, prob52, prob53, prob54;  
    
    Float_t pSig21, pSig22, pSig23, pSig31, pSig32, pSig33, pSig34, pSig41, pSig42, pSig43, pSig44, pSig51, pSig52, pSig53, pSig54;
    
    Float_t rar21, rar22, rar23, rar31, rar32, rar33, rar34, rar41, rar42, rar43, rar44, rar51, rar52, rar53, rar54;
    
    TLorentzVector *fECII = new TLorentzVector();
    
    gROOT->SetBatch(kTRUE);

    const Int_t nx = 2;
    const char *type1[nx] = {"True events","Background events"};
    
    TH1F *tt = new TH1F("Total number events( Duplicate events Excluded )","Total number events( Duplicate events Excluded )", 2,0,2);
    tt->GetYaxis()->SetTitle("count");
    tt->GetXaxis()->SetTitle(" ");

    const Int_t n = 8;
    const char *type2[n] = {"2 clusters(T)","3 clusters(T)","4 clusters(T)","5 clusters(T)","2 clusters(B)","3 clusters(B)","4 clusters(B)","5 clusters(B)"};
    
    TH1F *ttt = new TH1F("Identified events (Duplicates Excluded)","Identified events(Duplicates Excluded)", 7,0,7);
    ttt->GetYaxis()->SetTitle("count");
    ttt->GetXaxis()->SetTitle(" ");

   outputTree2->Branch("EventNumber", &EventNumber);
   outputTree2->Branch("PrimaryEnergy", &PriEnergy);
   outputTree2->Branch("ReEnergySum", &ReES);
   outputTree2->Branch("ECII", &fECII);
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
   outputTree2->Branch("EnergySum", &EnergySum);
   outputTree2->Branch("EnergySumUnc", &EnergySumunc);
   outputTree2->Branch("DiffPosition", &DiffPosition);
   outputTree2->Branch("AngularDistribution", &AngularDistribution);
   outputTree2->Branch("classID", &classID, "classID/I");
   outputTree2->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree2->Branch("cls21", &classifier21, "cls21/F");
   outputTree2->Branch("rarity21", &rar21, "rarity21/F");
  

   outputTree3->Branch("EventNumber", &EventNumber);
   outputTree3->Branch("PrimaryEnergy", &PriEnergy);
   outputTree3->Branch("ReEnergySum", &ReES);
   outputTree3->Branch("ECII", &fECII);
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
   outputTree3->Branch("EnergySum", &EnergySum);
   outputTree3->Branch("EnergySumUnc", &EnergySumunc);
   outputTree3->Branch("DiffPosition", &DiffPosition);
   outputTree3->Branch("AngularDistribution", &AngularDistribution);
   outputTree3->Branch("classID", &classID, "classID/I");
   outputTree3->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree3->Branch("cls31", &classifier31, "cls31/F");
   outputTree3->Branch("rarity31", &rar31, "rarity31/F");

   
   outputTree4->Branch("EventNumber", &EventNumber);
   outputTree4->Branch("PrimaryEnergy", &PriEnergy);
   outputTree4->Branch("ReEnergySum", &ReES);
   outputTree4->Branch("ECII", &fECII);
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
   outputTree4->Branch("EnergySum", &EnergySum);
   outputTree4->Branch("EnergySumUnc", &EnergySumunc);
   outputTree4->Branch("DiffPosition", &DiffPosition);
   outputTree4->Branch("AngularDistribution", &AngularDistribution);
   outputTree4->Branch("classID", &classID, "classID/I");
   outputTree4->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree4->Branch("cls41", &classifier41, "cls41/F");
   outputTree4->Branch("rarity41", &rar41, "rarity41/F");

   
   outputTree5->Branch("EventNumber", &EventNumber);
   outputTree5->Branch("PrimaryEnergy", &PriEnergy);
   outputTree5->Branch("ReEnergySum", &ReES);
   outputTree5->Branch("ECII", &fECII);
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
   outputTree5->Branch("EnergySum", &EnergySum);
   outputTree5->Branch("EnergySumUnc", &EnergySumunc);
   outputTree5->Branch("DiffPosition", &DiffPosition);
   outputTree5->Branch("AngularDistribution", &AngularDistribution);
   outputTree5->Branch("classID", &classID, "classID/I");
   outputTree5->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   outputTree5->Branch("cls51", &classifier51, "cls51/F");
   outputTree5->Branch("rarity51", &rar51, "rarity51/F");

   outputTree->Branch("EventNumber", &EventNumber);
   outputTree->Branch("PrimaryEnergy", &PriEnergy);
   outputTree->Branch("ReEnergySum", &ReES);
   outputTree->Branch("ECII", &fECII);
   outputTree->Branch("PosX_Scat", &PosX_scat);
   outputTree->Branch("PosX_ScatUnc", &PosX_scatunc);
   outputTree->Branch("PosY_Scat", &PosY_scat);
   outputTree->Branch("PosY_ScatUnc", &PosY_scatunc);
   outputTree->Branch("PosZ_Scat", &PosZ_scat);
   outputTree->Branch("PosZ_ScatUnc", &PosZ_scatunc);
   outputTree->Branch("Energy_Scat", &Energy_scat);
   outputTree->Branch("Energy_ScatUnc", &Energy_scatunc);
   outputTree->Branch("PosX_Abs", &PosX_abs);
   outputTree->Branch("PosX_AbsUnc", &PosX_absunc);
   outputTree->Branch("PosY_Abs", &PosY_abs);
   outputTree->Branch("PosY_AbsUnc", &PosY_absunc);
   outputTree->Branch("PosZ_Abs", &PosZ_abs);
   outputTree->Branch("PosZ_AbsUnc", &PosZ_absunc);
   outputTree->Branch("EnergyCluster_abs", &EnergyCluster_abs);
   outputTree->Branch("EnergyCluster_absUnc", &EnergyCluster_absunc );
   outputTree->Branch("Energy_Abs", &Energy_abs);
   outputTree->Branch("Energy_AbsUnc", &Energy_absunc);
   outputTree->Branch("EnergySum", &EnergySum);
   outputTree->Branch("EnergySumUnc", &EnergySumunc);
   outputTree->Branch("DiffPosition", &DiffPosition);
   outputTree->Branch("AngularDistribution", &AngularDistribution);
   outputTree->Branch("classID", &classID, "classID/I");
   outputTree->Branch("Multiplicity", &Multiplicity, "Multiplicity/I");
   
///////////////////////////////////////////////////////////////////////////////   
// --- Create the Reader object to correct energy sum of 
// predicted Compton events by the model using the produced weights 
// from the energy regression

   TMVA::Reader *reader21 = new TMVA::Reader( "!Color:!Silent" );


   TMVA::Reader *reader31 = new TMVA::Reader( "!Color:!Silent" );

   TMVA::Reader *reader41 = new TMVA::Reader( "!Color:!Silent" );

   TMVA::Reader *reader51 = new TMVA::Reader( "!Color:!Silent" );
   
   
   reader21->AddVariable( "Energy_Scat", &Energy_scat);
   
   reader21->AddVariable( "EnDiff", &EnDiff ); 

   
   reader31->AddVariable( "Energy_Scat", &Energy_scat);
  
   reader31->AddVariable( "EnDiff", &EnDiff );

   
   reader41->AddVariable( "Energy_Scat", &Energy_scat);
  
   reader41->AddVariable( "EnDiff", &EnDiff );
   
  
   reader51->AddVariable( "Energy_Scat", &Energy_scat );
  
   reader51->AddVariable( "EnDiff", &EnDiff );

   // load the weight files for the readers
   TString method =  "BDT method";
   
   reader21->BookMVA( method, "dataset2res/weights/TMVARegression_BDTG.weights.xml" );


   reader31->BookMVA( method, "dataset3res/weights/TMVARegression_BDTG.weights.xml" );

   
   reader41->BookMVA( method, "dataset4res/weights/TMVARegression_BDTG.weights.xml" );

   
   reader51->BookMVA( method, "dataset5res/weights/TMVARegression_BDTG.weights.xml" );



      // load the input file
   TFile *input(0);
   TString fname = "./tmva__applied-s2s3s4s5ST-9var.root";
   input = TFile::Open( fname );
   
   
   TTree* theTree = NULL;
/*
   for( int i=0; i<cut.size(); i++) {
        
        //cout << i << " : " << cuts[i] << endl;
       cout << i << " : " << cut.at(i) << endl;
         
   }
*/
 
    // loop through signal and background trees
   for( int treeNumber = 0; treeNumber < 4; ++treeNumber ) {
       
       TStopwatch sw;
       sw.Start();
       
       if( treeNumber == 0 ){
           
      
           theTree = (TTree*)input->Get("SigBkg2");
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
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
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "classID", &classID );
           theTree->SetBranchAddress("cls21", &classifier21);
           theTree->SetBranchAddress("rarity21", &rar21);
           theTree->SetBranchAddress( "ECII", &fECII );
       

           std::cout << "--- Processing: " << theTree->GetEntries() << " signal & background S2 events" << std::endl;    
           Int_t nEvent = theTree->GetEntries();
           //std::cout << "--- ... All events: " << nEvent << std::endl;
           for (Long64_t ievt= 0; ievt<nEvent; ievt++) {
                   
           
                if (ievt%10000 == 0){
                   
               
                    std::cout << "--- ... Processing events: " << ievt << std::endl;
                }
    
                theTree->GetEntry(ievt);
               
               
                  
                if (classifier21 > cut.at(0) ) {
                    
                //if (classifier21 > -0.0502465) {
                 
                    if(classID == 0){
                        count20++;
                        counts++;
                        ReES = (reader21->EvaluateRegression( method ))[0];
                        outputTree->Fill();
                        outputTree2->Fill();
                        
                    }else if(classID == 1){
                        count21++;
                        countb++;
                        if(fECII->Z() == 0){
/// Bad Compton events                            
                            ReES = (reader21->EvaluateRegression( method ))[0];
                            count2bb++;
                        } else if (fECII->Z() == 1){
/// Non-Compton                            
                            ReES = (reader21->EvaluateRegression( method ))[0];
                            count2non++;
                            count22++;
                        } else if (fECII->Z() == 11) {
/// Back-Scattering                            
                            ReES = (reader21->EvaluateRegression( method ))[0];
                            count2bs++;
                        }
                            
                        outputTree->Fill();
                        outputTree2->Fill();
                        
                    }
                }
            }

       
       } else if( treeNumber == 1 ){
           
           
           theTree = (TTree*)input->Get("SigBkg3");
          
           
           theTree->SetBranchAddress("EventNumber", &EventNumber);
           theTree->SetBranchAddress("eventID", &eventID);
           theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
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
           theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
           theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
           theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
           theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
           theTree->SetBranchAddress( "EnDiff", &EnDiff );
           theTree->SetBranchAddress( "EnergySum", &EnergySum );
           theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
           theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
           theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
           theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
           theTree->SetBranchAddress( "classID", &classID );
           theTree->SetBranchAddress("cls31", &classifier31);
           theTree->SetBranchAddress("rarity31", &rar31);
           theTree->SetBranchAddress( "ECII", &fECII );
          
           std::cout << "--- Processing: " << theTree->GetEntries() << " signal & background S3 events" << std::endl;
                
           Int_t nEvent = theTree->GetEntries();
           //std::cout << "--- ... All events: " << nEvent << std::endl;
           std::vector<Int_t> tmp;
           std::vector< std::pair<Int_t, Double_t>> vec;
           
           for (Long64_t ievt=0; ievt<nEvent; ievt++) {
           
                if (ievt%10000 == 0){
                   
               
                    std::cout << "--- ... Processing events: " << ievt << std::endl;
                }
    
                theTree->GetEntry(ievt);

                if (classifier31 > cut.at(1)) {
                //if (classifier31 > -0.0427619) {
                
///True Compton events 
                    if(classID == 2){
                        count30++;
                        counts++;
                        ReES = (reader31->EvaluateRegression( method ))[0];   
                        outputTree->Fill();
                        outputTree3->Fill();
                        
                    } else if(classID == 3){
                   
///Duplicate ones are removed from Bad Compton events   
                         if(fECII->Y() == fECII->Z()){
                             count33++;
                             continue;
///Bad Compton events                             
                         } else if(fECII->Z() == 0){
                             ReES = (reader31->EvaluateRegression( method ))[0];
                             outputTree->Fill();
                             outputTree3->Fill();
                             count31++;
                             countb++;
                             count3bb++;
///Non-Compton and back-scattering events                             
                         } else {
                             count32++;
                             tmp.push_back(ievt);
                             //cout << ievt << endl;
                             vec.push_back({ ievt, fECII->T()});
                             
                         }
                         
                         
     
                    }
                   
                }
           }
//Backgrounds with the same positions are compared and Higher probability of Non-Compton 
//and Back-scattering events are chosen as the real events 

           cout << " First(B-S3) : " << tmp.front() << " , " <<tmp.at(0) << " End(B-S3) : " << tmp.back()<< ", " << tmp.size() << " , Vector_size : " << vec.size() << endl;
           
           for(Long64_t i=0; i < vec.size()/2; ++i) {
               
               countp3++;
               countb++;               
               theTree->GetEntry(vec[i].first);
               double EN1 = fECII->X();
               //cout << EN1 <<endl << endl;
               
               theTree->GetEntry(vec[i + 1].first);
               double EN2 = fECII->X();
               //cout << EN2 <<endl << endl;
               
               //cout<< vec[i].first << " : " << EN1 << " ---- " << vec[i + 1].first << " : " << EN2 << endl << endl;
               
               if(EN1 == EN2) {
                   
                                         
                  if(vec[i].second > vec[i + 1].second) {
                                  
                     countp31++;
                       
                     theTree->GetEntry(vec[i].first);
                       
                     //cout << "P31 : " << vec[i].first << " , Prob31 : " << vec[i].second << endl << endl;
                     //cout << "Lower31 : " << vec[i + 1].first << " , LowerProb31 : " << vec[i + 1].second << endl << endl;
                     
                     if(fECII->Z() == 1 || fECII->Z() == 2 ){ 
                         
                         ReES = (reader31->EvaluateRegression( method ))[0];
                         count3non++;
                     }else if(fECII->Z() == 11 || fECII->Z() == 12){
                         
                         ReES = (reader31->EvaluateRegression( method ))[0];
                         count3bs++;
                     }
                   
                     outputTree->Fill();
                     outputTree3->Fill();
                                         
                  } else if(vec[i].second < vec[i + 1].second){
                                                                       
                      countp32++;
                      theTree->GetEntry(vec[i + 1].first);
                      //cout << "P32 :" << vec[i + 1].first << " , Prob32 : " << vec[i + 1].second << endl << endl;
                      //cout << "Lower32 :" << vec[i].first << " , lowerProb32 : " << vec[i].second << endl << endl;
                      
                      if(fECII->Z() == 1 || fECII->Z() == 2 ){  
                          
                         ReES = (reader31->EvaluateRegression( method ))[0];
                         count3non++;
                      }else if(fECII->Z() == 11 || fECII->Z() == 12){
                          
                         ReES = (reader31->EvaluateRegression( method ))[0];
                         count3bs++;
                      }
                   
                      outputTree->Fill();
                      outputTree3->Fill();
                      
                  }
                  
                  i++;
                  
               } else {
                   
                   countunique3++;
                   
                   theTree->GetEntry(vec[i].first);
                   //cout << " i : " << i << " , " << "Unique :" << vec[i].first  << ", Prob : " << vec[i].second << endl << endl;
                  
                   if(fECII->Z() == 1 || fECII->Z() == 2 ){ 
                       
                       ReES = (reader31->EvaluateRegression( method ))[0];
                       count3non++;
                   }else if(fECII->Z() == 11 || fECII->Z() == 12){
                       ReES = (reader31->EvaluateRegression( method ))[0];  
                       count3bs++;
                   }
                   
                   outputTree->Fill();
                   outputTree3->Fill();
                   
              }
               
           
           }          

        } else if( treeNumber == 2 ){
               
            theTree = (TTree*)input->Get("SigBkg4");
           
            theTree->SetBranchAddress("EventNumber", &EventNumber);
            theTree->SetBranchAddress("eventID", &eventID);
            theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
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
            theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
            theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
            theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
            theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
            theTree->SetBranchAddress( "EnDiff", &EnDiff );
            theTree->SetBranchAddress( "EnergySum", &EnergySum );
            theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
            theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
            theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
            theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
            theTree->SetBranchAddress( "classID", &classID );
            theTree->SetBranchAddress("cls41", &classifier41);
            theTree->SetBranchAddress("rarity41", &rar41);
            theTree->SetBranchAddress( "ECII", &fECII );
       
            std::cout << "--- Processing: " << theTree->GetEntries() << " signal & background S4 events" << std::endl;
           
       
            
            Int_t nEvent = theTree->GetEntries();
            //std::cout << "--- ... All events: " << nEvent << std::endl;
            std::vector<Int_t> tmp;
            std::vector< std::pair<Int_t, Double_t>> vec;
            for (Long64_t ievt= 0; ievt<nEvent; ievt++) {
           
                if (ievt%10000 == 0){
                   
               
                    std::cout << "--- ... Processing events: " << ievt << std::endl;
                }
    
                theTree->GetEntry(ievt);
            
                
                if (classifier41 > cut.at(2)) {
                //if (classifier41 > -0.0626804) {
/////////////////////////////////////////////////////////////////////////////////////
 
                    if(classID == 4){
                        count40++;
                        counts++;   
                        ReES = (reader41->EvaluateRegression( method ))[0];
                        outputTree->Fill();
                        outputTree4->Fill();
                        
                    }else if(classID == 5){
                        
                         
///Duplicate ones are removed from Bad Compton events   
                       if(fECII->Y() == fECII->Z()){
                           count43++; 
                           continue;
///Bad Compton events                             
                       } else if(fECII->Z() == 0){
                           ReES = (reader41->EvaluateRegression( method ))[0];  
                           outputTree->Fill();
                           outputTree4->Fill();
                           count41++;
                           countb++;
                           count4bb++;
///Non-Compton and back-scattering events  
                                                      
                       } else {
                           
                           count42++;
                           tmp.push_back(ievt);
                             //cout << ievt << endl;
                           vec.push_back({ ievt, fECII->T()});
                           
                          
                       }
                         
                    }
                }
            }
//Backgrounds with the same positions are compared and Higher probability 
//of Non-Compton and Back-scattering events are chosen as the real events
           cout << " First(B-S4) : " << tmp.front() << " , " <<tmp.at(0) << " End(B-S4) : " << tmp.back()<< ", " << tmp.size() << " , Vector_size : "<< vec.size() << endl;
           for(Long64_t i=0; i < vec.size()/3; ++i) { 
              
               countp4++;
               countb++;               
               theTree->GetEntry(vec[i].first);
               double EN1 = fECII->X();
               
               
               theTree->GetEntry(vec[i + 1].first);
               double EN2 = fECII->X();
               
               
               //cout<< vec[i].first << " : " << EN1 << " ---- " << vec[i + 1].first << " : " << EN2 << endl << endl;
               
               if(EN1 == EN2) {
                   
                   if(vec[i].second > vec[i + 1].second) {
                       
                       countp41++;
                       
                       
                       theTree->GetEntry(vec[i].first);
                       
                       //cout << "P41 : " << vec[i].first << " , Prob41 : " << vec[i].second << endl << endl;
                       //cout << "Lower41 : " << vec[i + 1].first << " , LowerProb41 : " << vec[i + 1].second << endl << endl;
                       
                       if(fECII->Z() == 1 || fECII->Z() == 2 || fECII->Z() == 3){ 
                           
                         ReES = (reader41->EvaluateRegression( method ))[0];
                         count4non++;
                       }else if(fECII->Z() == 11 || fECII->Z() == 12 || fECII->Z() == 13){
                         ReES = (reader41->EvaluateRegression( method ))[0]; 
                         count4bs++;
                       }
                   
                       outputTree->Fill();
                       outputTree4->Fill();
                       
                   } else if(vec[i].second < vec[i + 1].second){
                       countp42++;
                       theTree->GetEntry(vec[i + 1].first);
                       //cout << "P42 :" << vec[i + 1].first << " , Prob42 : " << vec[i + 1].second << endl << endl;
                       //cout << "Lower42 :" << vec[i].first << " , lowerProb42 : " << vec[i].second << endl << endl;
                       
                       if(fECII->Z() == 1 || fECII->Z() == 2 || fECII->Z() == 3){ 
                           
                          ReES = (reader41->EvaluateRegression( method ))[0];
                          count4non++;
                     
                       }else if(fECII->Z() == 11 || fECII->Z() == 12 || fECII->Z() == 13){
                          ReES = (reader41->EvaluateRegression( method ))[0]; 
                          count4bs++;
                       }
                   
                       outputTree->Fill();
                       outputTree4->Fill();
                       
                  }
                  
                  i++;
                  
                  
               } else {
                   countunique4++;
                   
                   theTree->GetEntry(vec[i].first);
                   //cout << " i : " << i << " , " << "Unique :" << vec[i].first  << ", Prob : " << vec[i].second << endl << endl;
                   
                   if(fECII->Z() == 1 || fECII->Z() == 2 || fECII->Z() == 3){  
                       ReES = (reader41->EvaluateRegression( method ))[0];
                     count4non++;
                   }else if(fECII->Z() == 11 || fECII->Z() == 12 || fECII->Z() == 13){
                         ReES = (reader41->EvaluateRegression( method ))[0];
                         count4bs++;
                   }
                   
                   outputTree->Fill();
                   outputTree4->Fill();
                   
               }
           }

           
        } else if( treeNumber == 3 ){
               
            theTree = (TTree*)input->Get("SigBkg5");
           
            theTree->SetBranchAddress("EventNumber", &EventNumber);
            theTree->SetBranchAddress("eventID", &eventID);
            theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
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
            theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
            theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
            theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
            theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
            theTree->SetBranchAddress( "EnDiff", &EnDiff );
            theTree->SetBranchAddress( "EnergySum", &EnergySum );
            theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
            theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
            theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
            theTree->SetBranchAddress( "AngularDistribution", &AngularDistribution );
            theTree->SetBranchAddress( "classID", &classID );
            theTree->SetBranchAddress("cls51", &classifier51);
            theTree->SetBranchAddress("rarity51", &rar51);
            theTree->SetBranchAddress( "ECII", &fECII );
           
       
            std::cout << "--- Processing: " << theTree->GetEntries() << " signal & background S5 events" << std::endl;
           
       
            
            Int_t nEvent = theTree->GetEntries();
            //std::cout << "--- ... All events: " << nEvent << std::endl;
            std::vector<Int_t> tmp;
            std::vector< std::pair<Int_t, Double_t>> vec;
            
            for (Long64_t ievt= 0; ievt<nEvent; ievt++) {
           
                if (ievt%10000 == 0){
                   
               
                    std::cout << "--- ... Processing events: " << ievt << std::endl;
                }
    
                theTree->GetEntry(ievt);
                
  
                if (classifier51 > cut.at(3)) {
                //if (classifier51 > -0.152635) {
 ////////////////////////////////////////////////////
 
                    if(classID == 6){
                        count50++;
                        counts++;   
                        ReES = (reader51->EvaluateRegression( method ))[0];
                        outputTree->Fill();
                        outputTree5->Fill();
                        
                    }else if(classID == 7){
                                 
///Duplicate ones are removed from Bad Compton events   
                       if(fECII->Y() == fECII->Z()){
                           count53++; 
                           continue;
///Bad Compton events                             
                       } else if(fECII->Z() == 0){
                           ReES = (reader51->EvaluateRegression( method ))[0];  
                           outputTree->Fill();
                           outputTree5->Fill();
                           count51++;
                           countb++;
                           count5bb++;
///Non-Compton and back-scattering events  
                                                      
                       } else {
                           
                           count52++;
                           tmp.push_back(ievt);
                             //cout << ievt << endl;
                           vec.push_back({ ievt, fECII->T()});
                          
                       }
                         
                    }
                }
            }
            
//Backgrounds with the same positions are compared and Higher probability
//of Non-Compton and Back-scattering events are chosen as the real events
           cout << " First(B-S5) : " << tmp.front() << " , " <<tmp.at(0) << " End(B-S5) : " << tmp.back()<< ", " << tmp.size() << " , Vector_size : "<< vec.size() << endl;
           
           for(Long64_t i=0; i < vec.size()/4; ++i) {
          
               countp5++;
               countb++;               
               theTree->GetEntry(vec[i].first);
               double EN1 = fECII->X();
               
               
               theTree->GetEntry(vec[i + 1].first);
               double EN2 = fECII->X();
            
               //cout<< vec[i].first << " : " << EN1 << " ---- " << vec[i + 1].first << " : " << EN2 << endl << endl;
               
               if(EN1 == EN2) {
                                      
                   if(vec[i].second > vec[i + 1].second) {
                       
                       countp51++;
                       
                       
                       theTree->GetEntry(vec[i].first);
                       
                       //cout << "P51 : " << vec[i].first << " , Prob51 : " << vec[i].second << endl << endl;
                       //cout << "Lower51 : " << vec[i + 1].first << " , LowerProb51 : " << vec[i + 1].second << endl << endl;
                       
                       if(fECII->Z() == 1 || fECII->Z() == 2 || fECII->Z() == 3 || fECII->Z() == 4){ 
                           
                          ReES = (reader51->EvaluateRegression( method ))[0];
                          count5non++;
                     
                       }else if(fECII->Z() == 11 || fECII->Z() == 12 || fECII->Z() == 13 || fECII->Z() == 14){
                          ReES = (reader51->EvaluateRegression( method ))[0];
                          count5bs++;
                       }
                   
                       outputTree->Fill();
                       outputTree5->Fill();
                       
                   } else if(vec[i].second < vec[i + 1].second){
                       
                       countp52++;
                       theTree->GetEntry(vec[i + 1].first);
                       //cout << "P52 :" << vec[i + 1].first << " , Prob52 : " << vec[i + 1].second << endl << endl;
                       //cout << "Lower52 :" << vec[i].first << " , lowerProb52 : " << vec[i].second << endl << endl;
                       
                       if(fECII->Z() == 1 || fECII->Z() == 2 || fECII->Z() == 3 || fECII->Z() == 4){ 
                           
                          ReES = (reader51->EvaluateRegression( method ))[0];
                          count5non++;
                     
                       }else if(fECII->Z() == 11 || fECII->Z() == 12 || fECII->Z() == 13 || fECII->Z() == 14){
                          ReES = (reader51->EvaluateRegression( method ))[0]; 
                          count5bs++;
                       }
                   
                       outputTree->Fill();
                       outputTree5->Fill();
                       
                  }
                  
                  i++;
                
                  
               } else {
                   countunique5++;
                   
                   theTree->GetEntry(vec[i].first);
                   //cout << " i : " << i << " , " << "Unique :" << vec[i].first  << ", Prob : " << vec[i].second << endl << endl;
                 
                   if(fECII->Z() == 1 || fECII->Z() == 2 || fECII->Z() == 3 || fECII->Z() == 4){  
                       
                       ReES = (reader51->EvaluateRegression( method ))[0];
                       count5non++;
                   }else if(fECII->Z() == 11 || fECII->Z() == 12 || fECII->Z() == 13 || fECII->Z() == 14){
                       ReES = (reader51->EvaluateRegression( method ))[0]; 
                       count5bs++;
                   }
                   
                   outputTree->Fill();
                   outputTree5->Fill();
                   
               } 
           }                 
            
        }
        
       // get elapsed time
       sw.Stop();
       std::cout << "--- End of event loop: "; sw.Print();
   }

   tt->SetStats(0);
   tt->SetFillColor(kRed);
   tt->SetLineColor(kRed);
   tt->SetLineWidth(1);
    
   tt->Fill(type1[0],counts);
   tt->Fill(type1[1],countb);
   
   
   tt->LabelsDeflate("X");
//////////////////////////////////////////////
   ttt->SetStats(0);
   ttt->SetFillColor(kRed);
   ttt->SetLineColor(kRed);
   ttt->SetLineWidth(1);
    
   ttt->Fill(type2[0],count20);
   ttt->Fill(type2[4],count21);
   ttt->Fill(type2[1],count30);
   ttt->Fill(type2[5],count31+count3non+count3bs);
   ttt->Fill(type2[2],count40);
   ttt->Fill(type2[6],count41+count4non+count4bs);
   ttt->Fill(type2[3],count50);
   ttt->Fill(type2[7],count51+count5non+count5bs);
   ttt->LabelsDeflate("X");
   

   input->Close();
   // write output tree

   outputFile->Write();
   outputFile->Close();
   std::cout << "--- Created root file: \"" << outfileName.Data() << "\" containing the MVA output histograms" << std::endl;
   
   cout << " s = 2 (T) : " << count20 << " \t " << " s = 2 (B) : " << count21 <<  " \t " << " s = 2 (RB) : " << count22 << /* " \t " << " s = 2 (BS) : " << count23 << */ " \t "<< " s = 3 (T) : " << count30 << " \t " << " s = 3 (B) : " << count31 << " s = 3 (Dup) : " << count33 <<" \t " << " s = 3 (RB) : " << count32 <</*" \t " << " s = 3 (BS) : " << count34 <<*/" \t "<< " s = 4 (T) : " << count40 << " \t " << " s = 4 (B) : " << count41 << " \t " <<" s = 4 (Dup) : " << count43 << " \t " <<" s = 4 (RB) : " << count42 << /*" \t " <<" s = 4 (BS) : " << count44 <<*/" \t "<< " s = 5 (T) : " << count50 << " \t " << " s = 5 (B) : " << count51 << " \t " <<" s = 5 (Dup) : " << count53 << " \t " <<" s = 5 (RB) : " << count52 <</* " \t " <<" s = 5 (BS) : " << count54 <<*/ endl << endl;

   cout << "Event Class 3 (Non-Compton & Back-Scattering) : " << countp31 << " + " << countp32 << " + " << countunique3 << " = " << countp31+countp32+countunique3 << " , Expected Number : " << countp3 << endl << endl; 
   
   
   cout << "Event Class 4 (Non-Compton & Back-Scattering) : " << countp41 << " + " << countp42 << " + " << countunique4 << " = " << countp41+countp42+countunique4 << " , Expected Number : " << countp4 << endl << endl; 
   
   
   cout << "Event Class 5 (Non-Compton & Back-Scattering) : " << countp51 << " + " << countp52 << " + " << countunique5 << " = " << countp51+countp52+countunique5 << " , Expected Number : " << countp5 << endl << endl; 
   
   
   
   cout << "Bad Compton : " << count2bb << " , " << count3bb  << " , " << count4bb << " , " << count5bb  << endl;
   
   cout << "Non-Compton : " << count2non << " , " << count3non  << " , " << count4non << " , " << count5non  << endl;
   
   cout << "Back-Scattering Compton : " << count2bs << " , " << count3bs  << " , " << count4bs << " , " << count5bs  << endl;
   
   //std::cout << "==> Application of readers is done! combined tree created" << std::endl << std::endl;

    
    
}

// ----------------------------------------------------------------------------------------------
// Convert function get the requested input file for our reconstruction framework
// ----------------------------------------------------------------------------------------------
void Convert() {
    
   //TFile *input = TFile::Open("./tmvapp-s2s3s4s5-STRegInAna.root");
   TFile *input = TFile::Open("./tmvapp-s2s3s4s5-ST-9var.root");
   
   Float_t Pos_eX, Pos_eY, Pos_eZ, Pos_pX, Pos_pY, Pos_pZ, RealEnergy_e, RealEnergy_p; 
   Int_t EventNumber; 
   Float_t PosX_scat, PosY_scat, PosZ_scat,Energy_scat, ReEnergy_scat;
   Float_t PosX_abs, PosY_abs, PosZ_abs, Energy_abs, ReEnergy_abs;
   
   Float_t PosX_scatunc, PosY_scatunc, PosZ_scatunc,Energy_scatunc;
   Float_t PosX_absunc, PosY_absunc, PosZ_absunc, Energy_absunc;
   Float_t EnergyCluster_abs, EnergyCluster_absunc, EnergySe, EnergySeunc;
   
   Int_t Multiplicity;
  
   Float_t DiffEnergy, RatioEnergy, DiffPosition, EnergySum, PriEnergy, EnergySumunc, ReES, ReESunc, CPriEn;
   
   Int_t   classID;
   Float_t weight;
   Float_t classifier0, classifier1;
   
   TLorentzVector *fECII = new TLorentzVector();
    
   TTree* theTree = (TTree*)input->Get("Tree");
   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
  
   
   theTree->SetBranchAddress("EventNumber", &EventNumber);
   theTree->SetBranchAddress("PrimaryEnergy", &PriEnergy);
   theTree->SetBranchAddress("ReEnergySum", &ReES);
   theTree->SetBranchAddress("ECII", &fECII);
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
   theTree->SetBranchAddress( "EnergyCluster_abs", &EnergyCluster_abs );
   theTree->SetBranchAddress("EnergyCluster_absUnc", &EnergyCluster_absunc );
   theTree->SetBranchAddress( "Energy_Abs", &Energy_abs );
   theTree->SetBranchAddress( "Energy_AbsUnc", &Energy_absunc );
   theTree->SetBranchAddress( "EnergySum", &EnergySum );
   theTree->SetBranchAddress( "EnergySumUnc", &EnergySumunc );
   theTree->SetBranchAddress("classID", &classID);
   theTree->SetBranchAddress( "Multiplicity", &Multiplicity );
   theTree->SetBranchAddress( "DiffPosition", &DiffPosition );
   
   //TFile *target  = new TFile( "tmvas2s3s4s5-STRegInAna.root","RECREATE" );
   TFile *target  = new TFile( "tmvas2s3s4s5-ST-9var.root","RECREATE" );
   //Int_t EventNumber;
   Double_t PosX_1, PosY_1, PosZ_1, Energy_1;
   Double_t PosX_2, PosY_2, PosZ_2, Energy_2;
   TVector3* PosScat = new TVector3();
   TVector3* PosAbs = new TVector3();
   Double_t DEnergy, REnergy, DPosition, EnergyS, PEnergy, ReEnSunc, ReEnS, ReEnP;
   //Int_t Multi;
   
   TTree* tree = new TTree("TreeSB", "TreeSB");
   tree->Branch("EventNumber", &EventNumber);
   tree->Branch("Energy_Primary", &PEnergy);
   tree->Branch("ReEnergy_Sum", &ReEnS);
   tree->Branch("PosX_Scat", &PosX_1);
   tree->Branch("PosY_Scat", &PosY_1);
   tree->Branch("PosZ_Scat", &PosZ_1);
   tree->Branch("Pos_Scat", &PosScat);
   tree->Branch("Energy_Scat", &Energy_1);
   tree->Branch("PosX_Abs", &PosX_2);
   tree->Branch("PosY_Abs", &PosY_2);
   tree->Branch("PosZ_Abs", &PosZ_2);
   tree->Branch("Pos_Abs", &PosAbs);
   tree->Branch("Energy_Abs", &Energy_2);
   tree->Branch("Energy_Sum", &EnergyS);
   tree->Branch("Multiplicity", &Multiplicity);
   tree->Branch("classID", &classID);
   
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      theTree->GetEntry(ievt);
      EventNumber = EventNumber;
      PEnergy = PriEnergy;
      ReEnS = ReES;
      ReEnSunc = EnergySumunc;
      
      PosX_1 = PosX_scat;
      PosY_1 = PosY_scat;
      PosZ_1 = PosZ_scat;
      PosScat->SetXYZ(PosX_1,PosY_1,PosZ_1);
      Energy_1 = Energy_scat;
      
          
      PosX_2 = PosX_abs;
      PosY_2 = PosY_abs;
      PosZ_2 = PosZ_abs;
      PosAbs->SetXYZ(PosX_2,PosY_2,PosZ_2);
      Energy_2 = Energy_abs;
      
      
      EnergyS = EnergySum;
      DEnergy = DiffEnergy;
     
      Multiplicity = Multiplicity;
      classID = classID;
      
      tree->Fill();
   }
   
   tree->Write();
      
   target->Close();
   
   std::cout << "==> conversion is done!" << std::endl << std::endl;
   
}

////////////////////////////////////////////////
/////////////////////////////////////////

void Application(){
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


   cout << "========================" << endl;
   cout << "--- maximize significance" << endl;
   
   MaxiSignificanceCutApp();
   
   
   
   cout << "========================" << endl;
   cout << "--- Conversion of Final Tree" << endl;
   Convert();
   cout << "========================" << endl;

  
}

int main( int argc, char** argv ) {
   Application();
} 
