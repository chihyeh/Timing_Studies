#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"


void TMVA_Timing()
{
    
    TFile* signalTree = new TFile("/Users/ms08962476/singularity/tev5mm_pythia6_zprime5tev_ww_with_Eta_cut_for_component_check_1_reco.root");
    TFile* backgroundTree = new TFile("/Users/ms08962476/singularity/tev5mm_pythia6_zprime5tev_qq_with_Eta_cut_for_component_check_1_reco.root");
    TTree *signal     = (TTree*)signalTree->Get("BDT_variables");
    TTree *background = (TTree*)backgroundTree->Get("BDT_variables");
    
    TString outfileName( "TMVA_for_timing_dR_PT_5TeV_reco.root" );
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
    
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

    TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
    dataloader->AddSignalTree(signal,1.);
    dataloader->AddBackgroundTree(background,1.);
  
/*
   dataloader->AddVariable("dR_Tr0T_HPt", 'F');
   dataloader->AddVariable("dR_Tr1T_HPt", 'F');
   dataloader->AddVariable("dR_Tr2T_HPt", 'F');
   dataloader->AddVariable("dR_Tr3T_HPt", 'F');
   dataloader->AddVariable("dR_Tr4T_HPt", 'F');
*/
   dataloader->AddVariable("dR_Tr0PT_HPt", 'F');
   dataloader->AddVariable("dR_Tr1PT_HPt", 'F');
   dataloader->AddVariable("dR_Tr2PT_HPt", 'F');
   dataloader->AddVariable("dR_Tr3PT_HPt", 'F');
   dataloader->AddVariable("dR_Tr4PT_HPt", 'F');

 /*
    dataloader->AddVariable("ID_Tr0T", 'F');
    dataloader->AddVariable("ID_Tr1T", 'F');
    dataloader->AddVariable("ID_Tr2T", 'F');
    dataloader->AddVariable("ID_Tr3T", 'F');
    dataloader->AddVariable("ID_Tr4T", 'F');
*/
    //====No mass soft drop====//
    //TCut mycuts= "j_tau21_b1<1 && j_tau21_b2<1 && j_c2_b1<1 && j_c2_b2<1 ";
    //TCut mycutb= "j_tau21_b1<1 && j_tau21_b2<1 && j_c2_b1<1 && j_c2_b2<1 ";
    //&& j_d2_a1_b2_mmdt<1000 && j_m2_b1_mmdt<1000 && j_m2_b2_mmdt<1000 && j_n2_b1_mmdt<1000 && j_n2_b2_mmdt<1000
    //====Have mass soft drop====//
    TCut mycuts= "dR_Tr0T_HPt<1 && dR_Tr1T_HPt<1 && dR_Tr2T_HPt<1 && dR_Tr3T_HPt<1 && dR_Tr4T_HPt<1 && dR_Tr0PT_HPt<1 && dR_Tr1PT_HPt<1 && dR_Tr2PT_HPt<1 && dR_Tr3PT_HPt<1 && dR_Tr4PT_HPt<1 && 0<dR_Tr0T_HPt && 0<dR_Tr1T_HPt && 0<dR_Tr2T_HPt && 0<dR_Tr3T_HPt && 0<dR_Tr4T_HPt && 0<dR_Tr0PT_HPt && 0<dR_Tr1PT_HPt && 0<dR_Tr2PT_HPt && 0<dR_Tr3PT_HPt && 0<dR_Tr4PT_HPt";
    TCut mycutb= "dR_Tr0T_HPt<1 && dR_Tr1T_HPt<1 && dR_Tr2T_HPt<1 && dR_Tr3T_HPt<1 && dR_Tr4T_HPt<1 && dR_Tr0PT_HPt<1 && dR_Tr1PT_HPt<1 && dR_Tr2PT_HPt<1 && dR_Tr3PT_HPt<1 && dR_Tr4PT_HPt<1 && 0<dR_Tr0T_HPt && 0<dR_Tr1T_HPt && 0<dR_Tr2T_HPt && 0<dR_Tr3T_HPt && 0<dR_Tr4T_HPt && 0<dR_Tr0PT_HPt && 0<dR_Tr1PT_HPt && 0<dR_Tr2PT_HPt && 0<dR_Tr3PT_HPt && 0<dR_Tr4PT_HPt";
    //TCut mycuts= "j_c2_b1_mmdt<100 && j_c2_b2_mmdt<100";
    //TCut mycutb= "j_c2_b1_mmdt<100 && j_c2_b2_mmdt<100";
    //=====Have mass soft drop and no mass soft drop=====//
    //TCut mycuts="j_tau21_b1<1 && j_tau21_b2<1 && j_c2_b1<1 && j_c2_b2<1 && j_d2_b1<1000 && j_d2_b2<1000 && j_d2_a1_b1<1000 && j_d2_a1_b2<1000 && j_m2_b1<1000 && j_m2_b2<1000 && j_n2_b1<1000 && j_n2_b2<1000 && j_tau21_b1_mmdt<1 && j_tau21_b2_mmdt<1 && j_c2_b1_mmdt<1 && j_c2_b2_mmdt<1 && j_d2_b1_mmdt<1000 && j_d2_b2_mmdt<1000 && j_d2_a1_b2_mmdt<1000 && j_m2_b1_mmdt<1000 && j_m2_b1_mmdt<1000 && j_m2_b2_mmdt<1000 && j_n2_b1_mmdt<1000 && j_n2_b2_mmdt<1000 && j_d2_a1_b1_mmdt<1000";
    //TCut mycutb="j_tau21_b1<1 && j_tau21_b2<1 && j_c2_b1<1 && j_c2_b2<1 && j_d2_b1<1000 && j_d2_b2<1000 && j_d2_a1_b1<1000 && j_d2_a1_b2<1000 && j_m2_b1<1000 && j_m2_b2<1000 && j_n2_b1<1000 && j_n2_b2<1000 && j_tau21_b1_mmdt<1 && j_tau21_b2_mmdt<1 && j_c2_b1_mmdt<1 && j_c2_b2_mmdt<1 && j_d2_b1_mmdt<1000 && j_d2_b2_mmdt<1000 && j_d2_a1_b2_mmdt<1000 && j_m2_b1_mmdt<1000 && j_m2_b1_mmdt<1000 && j_m2_b2_mmdt<1000 && j_n2_b1_mmdt<1000 && j_n2_b2_mmdt<1000 && j_d2_a1_b1_mmdt<1000";
    //TCut mycuts = "j_tau21_b1<1 && j_c2_b1<1 && j_mass_mmdt<800 && j_tau21_b1_mmdt<1 && j_tau21_b2_mmdt<1 && j_c2_b1_mmdt<1 && j_tau21_b2<1 && j_c2_b2<1 && j_c2_b2_mmdt<1 && jmass<500 && j_d2_b1<350 && j_d2_b2<350 "; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
    //TCut mycutb = "j_tau21_b1<1 && j_c2_b1<1 && j_mass_mmdt<800 && j_tau21_b1_mmdt<1 && j_tau21_b2_mmdt<1 && j_c2_b1_mmdt<1 && j_tau21_b2<1 && j_c2_b2<1 && j_c2_b2_mmdt<1 && jmass<500 && j_d2_b1<350 && j_d2_b2<350 "; // for example: TCut mycutb = "abs(var1)<0.5";

    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                           "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );
    //factory->BookMethod(TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
    //factory->BookMethod(TMVA::Types::kMLP,"MLP","H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator");
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
                        "!H:!V:NTrees=1042:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    
    outputFile->Write();
    outputFile->Close();
    TMVA::TMVAGui( outfileName );
    cout << "Thanks!" << endl;
}
