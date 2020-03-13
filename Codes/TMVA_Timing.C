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
    vector<int> Energies={5,10,20,40};
for(int energy = 0; energy<1 ; energy++)
{
            TFile* signalTree = new TFile(Form("/Users/ms08962476/singularity/TIming_Studies/tev%dmm_pythia6_zprime%dtev_ww_with_Eta_cut_for_component_check_truth.root",Energies[energy],Energies[energy]));
            TFile* backgroundTree = new TFile(Form("/Users/ms08962476/singularity/TIming_Studies/tev%dmm_pythia6_zprime%dtev_qq_with_Eta_cut_for_component_check_truth.root",Energies[energy],Energies[energy]));
            TTree *signal     = (TTree*)signalTree->Get("BDT_variables");
            TTree *background = (TTree*)backgroundTree->Get("BDT_variables");
            
            TString outfileName( Form("TMVA_for_timing_LogPT_PT_T_%dTeV.root",Energies[energy]) );
            TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
            
            TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                                       "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

            TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
            dataloader->AddSignalTree(signal,1.);
            dataloader->AddBackgroundTree(background,1.);

            dataloader->AddVariable("PT_Tr0PT_HPt", 'F');
            dataloader->AddVariable("PT_Tr1PT_HPt", 'F');
            dataloader->AddVariable("PT_Tr2PT_HPt", 'F');
            dataloader->AddVariable("PT_Tr3PT_HPt", 'F');
            dataloader->AddVariable("PT_Tr4PT_HPt", 'F');
    
            dataloader->AddVariable("PT_Tr0T_HPt", 'F');
            dataloader->AddVariable("PT_Tr1T_HPt", 'F');
            dataloader->AddVariable("PT_Tr2T_HPt", 'F');
            dataloader->AddVariable("PT_Tr3T_HPt", 'F');
            dataloader->AddVariable("PT_Tr4T_HPt", 'F');
    
       /*

         
            dataloader->AddVariable("ID_Tr0T", 'F');
            dataloader->AddVariable("ID_Tr1T", 'F');
            dataloader->AddVariable("ID_Tr2T", 'F');
            dataloader->AddVariable("ID_Tr3T", 'F');
            dataloader->AddVariable("ID_Tr4T", 'F');
        */
            //====No mass soft drop====//
            //TCut mycuts= "j_tau21_b1<5 && j_tau21_b2<5 && j_c2_b1<5 && j_c2_b2<5 ";
            //TCut mycutb= "j_tau21_b1<5 && j_tau21_b2<5 && j_c2_b1<5 && j_c2_b2<5 ";
            //&& j_d2_a1_b2_mmdt<5000 && j_m2_b1_mmdt<5000 && j_m2_b2_mmdt<5000 && j_n2_b1_mmdt<5000 && j_n2_b2_mmdt<5000
            //====Have mass soft drop====//
    
    
            //For the tracker
        /*
            TCut mycuts= "PT_Tr0PT_HPt_track<5 && PT_Tr1PT_HPt_track<5 && PT_Tr2PT_HPt_track<5 && PT_Tr3PT_HPt_track<5 && PT_Tr4PT_HPt_track<5 && 0<PT_Tr0PT_HPt_track && 0<PT_Tr1PT_HPt_track && 0<PT_Tr2PT_HPt_track && 0<PT_Tr3PT_HPt_track && 0<PT_Tr4PT_HPt_track";
            TCut mycutb= "PT_Tr0PT_HPt_track<5 && PT_Tr1PT_HPt_track<5 && PT_Tr2PT_HPt_track<5 && PT_Tr3PT_HPt_track<5 && PT_Tr4PT_HPt_track<5 &&  0<PT_Tr0PT_HPt_track && 0<PT_Tr1PT_HPt_track && 0<PT_Tr2PT_HPt_track && 0<PT_Tr3PT_HPt_track && 0<PT_Tr4PT_HPt_track";
        */
    TCut mycuts= "PT_Tr0PT_HPt<5 && PT_Tr1PT_HPt<5 && PT_Tr2PT_HPt<5 && PT_Tr3PT_HPt<5 && PT_Tr4PT_HPt<5 && -1<PT_Tr0PT_HPt && -1<PT_Tr1PT_HPt && -1<PT_Tr2PT_HPt && -1<PT_Tr3PT_HPt && -1<PT_Tr4PT_HPt && PT_Tr0T_HPt<5 && PT_Tr1T_HPt<5 && PT_Tr2T_HPt<5 && PT_Tr3T_HPt<5 && PT_Tr4T_HPt<5 && -1<PT_Tr0T_HPt && -1<PT_Tr1T_HPt && -1<PT_Tr2T_HPt && -1<PT_Tr3T_HPt && -1<PT_Tr4T_HPt";
    TCut mycutb= "PT_Tr0PT_HPt<5 && PT_Tr1PT_HPt<5 && PT_Tr2PT_HPt<5 && PT_Tr3PT_HPt<5 && PT_Tr4PT_HPt<5 && -1<PT_Tr0PT_HPt && -1<PT_Tr1PT_HPt && -1<PT_Tr2PT_HPt && -1<PT_Tr3PT_HPt && -1<PT_Tr4PT_HPt && PT_Tr0T_HPt<5 && PT_Tr1T_HPt<5 && PT_Tr2T_HPt<5 && PT_Tr3T_HPt<5 && PT_Tr4T_HPt<5 && -1<PT_Tr0T_HPt && -1<PT_Tr1T_HPt && -1<PT_Tr2T_HPt && -1<PT_Tr3T_HPt && -1<PT_Tr4T_HPt";

            //TCut mycuts= "j_c2_b1_mmdt<500 && j_c2_b2_mmdt<500";
            //TCut mycutb= "j_c2_b1_mmdt<500 && j_c2_b2_mmdt<500";
            //=====Have mass soft drop and no mass soft drop=====//
            //TCut mycuts="j_tau21_b1<5 && j_tau21_b2<5 && j_c2_b1<5 && j_c2_b2<5 && j_d2_b1<5000 && j_d2_b2<5000 && j_d2_a1_b1<5000 && j_d2_a1_b2<5000 && j_m2_b1<5000 && j_m2_b2<5000 && j_n2_b1<5000 && j_n2_b2<5000 && j_tau21_b1_mmdt<5 && j_tau21_b2_mmdt<5 && j_c2_b1_mmdt<5 && j_c2_b2_mmdt<5 && j_d2_b1_mmdt<5000 && j_d2_b2_mmdt<5000 && j_d2_a1_b2_mmdt<5000 && j_m2_b1_mmdt<5000 && j_m2_b1_mmdt<5000 && j_m2_b2_mmdt<5000 && j_n2_b1_mmdt<5000 && j_n2_b2_mmdt<5000 && j_d2_a1_b1_mmdt<5000";
            //TCut mycutb="j_tau21_b1<5 && j_tau21_b2<5 && j_c2_b1<5 && j_c2_b2<5 && j_d2_b1<5000 && j_d2_b2<5000 && j_d2_a1_b1<5000 && j_d2_a1_b2<5000 && j_m2_b1<5000 && j_m2_b2<5000 && j_n2_b1<5000 && j_n2_b2<5000 && j_tau21_b1_mmdt<5 && j_tau21_b2_mmdt<5 && j_c2_b1_mmdt<5 && j_c2_b2_mmdt<5 && j_d2_b1_mmdt<5000 && j_d2_b2_mmdt<5000 && j_d2_a1_b2_mmdt<5000 && j_m2_b1_mmdt<5000 && j_m2_b1_mmdt<5000 && j_m2_b2_mmdt<5000 && j_n2_b1_mmdt<5000 && j_n2_b2_mmdt<5000 && j_d2_a1_b1_mmdt<5000";
            //TCut mycuts = "j_tau21_b1<5 && j_c2_b1<5 && j_mass_mmdt<800 && j_tau21_b1_mmdt<5 && j_tau21_b2_mmdt<5 && j_c2_b1_mmdt<5 && j_tau21_b2<5 && j_c2_b2<5 && j_c2_b2_mmdt<5 && jmass<500 && j_d2_b1<350 && j_d2_b2<350 "; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<5";
            //TCut mycutb = "j_tau21_b1<5 && j_c2_b1<5 && j_mass_mmdt<800 && j_tau21_b1_mmdt<5 && j_tau21_b2_mmdt<5 && j_c2_b1_mmdt<5 && j_tau21_b2<5 && j_c2_b2<5 && j_c2_b2_mmdt<5 && jmass<500 && j_d2_b1<350 && j_d2_b2<350 "; // for example: TCut mycutb = "abs(var1)<0.5";

            dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                                   "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
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
}
