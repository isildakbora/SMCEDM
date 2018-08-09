#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainMVACat();

void trainMVACat()
{
  char name[1000];
  float XSEC[1];
  float NORM[1] = {10000};

  //TCut preselectionCut = "ht>400 && jetPt[5]>40 && (triggerBit[0] || triggerBit[2]) && nBJets>1 && nLeptons==0 && met<80";

  TFile *bkgSrc[1];

  bkgSrc[0] = TFile::Open("W_Jets_outfile.root");
  XSEC[0]   = 1.561e+4;

  TFile *sigSrc = TFile::Open("Signal_outfile.root");

  TTree *sigTree = (TTree*)sigSrc->Get("outTree"); 
  TTree *bkgTree[1];
  
  
  TFile *outf = new TFile("mva_SMCEDM.root", "RECREATE");
  TMVA::Factory* factory = new TMVA::Factory("factory_mva_SMCEDM",outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");

  TMVA::DataLoader loader("dataset");

  loader.AddSignalTree(sigTree, 1.);

  for(int k=0; k<1; k++) {
    bkgTree[k] = (TTree*)bkgSrc[k]->FindObjectAny("outTree");
    loader.AddBackgroundTree(bkgTree[k],XSEC[k]/NORM[k]);
  }
  
  const int NVAR = 21;
  TString VAR[NVAR] = {"nJets", "ht", "jetPt[0]", "jetPt[1]", "jetPt[2]", "jetPt[3]", "jetPt[4]", "sphericity", "aplanarity", "foxWolfram[0]", "foxWolfram[1]", "foxWolfram[2]", "foxWolfram[3]", "mTop", "yTop","ptTop", "met", "metPhi", "nLeptons", "mW", "nBJets"};
  //char TYPE[NVAR] = {"I", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "I", "F", "I"};

  for(int i=0; i<NVAR; i++) {
    loader.AddVariable(VAR[i]);
  }

  TCut mycuts;
  loader.PrepareTrainingAndTestTree(mycuts, "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V");

  factory->BookMethod(&loader, TMVA::Types::kBDT,"BDT", "!V:NTrees=200:MinNodeSize=20.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");


  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
  outf->Close();
}