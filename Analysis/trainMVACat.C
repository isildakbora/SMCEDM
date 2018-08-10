#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainMVACat();

void trainMVACat() {
    char name[1000];
    float bkg_XSEC[2] = {1.561e+4, 2.409e+3};
    float sig_XSEC[1];
    float bkg_NORM[2] = {10000, 10000};
    float sig_NORM[1] = {10000};

    //TCut preselectionCut = "ht>400 && jetPt[5]>40 && (triggerBit[0] || triggerBit[2]) && nBJets>1 && nLeptons==0 && met<80";

    TFile *bkgSrc[2];

    bkgSrc[0] = TFile::Open("W_Jets_outfile.root");
    bkgSrc[1] = TFile::Open("DY_Jets_outfile.root");

    TFile *sigSrc = TFile::Open("Signal_outfile.root");
    sig_XSEC[0] = 1.884e+2;

    TTree *sigTree = (TTree *)sigSrc->Get("outTree");
    TTree *bkgTree[1];


    TFile *outf = new TFile("mva_SMCEDM.root", "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("factory_mva_SMCEDM", outf, "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");

    TMVA::DataLoader loader("dataset");

    loader.AddSignalTree(sigTree, sig_XSEC[0] / sig_NORM[0]);

    for (int k = 0; k < 2; k++) {
        bkgTree[k] = (TTree *)bkgSrc[k]->FindObjectAny("outTree");
        loader.AddBackgroundTree(bkgTree[k], bkg_XSEC[k] / bkg_NORM[k]);
    }

    const int NVAR = 21;
    TString VAR[NVAR] = {"nJets", "ht", "jetPt[0]", "jetPt[1]", "jetPt[2]", "jetPt[3]", "jetPt[4]", "sphericity", "aplanarity", "foxWolfram[0]", "foxWolfram[1]", "foxWolfram[2]", "foxWolfram[3]", "mTop", "yTop", "ptTop", "met", "metPhi", "nLeptons", "mW", "nBJets"};
    //char TYPE[NVAR] = {"I", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "I", "F", "I"};

    for (int i = 0; i < NVAR; i++) {
        loader.AddVariable(VAR[i]);
    }

    TCut mycuts;
    loader.PrepareTrainingAndTestTree(mycuts, "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V");

    factory->BookMethod( &loader, TMVA::Types::kKNN, "KNN",
                         "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:CreateMVAPdfs:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

    factory->BookMethod(&loader, TMVA::Types::kBDT, "BDT", "!V:NTrees=200:MinNodeSize=20.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

    factory->BookMethod( &loader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    // Use Deep Neural-Network
    
    //////General layout.
    
    TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");
    
    //////Training strategies.
    TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
    TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString trainingStrategyString ("TrainingStrategy=");
    trainingStrategyString += training0 + "|" + training1 + "|" + training2;

    //////General Options.
    TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                        "WeightInitialization=XAVIERUNIFORM");
    dnnOptions.Append (":");
    dnnOptions.Append (layoutString);
    dnnOptions.Append (":");
    dnnOptions.Append (trainingStrategyString);

    TString cpuOptions = dnnOptions + ":Architecture=CPU";

    factory->BookMethod( &loader, TMVA::Types::kDNN, "DNN CPU", cpuOptions);

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    outf->Close();
}