#define SMCEDM_Analysis_cxx
#include "SMCEDM_Analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <util.h>

void SMCEDM_Analysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L SMCEDM_Analysis.C
//      root> SMCEDM_Analysis t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  TFile* outfile_ = new TFile("Signal_outfile.root","RECREATE");

//--- book the tree -----------------------
  TTree* outTree_ = new TTree("outTree","outTree");

  //outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  int dummy_index(1), idx_1(0), idx_2(0), nJets_, nLeptons_, nBJets_;
  float ht_, met_, metPhi_, metSig_, sphericity_, foxWolfram_[4], aplanarity_, mW_, mTop_, ptTop_, yTop_;
  float mW_ref = 80;
  float jetPt_[5];

  TLorentzVector aux_jetP4(0,0,0,0);
  TLorentzVector had_W_candidateP4(0,0,0,0);
  TLorentzVector had_top_candidateP4(0,0,0,0);

  std::vector<float> delta_mW_candidate;
  std::vector<float> deltaR;
  std::vector<std::pair <int,int>> light_jet_pair_index;

  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  outTree_->Branch("nLeptons"             ,&nLeptons_          ,"nLeptons_/I");
  outTree_->Branch("nBJets"               ,&nBJets_            ,"nBJets_/I");
  outTree_->Branch("jetPt"                ,&jetPt_             ,"jetPt_[5]/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metPhi"               ,&metPhi_            ,"metPhi_/F");
  /*outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");*/
  outTree_->Branch("sphericity"           ,&sphericity_        ,"sphericity_/F");
  outTree_->Branch("aplanarity"           ,&aplanarity_        ,"aplanarity_/F");
  outTree_->Branch("foxWolfram"           ,&foxWolfram_        ,"foxWolfram_[4]/F");
  outTree_->Branch("mW"                   ,&mW_                ,"mW_/F");
  outTree_->Branch("mTop"                 ,&mTop_              ,"mTop_/F");
  outTree_->Branch("ptTop"                ,&ptTop_             ,"ptTop/F");
  outTree_->Branch("yTop"                 ,&yTop_              ,"yTop_/F");
  /*outTree_->Branch("mTTbar"               ,&mTTbar_            ,"mTTbar_/F");
  outTree_->Branch("yTTbar"               ,&yTTbar_            ,"yTTbar_/F");
  outTree_->Branch("ptTTbar"              ,&ptTTbar_           ,"ptTTbar_/F");
  outTree_->Branch("dRbbTop"              ,&dRbbTop_           ,"dRbbTop_/F");*/
  //------------------------------------------------------------------
  // Define histograms

  auto h_had_W_candidate_mass =  new TH1F("had_W_candidate_mass","had_W_candidate_mass", 200, 0, 200);

  //------------------------------------------------------------------

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      nJets_    = Jet_size;
      nLeptons_ = Electron_size + Muon_size;

      nBJets_ = 0;

      // sort jetpT in a descending order
      std::sort(std::begin(Jet_PT), std::end(Jet_PT), greater<float>());
      
      // make sure that jetPt_ array is zeroed 
      std::fill(std::begin(jetPt_), std::end(jetPt_), 0.);

      for (int i = 0; i < Jet_size; ++i){
        if(i < 5)
          jetPt_[i] = Jet_PT[i];

      	if (Jet_BTag[i] == 1)
          nBJets_++;
      }

      ht_     = ScalarHT_HT[0];
      met_    = MissingET_MET[0];
      metPhi_ = MissingET_Phi[0];

    /////////////Calculate Event Shapes/////////////
    float sumE(0.0),sumP2(0.0),sumPxx(0.0),sumPxy(0.0),sumPxz(0.0),sumPyy(0.0),sumPyz(0.0),sumPzz(0.0);

    vector<TLorentzVector> jetP4;
    vector<TLorentzVector> non_b_jetsP4;
    vector<TLorentzVector> b_jetsP4;

  	for (int i = 0; i < Jet_size; ++i){ 
      
      aux_jetP4.SetPtEtaPhiM(Jet_PT[i], Jet_Eta[i], Jet_Phi[i],Jet_Mass[i]);
      
      if(Jet_BTag[i] == 0)
        non_b_jetsP4.push_back(aux_jetP4);
      else
        b_jetsP4.push_back(aux_jetP4);

      jetP4.push_back(aux_jetP4);

  		sumE   += aux_jetP4.Energy();
  		sumP2  += aux_jetP4.P()  * aux_jetP4.P();
  		sumPxx += aux_jetP4.Px() * aux_jetP4.Px();
    	sumPxy += aux_jetP4.Px() * aux_jetP4.Py();
    	sumPxz += aux_jetP4.Px() * aux_jetP4.Pz();
    	sumPyy += aux_jetP4.Py() * aux_jetP4.Py();
    	sumPyz += aux_jetP4.Py() * aux_jetP4.Pz();
    	sumPzz += aux_jetP4.Pz() * aux_jetP4.Pz();
  	}

    if(sumP2 > 0){
      float sumPij[4] = {0.0,0.0,0.0,0.0};
      for(int k1=0;k1<jetP4.size();k1++) {
        for(int k2=k1+1;k2<jetP4.size();k2++) {
          float cosDPhi = cos(jetP4[k1].DeltaPhi(jetP4[k2]));
          float factor  = jetP4[k1].P()*jetP4[k2].P()/(sumE * sumE);
          for(int i=0;i<4;i++) {
            sumPij[i] += factor*ROOT::Math::legendre(i,cosDPhi);//legendre 0 = 1
          } 
        }
      }
      for(int i=0;i<4;i++){
        foxWolfram_[i] = sumPij[i];
      }
      //---- compute sphericity -----------------------
      float Txx = sumPxx/sumP2;
      float Tyy = sumPyy/sumP2;
      float Tzz = sumPzz/sumP2;  
      float Txy = sumPxy/sumP2; 
      float Txz = sumPxz/sumP2; 
      float Tyz = sumPyz/sumP2;
      TMatrixDSym T(3);
      T(0,0) = Txx;
      T(0,1) = Txy;
      T(0,2) = Txz;
      T(1,0) = Txy;
      T(1,1) = Tyy;
      T(1,2) = Tyz;
      T(2,0) = Txz;
      T(2,1) = Tyz;
      T(2,2) = Tzz;
      TMatrixDSymEigen TEigen(T);
      TVectorD eigenValues(TEigen.GetEigenValues());
      sphericity_ = 1.5*(eigenValues(1)+eigenValues(2));
      aplanarity_ = 1.5*eigenValues(2);
    }
    //////////////////////////

    //////Calculate W (from light quarks)//////////
    mW_=0.;

    if(non_b_jetsP4.size() > 1)
    {
      for (int i = 0; i < non_b_jetsP4.size(); ++i)
      { 
        for (int j = i+1; j < non_b_jetsP4.size(); ++j)
        { 
          // meanwhile store the indices of jet pairs
          light_jet_pair_index.push_back(std::make_pair (i,j));

          /// find the invariant mass difference between all possible non-b dijets and the mass of W (mW_ref)
          delta_mW_candidate.push_back(abs(((non_b_jetsP4[i]+non_b_jetsP4[j]).M())- mW_ref));
        }
      }

      /*for (int i = 0; i < delta_mW_candidate.size(); ++i){std::cout << mW_ref + delta_mW_candidate[i] << ",";}std::cout << std::endl;
      for (int i = 0; i < light_jet_pair_index.size(); ++i){std::cout << "(" <<light_jet_pair_index[i].first << "," <<light_jet_pair_index[i].second << ")"<< ",";}std::cout << std::endl;*/

      // 
      auto min_delta_mW = std::min_element(delta_mW_candidate.begin(),delta_mW_candidate.end());
      auto light_jet_pair_idx = std::distance(std::begin(delta_mW_candidate), min_delta_mW);
      idx_1 = light_jet_pair_index[light_jet_pair_idx].first;
      idx_2 = light_jet_pair_index[light_jet_pair_idx].second;

      // construct the hadronic W four-momenta from two light quarks
      had_W_candidateP4 = non_b_jetsP4[idx_1] + non_b_jetsP4[idx_2];
      mW_ = had_W_candidateP4.M();

      

      // Construct hadronic top quark
      // find the W-bjet system that gives the minimum deltaR
      for (int i = 0; i < b_jetsP4.size(); ++i)
      {
        deltaR.push_back(had_W_candidateP4.DeltaR(b_jetsP4[i]));
      }
      auto min_deltaR = std::min_element(deltaR.begin(),deltaR.end());
      auto idx_min_deltaR = std::distance(std::begin(deltaR), min_deltaR);
      if(nBJets_ > 0){
        had_top_candidateP4 =  had_W_candidateP4 + b_jetsP4[idx_min_deltaR];
        mTop_  = had_top_candidateP4.M();
        ptTop_ = had_top_candidateP4.Pt();
        yTop_  = had_top_candidateP4.Rapidity();
      }

      ///Clean the objects
      delta_mW_candidate.clear();
      light_jet_pair_index.clear();
      deltaR.clear();
    }


    /// Define basic selection criteria
    bool mW_cut = abs(mW_ - mW_ref) < 15.;
    bool n_light_quark = (nJets_ - nBJets_) > 1;
    
    /// Fill the tree after basic selections
    if(mW_cut*n_light_quark)
      outTree_->Fill();
   }

   outfile_->Write();
   outfile_->Close();
}