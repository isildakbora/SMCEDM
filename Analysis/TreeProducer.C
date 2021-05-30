#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "Math/SpecFuncMathMore.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "KKousour/TopAnalysis/plugins/TTHFlatTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

using namespace std;
using namespace reco;

TTHFlatTreeProducer::TTHFlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag>                      ("jets");
  srcMuons_           = cfg.getParameter<edm::InputTag>                      ("muons");
  srcElectrons_       = cfg.getParameter<edm::InputTag>                      ("electrons");
  srcMET_             = cfg.getParameter<edm::InputTag>                      ("met");
  srcRho_             = cfg.getParameter<edm::InputTag>                      ("rho"); 
  srcVtx_             = cfg.getParameter<edm::InputTag>                      ("vertices");
  srcQGL_             = cfg.getParameter<edm::InputTag>                      ("qgtagger"); 
  srcBtag_            = cfg.getParameter<std::string>                        ("btagger");
  kinfit_             = cfg.getParameter<std::string>                        ("kinfit");
  xmlFile_            = cfg.getParameter<std::string>                        ("xmlFile");
  nJetsMin_           = cfg.getParameter<int>                                ("nJetsMin");
  nBJetsMin_          = cfg.getParameter<int>                                ("nBJetsMin");
  etaMax_             = cfg.getParameter<double>                             ("etaMax");
  ptMin_              = cfg.getParameter<double>                             ("ptMin");
  htMin_              = cfg.getParameter<double>                             ("htMin");
  btagMinThreshold_   = cfg.getParameter<double>                             ("btagMinThreshold");
  btagMaxThreshold_   = cfg.getParameter<double>                             ("btagMaxThreshold");
  srcPU_              = cfg.getUntrackedParameter<std::string>               ("pu","");
  srcGenParticles_    = cfg.getUntrackedParameter<edm::InputTag>             ("genparticles",edm::InputTag("")); 
  triggerNames_       = cfg.getParameter<std::vector<std::string> >          ("triggerNames");
  triggerResults_     = cfg.getParameter<edm::InputTag>                      ("triggerResults");
  triggerPrescales_   = cfg.getParameter<edm::InputTag>                      ("triggerPrescales");
}
//////////////////////////////////////////////////////////////////////////////////////////
void TTHFlatTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<triggerNames_.size();i++) {
    triggerNamesHisto_->Fill(triggerNames_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetBit(TH1::kCanRebin);

  cutFlowHisto_ = fs_->make<TH1F>("CutFlow","CutFlow",1,0,1);
  cutFlowHisto_->SetBit(TH1::kCanRebin);
  
  //--- book the pileup histogram -----------
  puHisto_ = fs_->make<TH1F>("pileup","pileup",40,0,40);
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  outTree_->Branch("nLeptons"             ,&nLeptons_          ,"nLeptons_/I");
  outTree_->Branch("nBJets"               ,&nBJets_            ,"nBJets_/I");
  outTree_->Branch("status"               ,&status_            ,"status_/I");
  outTree_->Branch("prob"                 ,&prob_              ,"prob_/F");
  outTree_->Branch("chi2"                 ,&chi2_              ,"chi2_/F");
  outTree_->Branch("mva"                  ,&mva_               ,"mva_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("htBtag"               ,&htBtag_            ,"htBtag_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metPhi"               ,&metPhi_            ,"metPhi_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("sphericity"           ,&sphericity_        ,"sphericity_/F");
  outTree_->Branch("aplanarity"           ,&aplanarity_        ,"aplanarity_/F");
  outTree_->Branch("foxWolfram"           ,&foxWolfram_        ,"foxWolfram_[4]/F");

  outTree_->Branch("mbbAve"               ,&mbbAve_            ,"mbbAve_/F");
  outTree_->Branch("mbbMin"               ,&mbbMin_            ,"mbbMin_/F");
  outTree_->Branch("dRbbAve"              ,&dRbbAve_           ,"dRbbAve_/F");
  outTree_->Branch("dRbbMin"              ,&dRbbMin_           ,"dRbbMin_/F");
  outTree_->Branch("btagAve"              ,&btagAve_           ,"btagAve_/F");
  outTree_->Branch("btagMax"              ,&btagMax_           ,"btagMax_/F");
  outTree_->Branch("btagMin"              ,&btagMin_           ,"btagMin_/F");
  outTree_->Branch("qglAve"               ,&qglAve_            ,"qglAve_/F");
  outTree_->Branch("qglMin"               ,&qglMin_            ,"qglMin_/F");
  outTree_->Branch("qglMedian"            ,&qglMedian_         ,"qglMedian_/F");

  outTree_->Branch("mW"                   ,&mW_                ,"mW_[2]/F");
  outTree_->Branch("mTop"                 ,&mTop_              ,"mTop_[2]/F");
  outTree_->Branch("ptTop"                ,&ptTop_             ,"ptTop_[2]/F");
  outTree_->Branch("yTop"                 ,&yTop_              ,"yTop_[2]/F");
  outTree_->Branch("mTTbar"               ,&mTTbar_            ,"mTTbar_/F");
  outTree_->Branch("yTTbar"               ,&yTTbar_            ,"yTTbar_/F");
  outTree_->Branch("ptTTbar"              ,&ptTTbar_           ,"ptTTbar_/F");
  outTree_->Branch("dRbbTop"              ,&dRbbTop_           ,"dRbbTop_/F");
  //------------------------------------------------------------------
  isBtag_         = new std::vector<bool>;
  flavor_         = new std::vector<int>;
  pt_             = new std::vector<float>;
  btag_           = new std::vector<float>;  
  qgl_            = new std::vector<float>;
  eta_            = new std::vector<float>;
  phi_            = new std::vector<float>;
  mass_           = new std::vector<float>;
  energy_         = new std::vector<float>;
  chf_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;
  puMva_          = new std::vector<float>;
  lId_            = new std::vector<int>;
  lPt_            = new std::vector<float>;
  lEta_           = new std::vector<float>;
  lPhi_           = new std::vector<float>;
  lE_             = new std::vector<float>;
  lIso_           = new std::vector<float>;
  outTree_->Branch("jetIsBtag"            ,"vector<bool>"      ,&isBtag_);
  outTree_->Branch("jetFlavor"            ,"vector<int>"       ,&flavor_);
  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetBtag"              ,"vector<float>"     ,&btag_);  
  outTree_->Branch("jetQGL"               ,"vector<float>"     ,&qgl_);  
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);
  outTree_->Branch("jetPuMva"             ,"vector<float>"     ,&puMva_);
  outTree_->Branch("lepId"                ,"vector<int>"       ,&lId_);
  outTree_->Branch("lepPt"                ,"vector<float>"     ,&lPt_);
  outTree_->Branch("lepEta"               ,"vector<float>"     ,&lEta_);
  outTree_->Branch("lepPhi"               ,"vector<float>"     ,&lPhi_);
  outTree_->Branch("lepEnergy"            ,"vector<float>"     ,&lE_);
  outTree_->Branch("lepIso"               ,"vector<float>"     ,&lIso_);
  //------------------------------------------------------------------
  discr_ = new DiscriminatorMVA("KKousour/TopAnalysis/data/"+xmlFile_);
  triggerBit_ = new std::vector<bool>;
  triggerPre_ = new std::vector<int>;
  outTree_->Branch("triggerBit"           ,"vector<bool>"      ,&triggerBit_);
  outTree_->Branch("triggerPre"           ,"vector<int>"       ,&triggerPre_);
  //------------------- MC ---------------------------------
  outTree_->Branch("decay"         ,&decay_        ,"decay_/I");
  outTree_->Branch("HToBB"         ,&HToBB_        ,"HToBB_/O");
  outTree_->Branch("npu"           ,&npu_          ,"npu_/I");
  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void TTHFlatTreeProducer::endJob() 
{  
  delete isBtag_;
  delete flavor_;
  delete pt_;
  delete btag_;
  delete qgl_;
  delete eta_;
  delete phi_;
  delete mass_;
  delete energy_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  delete puMva_;
  delete discr_;
  delete triggerBit_;
  delete triggerPre_;
  delete lId_;
  delete lIso_;
  delete lPt_;
  delete lEta_;
  delete lPhi_;
  delete lE_;
}
//////////////////////////////////////////////////////////////////////////////////////////
bool TTHFlatTreeProducer::isGoodJet(const pat::Jet &jet)
{
  bool res = true; // by default is good, unless fails a cut bellow
  float chf = jet.chargedHadronEnergyFraction();
  float nhf = jet.neutralHadronEnergyFraction();
  float phf = jet.photonEnergyFraction();
  float muf = jet.muonEnergyFraction();
  float elf = jet.electronEnergyFraction();
  int chm   = jet.chargedHadronMultiplicity();
  int npr   = jet.neutralMultiplicity()+jet.chargedMultiplicity();
  float eta = fabs(jet.eta());
  float pt  = jet.pt();
  bool idL = (npr>1 && phf<0.99 && nhf<0.99);
  //bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
  bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || eta>2.4));
  if (!idT) res = false;
  if (pt < ptMin_) res = false;
  if (eta > etaMax_) res = false;
  return res;
}
//////////////////////////////////////////////////////////////////////////////////////////
bool TTHFlatTreeProducer::isGoodMuon(const pat::Muon &mu,const reco::Vertex &vtx,float rho)
{
  bool res = true; // by default is good, unless fails a cut bellow
  if(mu.pt() < 10) res = false;
  if(fabs(mu.eta()) > 2.4) res = false;
  if(!mu.isTightMuon(vtx)) res = false;
  // --- isolation --- those not used are commented out
  if(res && LeptonRelIso((reco::Candidate*)&mu,rho) > 0.15) res = false;
  return res;
}
//////////////////////////////////////////////////////////////////////////////////////////
float TTHFlatTreeProducer::MuonRelIso(const reco::Candidate *cand,float rho)
{
  pat::Muon mu = *((pat::Muon*)cand);
  float relIsoWithEA = 0.001;
  const int nEtaBins = 5;
  const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
  const float effectiveAreaValues[nEtaBins] = {0.0913, 0.0765, 0.0546, 0.0728, 0.1177};
  reco::MuonPFIsolation pfIso = mu.pfIsolationR03();
  // Find eta bin first. If eta>2.5, the last eta bin is used.
  int etaBin = 0;
  while(etaBin < nEtaBins-1 && abs(mu.eta()) > etaBinLimits[etaBin+1]) ++etaBin;
  float area = effectiveAreaValues[etaBin];
  relIsoWithEA = (float)(pfIso.sumChargedHadronPt+max(float(0.0),pfIso.sumNeutralHadronEt+pfIso.sumPhotonEt-rho*area))/mu.pt();
  return relIsoWithEA;
}
//////////////////////////////////////////////////////////////////////////////////////////
bool TTHFlatTreeProducer::isGoodElectron(const pat::Electron &el,const reco::Vertex &vtx,float rho)
{
  bool res = true; // by default is good, unless fails a cut bellow
  bool isEBEEGap = fabs(el.superCluster()->eta()) > 1.4442 && fabs(el.superCluster()->eta()) < 1.5660 ? 1 : 0;
  if(el.pt() < 10) res = false;
  if(fabs(el.eta()) > 2.4 && res == true) res = false;
  if(isEBEEGap && res==true) res = false;
  bool isEB = fabs(el.superCluster()->eta()) < 1.4442 ? 1 : 0;
  bool isEE = fabs(el.superCluster()->eta()) > 1.5660 ? 1 : 0;
  if(res) {
    float trackMomentumAtVtx = (float)sqrt(el.trackMomentumAtVtx().mag2());
    float ecalEnergy = (float)el.ecalEnergy();
    float full5x5_sigmaIetaIeta = (float)el.full5x5_sigmaIetaIeta();
    float dEtaIn = (float)el.deltaEtaSuperClusterTrackAtVtx();
    float dPhiIn = (float)el.deltaPhiSuperClusterTrackAtVtx();
    float HoE = (float)el.hadronicOverEm();
    float ooEmooP = (float)fabs(1/ecalEnergy - 1/trackMomentumAtVtx);
    float d0 = (float)el.gsfTrack()->dxy(vtx.position());
    float dz = (float)el.gsfTrack()->dz(vtx.position());
    int expectedMissingInnerHits = el.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    bool passConversionVeto = el.passConversionVeto();
    if(isEB) {
      if(res && full5x5_sigmaIetaIeta > 0.010557) res = false;
      if(res && fabs(dEtaIn) > 0.012442) res = false;
      if(res && fabs(dPhiIn) > 0.072624) res = false;
      if(res && HoE > 0.121476) res = false;
      if(res && ooEmooP > 0.221803) res = false;
      if(res && fabs(d0) > 0.022664) res = false;
      if(res && fabs(dz) > 0.173670) res = false;
      if(res && expectedMissingInnerHits >= 2 ) res = false;
      if(res && passConversionVeto == false ) res = false;
    }
    if(isEE) {
      if(res && full5x5_sigmaIetaIeta > 0.032602) res = false;
      if(res && fabs(dEtaIn) > 0.010654) res = false;
      if(res && fabs(dPhiIn) > 0.145129) res = false;
      if(res && HoE > 0.131862) res = false;
      if(res && ooEmooP > 0.142283) res = false;
      if(res && fabs(d0) > 0.097358) res = false;
      if(res && fabs(dz) > 0.198444) res = false;
      if(res && expectedMissingInnerHits >= 2 ) res = false;
      if(res && passConversionVeto == false ) res = false;
    }
  }
  if(res && LeptonRelIso((reco::Candidate*)&el,rho) > 0.15) res = false;
  return res;
}
//////////////////////////////////////////////////////////////////////////////////////////
float TTHFlatTreeProducer::ElectronRelIso(const reco::Candidate *cand,float rho)
{
  pat::Electron el = *((pat::Electron*)cand); 
  float relIsoWithEA = 0;
  const int nEtaBins = 5;
  const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
  const float effectiveAreaValues[nEtaBins] = {0.1013, 0.0988, 0.0572, 0.0842, 0.1530};
  reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
  float etaSC = el.superCluster()->eta();
  // Find eta bin first. If eta>2.5, the last eta bin is used.
  int etaBin = 0;
  while(etaBin < nEtaBins-1 && abs(etaSC) > etaBinLimits[etaBin+1]) ++etaBin;
  float area = effectiveAreaValues[etaBin];
  relIsoWithEA = (float)(pfIso.sumChargedHadronPt+max(float(0.0),pfIso.sumNeutralHadronEt+pfIso.sumPhotonEt-rho*area))/el.pt();
  return relIsoWithEA;
}
//////////////////////////////////////////////////////////////////////////////////////////
void TTHFlatTreeProducer::computeEventShapes(vector<const reco::Candidate *> myObj)
{
  float sumE(0.0),sumP2(0.0),sumPxx(0.0),sumPxy(0.0),sumPxz(0.0),sumPyy(0.0),sumPyz(0.0),sumPzz(0.0);
  vector<TLorentzVector> vP4;
  for(auto & obj: myObj) {
    vP4.push_back(TLorentzVector(obj->px(),obj->py(),obj->pz(),obj->energy())); 
    sumE   += obj->energy();
    sumP2  += obj->p()  * obj->p();
    sumPxx += obj->px() * obj->px();
    sumPxy += obj->px() * obj->py();
    sumPxz += obj->px() * obj->pz();
    sumPyy += obj->py() * obj->py();
    sumPyz += obj->py() * obj->pz();
    sumPzz += obj->pz() * obj->pz();
  } 
  if (sumP2 > 0) {
    //---- loop over jet pairs and compute Fox-Wolfram moments ----
    float sumPij[4] = {0.0,0.0,0.0,0.0};
    for(unsigned int k1=0;k1<vP4.size();k1++) {
      for(unsigned int k2=k1+1;k2<vP4.size();k2++) {
        float cosDPhi = cos(deltaPhi(vP4[k1].Phi(),vP4[k2].Phi()));
        float factor  = vP4[k1].P()*vP4[k2].P()/pow(sumE,2);
        for(int i=0;i<4;i++) {
          sumPij[i] += factor*ROOT::Math::legendre(i,cosDPhi);//legendre 0 = 1
        } 
      }
    }
    for(int i=0;i<4;i++) {
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
}
//////////////////////////////////////////////////////////////////////////////////////////
void TTHFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByLabel(srcJets_,jets);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(srcMuons_,muons);

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByLabel(srcElectrons_,electrons);

  edm::Handle<pat::METCollection>  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<edm::ValueMap<float> > qgtagger;
  iEvent.getByLabel(srcQGL_,qgtagger);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(srcVtx_,recVtxs);  

  //---------- mc -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > pupInfo;
  edm::Handle<GenParticleCollection> genParticles;
  
  if (!iEvent.isRealData()) { 
    iEvent.getByLabel(edm::InputTag(srcGenParticles_),genParticles);
    bool WPlusLep(false),WMinusLep(false);
    for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
      const GenParticle &p = (*genParticles)[ip];
      if (p.pdgId() == 24) {
        for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
          int daughterID = p.daughter(k)->pdgId();
          if (daughterID == -11 || daughterID == -13 || daughterID == -15) {
            WPlusLep = true;
          }
        }
      }
      if (p.pdgId() == -24) {
        for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
          int daughterID = p.daughter(k)->pdgId();
          if (daughterID == 11 || daughterID == 13 || daughterID == 15) {
            WMinusLep = true;
          }
        }
      }
      if (p.pdgId() == 25) {
        for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
          int daughterID = p.daughter(k)->pdgId();
          if (fabs(daughterID) == 5) {
            HToBB_ = true;
          }
        }
      }
    }// end of particle loop
    if (WPlusLep && WMinusLep)   decay_ = 0;
    if (WPlusLep && !WMinusLep)  decay_ = 1;
    if (!WPlusLep && WMinusLep)  decay_ = 1;
    if (!WPlusLep && !WMinusLep) decay_ = 2;
    
    iEvent.getByLabel(srcPU_,pupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PUI;
    for(PUI = pupInfo->begin(); PUI != pupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() == 0) {
        npu_ = PUI->getTrueNumInteractions();
      }
    }
    puHisto_->Fill(npu_);   
  }
  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(triggerResults_,triggerResults);  

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByLabel(triggerPrescales_,triggerPrescales);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);  
  for(unsigned int k=0;k<triggerNames_.size();k++) {
    bool bit(false);
    int pre(1);
    for(unsigned int itrig=0;itrig<triggerResults->size();itrig++) {
      string trigger_name = string(names.triggerName(itrig));
      //--- erase the last character, i.e. the version number----
      trigger_name.pop_back();
      if (trigger_name == triggerNames_[k]) {
        bit = triggerResults->accept(itrig);
        pre = triggerPrescales->getPrescaleForIndex(itrig);
        if (bit) {
          triggerPassHisto_->Fill(triggerNames_[k].c_str(),1);
        } 
      }
    }
    triggerBit_->push_back(bit);   
    triggerPre_->push_back(pre);
  }   
  vector<const reco::Candidate *> myLeptons;
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  if (cut_vtx) {
    //----- loop over leptons --------------------
    for (const pat::Muon &mu : *muons) {
      if (isGoodMuon(mu,(*recVtxs)[0],*rho)) myLeptons.push_back(&mu);
    }
    for (const pat::Electron &el : *electrons) {
      if (isGoodElectron(el,(*recVtxs)[0],*rho)) myLeptons.push_back(&el);
    }
    std::sort(myLeptons.begin(),myLeptons.end(),[](const reco::Candidate *a,const reco::Candidate *b){return a->pt() > b->pt();});
    nLeptons_ = (int)myLeptons.size();
    for(int ii = 0 ; ii < nLeptons_; ii++) {
      lPt_->push_back(myLeptons[ii]->pt());
      lEta_->push_back(myLeptons[ii]->eta());
      lPhi_->push_back(myLeptons[ii]->phi());
      lE_->push_back(myLeptons[ii]->energy());
      lId_->push_back(myLeptons[ii]->pdgId());
      lIso_->push_back(LeptonRelIso(myLeptons[ii],*rho));
    }
  }// if vtx
  //----- PF jets ------------------------------
  vector<const reco::Candidate *> myJets;
  nJets_  = 0;
  nBJets_ = 0;
  ht_     = 0.0;
  htBtag_ = 0.0;
    
  vector<float> vqgl;
  vector<int> vBIdx; 
  vector<TLorentzVector> vBP4,vBP4noTop;
  for(pat::JetCollection::const_iterator ijet =jets->begin();ijet != jets->end(); ++ijet) {  
    if (isGoodJet(*ijet)) {
      float btag= ijet->bDiscriminator(srcBtag_.c_str());
      bool isBtag = (btag >= btagMinThreshold_ && btag < btagMaxThreshold_);
      bool isLeptonMatched = false;
      float DRmax = 0.4;
      for(auto & lep: myLeptons) if( deltaR(lep->eta(),lep->phi(),ijet->eta(),ijet->phi()) < DRmax ) isLeptonMatched = true;
      if (!isLeptonMatched) {
        myJets.push_back(&*ijet);
        flavor_        ->push_back(ijet->partonFlavour());
        chf_           ->push_back(ijet->chargedHadronEnergyFraction());
        nhf_           ->push_back(ijet->neutralHadronEnergyFraction());
        phf_           ->push_back(ijet->photonEnergyFraction());
        elf_           ->push_back(ijet->electronEnergyFraction());
        muf_           ->push_back(ijet->muonEnergyFraction());
        pt_            ->push_back(ijet->pt());
        phi_           ->push_back(ijet->phi());
        eta_           ->push_back(ijet->eta());
        mass_          ->push_back(ijet->mass());
        energy_        ->push_back(ijet->energy());
        btag_          ->push_back(btag);
        edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jets, ijet-jets->begin())); 
        qgl_           ->push_back((*qgtagger)[jetRef]);
        puMva_         ->push_back(ijet->userFloat("pileupJetId:fullDiscriminant"));
        isBtag_        ->push_back(isBtag);
        ht_ += ijet->pt();
        nJets_++;
        if (isBtag) {
          vBP4.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
          vBIdx.push_back(nJets_-1);
          htBtag_ += ijet->pt();
          nBJets_++;
        }
        else {
          if ((*qgtagger)[jetRef] >= 0) {
            vqgl.push_back((*qgtagger)[jetRef]); 
          }
        } 
      }
    }
  }// jet loop
  //--- compute event-shape variables --------
  computeEventShapes(myJets); 
  //--- find the median of the qgl values
  if (vqgl.size() > 0) {
    qglMin_    = TMath::MinElement(vqgl.size(),&vqgl[0]);
    qglAve_    = TMath::Mean(vqgl.size(),&vqgl[0]);
    qglMedian_ = TMath::Median(vqgl.size(),&vqgl[0]);  
  }
    
  //---- loop over btagged jets and compute variables ----  
  int nPairs(0);
  float sumMbb(0),sumDRbb(0),sumBtag(0);
  dRbbMin_ = 100;
  btagMax_ = -100;
  btagMin_ = 100; 
  for(unsigned int k1=0;k1<vBP4.size();k1++) {
    sumBtag += (*btag_)[vBIdx[k1]];
    if ((*btag_)[vBIdx[k1]] < btagMin_) {
      btagMin_ = (*btag_)[vBIdx[k1]];
    }
    if ((*btag_)[vBIdx[k1]] > btagMax_) {
      btagMax_ = (*btag_)[vBIdx[k1]];
    } 
    for(unsigned int k2=k1+1;k2<vBP4.size();k2++) {
      float m = (vBP4[k1]+vBP4[k2]).M();
      float dR = deltaR(vBP4[k1].Eta(),vBP4[k1].Phi(),vBP4[k2].Eta(),vBP4[k2].Phi());
      sumMbb  += m;
      sumDRbb += dR;
      if (dR < dRbbMin_) {
        dRbbMin_ = dR;
        mbbMin_  = m; 
      }
      nPairs++;
    }
  }
  if (vBP4.size() > 0) {
    btagAve_ = sumBtag/vBP4.size();
  }
  if (nPairs > 0) {
    mbbAve_  = sumMbb/nPairs;
    dRbbAve_ = sumDRbb/nPairs;
  }
        
  rho_    = *rho;
  met_    = (*met)[0].et();
  metPhi_ = (*met)[0].phi();
  if ((*met)[0].sumEt() > 0) {
    metSig_ = (*met)[0].et()/(*met)[0].sumEt();
  }
  nVtx_   = recVtxs->size();
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();
  if (kinfit_ != "") {
    edm::Handle<std::vector<double> > vchi2;
    iEvent.getByLabel(edm::InputTag(kinfit_,"Chi2"),vchi2);

    edm::Handle<std::vector<double> > vprob;
    iEvent.getByLabel(edm::InputTag(kinfit_,"Prob"),vprob);

    edm::Handle<std::vector<int> > vstatus;
    iEvent.getByLabel(edm::InputTag(kinfit_,"Status"),vstatus);

    edm::Handle<std::vector<pat::Particle> > partonsB;
    iEvent.getByLabel(edm::InputTag(kinfit_,"PartonsB"),partonsB);

    edm::Handle<std::vector<pat::Particle> > partonsBbar;
    iEvent.getByLabel(edm::InputTag(kinfit_,"PartonsBBar"),partonsBbar);

    edm::Handle<std::vector<pat::Particle> > partonsQ;
    iEvent.getByLabel(edm::InputTag(kinfit_,"PartonsLightQ"),partonsQ);

    edm::Handle<std::vector<pat::Particle> > partonsQbar;
    iEvent.getByLabel(edm::InputTag(kinfit_,"PartonsLightQBar"),partonsQbar);
  
    edm::Handle<std::vector<pat::Particle> > partonsP;
    iEvent.getByLabel(edm::InputTag(kinfit_,"PartonsLightP"),partonsP);

    edm::Handle<std::vector<pat::Particle> > partonsPbar;
    iEvent.getByLabel(edm::InputTag(kinfit_,"PartonsLightPBar"),partonsPbar); 
      
    //---- KinFit information -----------------------------
    status_   = (*vstatus)[0];
    chi2_     = (*vchi2)[0];
    prob_     = (*vprob)[0];
    dRbbTop_  = deltaR((*partonsB)[0].p4(),(*partonsBbar)[0].p4());
    mW_[0]    = ((*partonsQ)[0].p4()+(*partonsQbar)[0].p4()).mass();
    mW_[1]    = ((*partonsP)[0].p4()+(*partonsPbar)[0].p4()).mass();
    mTop_[0]  = ((*partonsQ)[0].p4()+(*partonsQbar)[0].p4()+(*partonsB)[0].p4()).mass();
    mTop_[1]  = ((*partonsP)[0].p4()+(*partonsPbar)[0].p4()+(*partonsBbar)[0].p4()).mass();
    ptTop_[0] = ((*partonsQ)[0].p4()+(*partonsQbar)[0].p4()+(*partonsB)[0].p4()).pt();
    ptTop_[1] = ((*partonsP)[0].p4()+(*partonsPbar)[0].p4()+(*partonsBbar)[0].p4()).pt();
    yTop_[0]  = ((*partonsQ)[0].p4()+(*partonsQbar)[0].p4()+(*partonsB)[0].p4()).Rapidity();
    yTop_[1]  = ((*partonsP)[0].p4()+(*partonsPbar)[0].p4()+(*partonsBbar)[0].p4()).Rapidity();
    LorentzVector p4TTbar(0,0,0,0);
    p4TTbar = (*partonsQ)[0].p4()+(*partonsQbar)[0].p4()+(*partonsB)[0].p4()+(*partonsP)[0].p4()+(*partonsPbar)[0].p4()+(*partonsBbar)[0].p4();
    mTTbar_     = p4TTbar.mass();
    yTTbar_     = p4TTbar.Rapidity();
    ptTTbar_    = p4TTbar.pt();
  }// if kinfit
  
  cutFlowHisto_->Fill("All",1);
  if (nJets_ >= nJetsMin_) {
    cutFlowHisto_->Fill("nJets",1);
    if (nBJets_ >= nBJetsMin_) {
      cutFlowHisto_->Fill("nBJets",1);
      if (ht_ > htMin_) {
        cutFlowHisto_->Fill("ht",1);
        if (nJets_ > 5 && nBJets_ > 1) {
          mva_ = discr_->eval(status_,nBJets_,nJets_,ht_,(*pt_)[2],(*pt_)[3],(*pt_)[4],(*pt_)[5],mbbMin_,dRbbMin_,sphericity_,aplanarity_,foxWolfram_[0],foxWolfram_[1],foxWolfram_[2],foxWolfram_[3],mTop_[0],ptTTbar_,mTTbar_,dRbbTop_,chi2_);
        }
        outTree_->Fill();     
      }
    }  
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void TTHFlatTreeProducer::initialize()
{
  mbbAve_         = -999;
  mbbMin_         = -999;
  dRbbAve_        = -999;
  dRbbMin_        = -999;
  qglAve_         = -999;
  qglMin_         = -999;
  qglMedian_      = -999;
  btagAve_        = -999;
  btagMax_        = -999;
  btagMin_        = -999;
  status_         = -999;
  prob_           = -999;
  chi2_           = -999;
  mva_            = -999;
  dRbbTop_        = -999;
  mW_[0]          = -999;
  mW_[1]          = -999;
  mTop_[0]        = -999;
  mTop_[1]        = -999;
  ptTop_[0]       = -999;
  ptTop_[1]       = -999;
  yTop_[0]        = -999;
  yTop_[1]        = -999;
  mTTbar_         = -999;
  yTTbar_         = -999;
  ptTTbar_        = -999;
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  nJets_          = -999;
  nLeptons_       = -999;
  nBJets_         = -999;
  rho_            = -999;
  met_            = -999;
  metPhi_         = -999;
  metSig_         = -999;
  ht_             = -999;
  htBtag_         = -999;
  sphericity_     = -999;
  aplanarity_     = -999;
  for(int i=0;i<4;i++) {
    foxWolfram_[i] = -999;
  }
  flavor_         ->clear();
  pt_             ->clear();
  eta_            ->clear();
  phi_            ->clear();
  mass_           ->clear();
  energy_         ->clear();
  chf_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  btag_           ->clear();
  qgl_            ->clear();
  puMva_          ->clear();
  isBtag_         ->clear();
  triggerBit_     ->clear();
  triggerPre_     ->clear();
  lId_            ->clear();
  lIso_           ->clear();
  lPt_            ->clear();
  lEta_           ->clear();
  lPhi_           ->clear();
  lE_             ->clear();
  //----- MC -------
  HToBB_ = false;
  decay_ = -999;
  npu_ = -999;
}
//////////////////////////////////////////////////////////////////////////////////////////
TTHFlatTreeProducer::~TTHFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(TTHFlatTreeProducer);















