#include "SMCEDM_Analysis.h"
#include "progressbar.hpp"

void SMCEDM_Analysis(std::string infile_name, std::string outfile_name, int n_btag_threshold, int n_light_jet_threshold, float MET_threshold) {

	//gSystem->RedirectOutput("/dev/null");
	int nthreads = 24;
	ROOT::EnableImplicitMT(nthreads);

	//auto file = new TFile("delphes.root");
	auto file = new TFile(infile_name.c_str());

	TTreeReader myReader("Delphes", file);

	auto outfile = new TFile(outfile_name.c_str(), "RECREATE");

	//******* Create the TTree and branches
	TTree *outtree = new TTree("outtree", "outtree");

	int njets;
	int nleptons;
	int nbjets;

	float scalar_ht;
	float met;
	float met_phi;
	float aplanarity;
	float sphericity;

	float weight;
	float jet_pt_1, jet_pt_2, jet_pt_3, jet_pt_4;
	float fox_wolfram_1, fox_wolfram_2, fox_wolfram_3, fox_wolfram_4;
	float w_pt, w_eta, w_phi;

	float operator_4, operator_9, operator_10, operator_12, operator_14;

	std::vector<float> aux_jet_pt;
	std::vector<float> aux_fox_wolfram;

	//******* Set branches
	outtree->Branch("br_weight",        &weight,           "weight/F"       );
	outtree->Branch("br_njets",         &njets,            "njets/I"        );
	outtree->Branch("br_nleptons",      &nleptons,         "nleptons/I"     );
	outtree->Branch("br_nbjets",        &nbjets,           "nbjets/I"       );
	outtree->Branch("br_scalar_ht",     &scalar_ht,        "scalar_ht/F"    );
	outtree->Branch("br_jet_pt_1",      &jet_pt_1,         "jet_pt_1/F"     );
	outtree->Branch("br_jet_pt_2",      &jet_pt_2,         "jet_pt_2/F"     );
	outtree->Branch("br_jet_pt_3",      &jet_pt_3,         "jet_pt_3/F"     );
	outtree->Branch("br_jet_pt_4",      &jet_pt_4,         "jet_pt_4/F"     );
	outtree->Branch("br_met",           &met,              "met/F"          );
	outtree->Branch("br_met_phi",       &met_phi,          "met_phi/F"      );
	outtree->Branch("br_sphericity",    &sphericity,       "sphericity/F"   );
	outtree->Branch("br_aplanarity",    &aplanarity,       "aplanarity/F"   );
	outtree->Branch("br_fox_wolfram_1", &fox_wolfram_1,    "fox_wolfram_1/F");
	outtree->Branch("br_fox_wolfram_2", &fox_wolfram_2,    "fox_wolfram_2/F");
	outtree->Branch("br_fox_wolfram_3", &fox_wolfram_3,    "fox_wolfram_3/F");
	outtree->Branch("br_fox_wolfram_4", &fox_wolfram_4,    "fox_wolfram_4/F");
	outtree->Branch("br_w_pt",          &w_pt,             "w_pt/F"     );
	outtree->Branch("br_w_eta",         &w_eta,            "w_eta/F"     );
	outtree->Branch("br_w_phi",         &w_phi,            "w_phi/F"     );
	outtree->Branch("br_operator_4",    &operator_4,       "operator_4/F");
	outtree->Branch("br_operator_9",    &operator_9,       "operator_9/F");
	outtree->Branch("br_operator_10",   &operator_10,      "operator_10/F");
	outtree->Branch("br_operator_12",   &operator_12,      "operator_12/F");
	outtree->Branch("br_operator_14",   &operator_14,      "operator_14/F");

	//******* Book histograms
	auto h_weight        = new TH1F("weight",        "weight",        100, -1, 1);
	auto h_njets         = new TH1F("njets",         "njets",         20,  0,  20);
	auto h_nleptons      = new TH1F("nleptons",      "nleptons",      20,  0,  20);
	auto h_nbjets        = new TH1F("nbjets",        "nbjets",        20,  0,  20);
	auto h_scalar_HT     = new TH1F("scalar_ht",     "scalar_ht",     200, 0,  2000);
	auto h_jet_pt_1      = new TH1F("jet_pt_1",      "jet_pt_1",      100, 0,  1000);
	auto h_jet_pt_2      = new TH1F("jet_pt_2",      "jet_pt_2",      50,  0,  500);
	auto h_jet_pt_3      = new TH1F("jet_pt_3",      "jet_pt_3",      50,  0,  500);
	auto h_jet_pt_4      = new TH1F("jet_pt_4",      "jet_pt_4",      50,  0,  500);
	auto h_MET           = new TH1F("met",           "met",           50,  0,  500);
	auto h_MET_Phi       = new TH1F("met_phi",       "met_phi",       100, -TMath::Pi(), TMath::Pi());
	auto h_sphericity    = new TH1F("sphericity",    "sphericity",    100, 0,   1);
	auto h_aplanarity    = new TH1F("aplanarity",    "aplanarity",    100, 0,   0.5);
	auto h_Fox_Wolfram_1 = new TH1F("fox_wolfram_1", "fox_wolfram_1", 100, 0.5, 1);
	auto h_Fox_Wolfram_2 = new TH1F("fox_wolfram_2", "fox_wolfram_2", 100, 0,   1);
	auto h_Fox_Wolfram_3 = new TH1F("fox_wolfram_3", "fox_wolfram_3", 100, 0,   1);
	auto h_Fox_Wolfram_4 = new TH1F("fox_wolfram_4", "fox_wolfram_4", 100, 0,   1);

	TTreeReaderArray<float>          ra_event_weight(myReader, "Weight.Weight");

	//******* Jet Definitions
	TTreeReaderValue<int>            rv_Jet_size(myReader, "Jet_size");
	TTreeReaderArray<float>          ra_Jet_pT(myReader,   "Jet.PT");
	TTreeReaderArray<float>          ra_Jet_Eta(myReader,  "Jet.Eta");
	TTreeReaderArray<float>          ra_Jet_Phi(myReader,  "Jet.Phi");
	TTreeReaderArray<float>          ra_Jet_Mass(myReader, "Jet.Mass");
	TTreeReaderArray<unsigned int>   ra_Jet_BTag(myReader, "Jet.BTag");

	//******* Electron Definitions
	TTreeReaderValue<int>            rv_Electron_size(myReader,   "Electron_size");
	TTreeReaderArray<float>          ra_Electron_pT(myReader,     "Electron.PT");
	TTreeReaderArray<float>          ra_Electron_Eta(myReader,    "Electron.Eta");
	TTreeReaderArray<float>          ra_Electron_Phi(myReader,    "Electron.Phi");
	TTreeReaderArray<float>          ra_Electron_Energy(myReader, "Electron.T");
	TTreeReaderArray<int>            ra_Electron_Charge(myReader, "Electron.Charge");

	//******* Muon Definitions
	TTreeReaderValue<int>            rv_Muon_size(myReader,   "Muon_size");
	TTreeReaderArray<float>          ra_Muon_pT(myReader,     "Muon.PT");
	TTreeReaderArray<float>          ra_Muon_Eta(myReader,    "Muon.Eta");
	TTreeReaderArray<float>          ra_Muon_Phi(myReader,    "Muon.Phi");
	TTreeReaderArray<float>          ra_Muon_Energy(myReader, "Muon.T");
	TTreeReaderArray<int>            ra_Muon_Charge(myReader, "Muon.Charge");

	//******* MET Definitions
	TTreeReaderArray<float>          ra_MissingET_MET(myReader,   "MissingET.MET");
	TTreeReaderArray<float>          ra_MissingET_Phi(myReader,   "MissingET.Phi");

	//******* HT Definitions
	TTreeReaderArray<float>          ra_scalar_HT(myReader,   "ScalarHT.HT");

	auto n_events = myReader.GetEntries(1);

	int n_btag;
	int n_btag_cut = 0;
	int n_light_jet;
	int n_jet_cut = 0;
	int n_lepton_cut = 0;
	int n_MET_cut  = 0;
	int n_events_passed = 0;

	//******* initialize the progress bar
	const int limit = myReader.GetEntries(1);
	ProgressBar progressBar(limit, 70);


	//******* Loop over Events
	for (int i_event = 0; i_event < n_events; ++i_event) {
		
		++progressBar;
		//******* display the bar
		progressBar.display();

		n_btag = 0;

		bool pass_lepton_criteria(1), pass_b_jet_criteria(1), pass_MET_criteria(1), pass_njet_criteria(1);

		myReader.SetEntry(i_event);
		weight = ra_event_weight[0];

		//******* put all electrons to a container (vector of leptons object)
		std::vector<leptons> v_electrons;
		leptons dummy_lepton;
		leptons the_lepton;

		//******* count the number of negative and positive electrons *******//
		for (int i_electron = 0; i_electron < *rv_Electron_size; ++i_electron) {
			dummy_lepton.setPt(ra_Electron_pT.At(i_electron));
			dummy_lepton.setEta(ra_Electron_Eta.At(i_electron));
			dummy_lepton.setPhi(ra_Electron_Phi.At(i_electron));
			dummy_lepton.setM(0.511 / 1000);
			dummy_lepton.setQ(ra_Electron_Charge.At(i_electron));
			dummy_lepton.setPDGID(11);

			v_electrons.push_back(dummy_lepton);
		}

		//******* put all muons to a container (vector of leptons object)
		std::vector<leptons> v_muons;
		leptons dummy_muon;
        
		//******* count the number of negative and positive muons *******//
		for (int i_muon = 0; i_muon < *rv_Muon_size; ++i_muon) {
			dummy_lepton.setPt(ra_Muon_pT.At(i_muon));
			dummy_lepton.setEta(ra_Muon_Eta.At(i_muon));
			dummy_lepton.setPhi(ra_Muon_Phi.At(i_muon));
			dummy_lepton.setM(105 / 1000);
			dummy_lepton.setQ(ra_Muon_Charge.At(i_muon));
			dummy_lepton.setPDGID(13);
			v_muons.push_back(dummy_lepton);
		}

		//******* put all jets to a container (vector of jets object)
		std::vector<jets> v_jets;
		jets dummy_jet;
        
		//******* count the number of negative and positive muons *******//
		for (int i_jet = 0; i_jet < *rv_Jet_size; ++i_jet) {
			dummy_jet.setPt(ra_Jet_pT.At(i_jet));
			dummy_jet.setEta(ra_Jet_Eta.At(i_jet));
			dummy_jet.setPhi(ra_Jet_Phi.At(i_jet));
			dummy_jet.setM(ra_Jet_Mass.At(i_jet));
			dummy_jet.setBTag(ra_Jet_BTag.At(i_jet));

			v_jets.push_back(dummy_jet);
		}

		if (v_electrons.size() == 1 && v_muons.size() == 0) {
			pass_lepton_criteria = 1;
			the_lepton = v_electrons.at(0);
		}

		else if (v_muons.size() == 1 && v_electrons.size() == 0) {
			pass_lepton_criteria = 1;
			the_lepton = v_muons.at(0);
		}

		else {
			pass_lepton_criteria = 0;
		}

		auto v_light_jets = create_light_jet_collection(v_jets);
		auto v_b_jets     = create_b_jet_collection(v_jets);

		n_btag = v_b_jets.size();

		if (v_light_jets.size() < n_light_jet_threshold) {
			pass_njet_criteria = 0;
			n_jet_cut++;
		}

		if (n_btag < n_btag_threshold) {
			pass_b_jet_criteria = 0;
			n_btag_cut++;
		}

		if (ra_MissingET_MET[0] < MET_threshold) {
			pass_MET_criteria = 0;
			n_MET_cut++;
		}

		//******* EVENT SELECTION (skips event if not matching all the criteria below)
		if (pass_lepton_criteria == 0 || pass_MET_criteria == 0 || pass_njet_criteria == 0 || pass_b_jet_criteria == 0) {
			continue;
		}

		double W_mass = 80.379;
        
		auto v_W_candidates = create_W_candidates(v_light_jets);
		std::sort(v_W_candidates.begin(), v_W_candidates.end(),[&](candidate first, candidate second) {return abs(first.getM()- W_mass) < abs(second.getM()- W_mass);});

		if (n_btag > 1) {
			auto v_had_top_candidates = create_had_top_candidates(v_b_jets, v_W_candidates);
		}
        
		aux_jet_pt = create_jet_pt_vector(v_jets, 4);

		njets     = v_jets.size();
		nbjets    = v_b_jets.size();
		nleptons  = v_electrons.size() + v_muons.size();
		scalar_ht  = ra_scalar_HT[0];

		aplanarity      = get_event_shape_variable(v_jets, "aplanarity");
		sphericity      = get_event_shape_variable(v_jets, "sphericity");
		aux_fox_wolfram = get_fox_wolfram_parameters(v_jets, 4);

		jet_pt_1 = aux_jet_pt.at(0);
		jet_pt_2 = aux_jet_pt.at(1);
		jet_pt_3 = aux_jet_pt.at(2);
		jet_pt_4 = aux_jet_pt.at(3);

		fox_wolfram_1 = aux_fox_wolfram.at(0);
		fox_wolfram_2 = aux_fox_wolfram.at(1);
		fox_wolfram_3 = aux_fox_wolfram.at(2);
		fox_wolfram_4 = aux_fox_wolfram.at(3);

		met     = ra_MissingET_MET[0];
		met_phi = ra_MissingET_Phi[0];

		w_pt  = v_W_candidates.at(0).getPt();
		w_eta = v_W_candidates.at(0).getEta();
		w_phi = v_W_candidates.at(0).getPhi();

		operator_4  = calculate_operator_4(v_b_jets,  the_lepton, v_light_jets.at(0));
		operator_9  = calculate_operator_9(v_b_jets,  the_lepton, v_light_jets.at(0));
		operator_10 = calculate_operator_10(v_b_jets, the_lepton, v_light_jets.at(0));
		operator_12 = calculate_operator_12(v_b_jets, the_lepton, v_light_jets.at(0));
		operator_14 = calculate_operator_14(v_b_jets, the_lepton, v_light_jets.at(0));

		//std::cout << "op4:" << operator_4 << "\top9:" << operator_9 << "\top10:" << operator_10 << "\op12:" << operator_12 << "\top14:" << operator_14 << std::endl;

		h_weight->Fill(weight);
		h_njets->Fill(njets);
		h_nleptons->Fill(nleptons);
		h_nbjets->Fill(nbjets);
		h_scalar_HT->Fill(scalar_ht);
		h_jet_pt_1->Fill(jet_pt_1);
		h_jet_pt_2->Fill(jet_pt_2);
		h_jet_pt_3->Fill(jet_pt_3);
		h_jet_pt_4->Fill(jet_pt_4);
		h_MET->Fill(met);
		h_MET_Phi->Fill(met);
		h_sphericity->Fill(sphericity);
		h_aplanarity->Fill(aplanarity);
		h_Fox_Wolfram_1->Fill(fox_wolfram_1);
		h_Fox_Wolfram_2->Fill(fox_wolfram_2);
		h_Fox_Wolfram_3->Fill(fox_wolfram_3);
		h_Fox_Wolfram_4->Fill(fox_wolfram_4);

		outtree->Fill();
		n_events_passed++;
	}

	outfile->cd();
	h_weight->Write();
	h_njets->Write();
	h_nleptons->Write();
	h_nbjets->Write();
	h_scalar_HT->Write();
	h_jet_pt_1->Write();
	h_jet_pt_2->Write();
	h_jet_pt_3->Write();
	h_jet_pt_4->Write();
	h_MET->Write();
	h_MET_Phi->Write();
	h_sphericity->Write();
	h_aplanarity->Write();
	h_Fox_Wolfram_1->Write();
	h_Fox_Wolfram_2->Write();
	h_Fox_Wolfram_3->Write();
	h_Fox_Wolfram_4->Write();
	outtree->Write();
	outfile->Close();

	std::cout << "n_jet_cut: " << n_jet_cut << std::endl;
	std::cout << "n_btag_cut: " << n_btag_cut << std::endl;
	std::cout << "n_lepton_cut: " << n_lepton_cut << std::endl;
	std::cout << "n_MET_cut: " << n_MET_cut << std::endl;
	std::cout << "n_events_passed/n_events: " << n_events_passed << "/" << n_events << std::endl;
}

int main(int argc, char **argv) {
	/*if (argc != 5)
	   {
	   std::cout << "Use executable with following arguments: ./SMCEDM_Analysis input event_begin event_end n_btag n_light_jet MET_threshold" << std::endl;
	   return -1;
	   }*/
	std::string infile_name   = argv[1];
	std::string outfile_name  = argv[2];
	int n_btag_threshold      = atoi(argv[3]);
	int n_light_jet_threshold = atoi(argv[4]);
	float MET_threshold       = atof(argv[5]);

	SMCEDM_Analysis(infile_name, outfile_name, n_btag_threshold, n_light_jet_threshold, MET_threshold);
}