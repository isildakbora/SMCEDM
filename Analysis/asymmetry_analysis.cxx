#include "SMCEDM_Analysis.h"
#include "progressbar.hpp"

void asymmetry_analysis(std::string infile_name, int n_btag_threshold, int n_light_jet_threshold, float MET_threshold) {

	//gSystem->RedirectOutput("/dev/null");
	int nthreads = 24;
	ROOT::EnableImplicitMT(nthreads);

	//auto file = new TFile("delphes.root");
	auto file = new TFile(infile_name.c_str());

	TTreeReader myReader("Delphes", file);

	// Create the TTree and branches
	TTree *outtree = new TTree("outtree", "outtree");

	int njets;
	int nleptons;
	int nbjets;
	float met;

	std::vector<float> aux_jet_pt;
	std::vector<float> aux_fox_wolfram;

	//Jet Definitions
	TTreeReaderValue<int>            rv_Jet_size(myReader, "Jet_size");
	TTreeReaderArray<float>          ra_Jet_pT(myReader,   "Jet.PT");
	TTreeReaderArray<float>          ra_Jet_Eta(myReader,  "Jet.Eta");
	TTreeReaderArray<float>          ra_Jet_Phi(myReader,  "Jet.Phi");
	TTreeReaderArray<float>          ra_Jet_Mass(myReader, "Jet.Mass");
	TTreeReaderArray<unsigned int>   ra_Jet_BTag(myReader, "Jet.BTag");

	//Electron Definitions
	TTreeReaderValue<int>            rv_Electron_size(myReader,  "Electron_size");
	TTreeReaderArray<float>          ra_Electron_pT(myReader,    "Electron.PT");
	TTreeReaderArray<float>          ra_Electron_Eta(myReader,    "Electron.Eta");
	TTreeReaderArray<float>          ra_Electron_Phi(myReader,    "Electron.Phi");
	TTreeReaderArray<float>          ra_Electron_Energy(myReader,   "Electron.T");
	TTreeReaderArray<int>            ra_Electron_Charge(myReader, "Electron.Charge");

	//Muon Definitions
	TTreeReaderValue<int>            rv_Muon_size(myReader,  "Muon_size");
	TTreeReaderArray<float>          ra_Muon_pT(myReader,    "Muon.PT");
	TTreeReaderArray<float>          ra_Muon_Eta(myReader,    "Muon.Eta");
	TTreeReaderArray<float>          ra_Muon_Phi(myReader,    "Muon.Phi");
	TTreeReaderArray<float>          ra_Muon_Energy(myReader,   "Muon.T");
	TTreeReaderArray<int>            ra_Muon_Charge(myReader, "Muon.Charge");

	//MET Definitions
	TTreeReaderArray<float>          ra_MissingET_MET(myReader,   "MissingET.MET");
	TTreeReaderArray<float>          ra_MissingET_Phi(myReader,   "MissingET.Phi");

	//HT Definitions
	TTreeReaderArray<float>          ra_scalar_HT(myReader,   "ScalarHT.HT");

	auto n_events = myReader.GetEntries(1);

	int n_btag;
	int n_btag_cut = 0;
	int n_light_jet;
	int n_events_passed = 0;

	double operator_value;
	double asymmetry_value_operator_value = 0;

	// initialize the progress bar
	const int limit = myReader.GetEntries(1);
	ProgressBar progressBar(limit, 70);


	// Loop over Events
	for (int i_event = 0; i_event < n_events; ++i_event) {
		++progressBar;
		// display the bar
		progressBar.display();

		n_btag = 0;

		bool pass_lepton_criteria(1), pass_b_jet_criteria(1), pass_MET_criteria(1), pass_njet_criteria(1);

		myReader.SetEntry(i_event);

		// put all electrons to a container (vector of leptons object)
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

		// put all muons to a container (vector of leptons object)
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

		// put all jets to a container (vector of jets object)
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
		}

		if (n_btag < n_btag_threshold) {
			pass_b_jet_criteria = 0;
		}

		if (ra_MissingET_MET[0] < MET_threshold) {
			pass_MET_criteria = 0;
		}

		// EVENT SELECTION (skips event if not matching all the criteria below)
		if (pass_lepton_criteria == 0 || pass_MET_criteria == 0 || pass_njet_criteria == 0 || pass_b_jet_criteria == 0) {
			continue;
		} else {
			n_events_passed += 1;
		}

		operator_value = calculate_operator_2(v_b_jets, the_lepton, v_light_jets.at(0));

		//std::cout << "operator_value value:" << operator_value << std::endl;
		asymmetry_value_operator_value += operator_value / abs(operator_value);
	}

	std::cout << asymmetry_value_operator_value / n_events_passed << std::endl;
}

int main(int argc, char **argv) {
	/*if (argc != 5)
	   {
	   std::cout << "Use executable with following arguments: ./SMCEDM_Analysis input event_begin event_end n_btag n_light_jet MET_threshold" << std::endl;
	   return -1;
	   }*/
	std::string infile_name   = argv[1];
	int n_btag_threshold      = atoi(argv[2]);
	int n_light_jet_threshold = atoi(argv[3]);
	float MET_threshold       = atof(argv[4]);
	asymmetry_analysis(infile_name, n_btag_threshold, n_light_jet_threshold, MET_threshold);
}
