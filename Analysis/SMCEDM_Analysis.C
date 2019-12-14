#include "leptons.h"
#include "jets.h"
#include "candidate.h"


void semi_leptonic()
{
	int nthreads = 4;
	ROOT::EnableImplicitMT(nthreads);

	//auto file = new TFile("tag_1_delphes_events.root");
	auto file = new TFile("semi_lep.root");
	TTreeReader myReader("Delphes", file);

	auto outfile = new TFile("outfile_name.root", "RECREATE");

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

	// Histogram Definitions
	auto histo_bjet_pT      = new TH1F("histo_bjet_pT", "histo_bjet_pT", 200, 0, 2000);
	auto histo_light_jet_pT = new TH1F("histo_light_jet_pT", "histo_light_jet_pT", 200, 0, 2000);
	auto histo_had_top_mass = new TH1F("histo_had_top_mass", "histo_had_top_mass", 1000, 0, 1000);
	auto histo_had_top_eta  = new TH1F("histo_had_top_eta", "histo_had_top_eta", 100, -5, 5);
	auto histo_had_top_phi  = new TH1F("histo_had_top_phi", "histo_had_top_phi", 30, -3.14, 3.14);


	auto n_events = myReader.GetEntries(1);
	std::cout << "n_events:" << n_events << std::endl;

	int n_btag;
	int n_btag_cut = 0;	
	int n_light_jet;
	int n_light_jet_cut = 0;
	int n_lepton_cut = 0;
	int n_MET_cut  = 0;
	int n_events_passed = 0;
	float MET_threshold = 0.0;
	int n_positive_electrons, n_negative_electrons, n_positive_muons, n_negative_muons;

	// Loop over Events
	for (int i_event = 0; i_event < n_events; ++i_event)
	{
		n_btag = 0;
		n_positive_electrons = 0;
		n_negative_electrons = 0;
		n_positive_muons = 0;
		n_negative_muons = 0;

		bool bool_lepton_cut(1), bool_b_jet_cut(1), bool_MET_cut(1), bool_light_jet_cut(1);

		myReader.SetEntry(i_event);

		// put all electrons to a container (vector of leptons object)
		std::vector<leptons> v_electrons;
		leptons dummy_electron;
		//******* count the number of negative and positive electrons *******//
		for (int i_electron = 0; i_electron < *rv_Electron_size; ++i_electron)
		{
			dummy_electron.setPt(ra_Electron_pT.At(i_electron));
			dummy_electron.setEta(ra_Electron_Eta.At(i_electron));
			dummy_electron.setPhi(ra_Electron_Phi.At(i_electron));
			dummy_electron.setM(0.511/1000);
			dummy_electron.setQ(ra_Electron_Charge.At(i_electron));
			dummy_electron.setPDGID(11);

			v_electrons.push_back(dummy_electron);
		}
		
		// put all muons to a container (vector of leptons object)
		std::vector<leptons> v_muons;
		leptons dummy_muon;
		//******* count the number of negative and positive muons *******//
		for (int i_muon = 0; i_muon < *rv_Muon_size; ++i_muon)
		{
			dummy_muon.setPt(ra_Muon_pT.At(i_muon));
			dummy_muon.setEta(ra_Muon_Eta.At(i_muon));
			dummy_muon.setPhi(ra_Muon_Phi.At(i_muon));
			dummy_muon.setM(105/1000);
			dummy_muon.setQ(ra_Muon_Charge.At(i_muon));
			dummy_muon.setPDGID(13);

			v_muons.push_back(dummy_muon);
		}

		// put all jets to a container (vector of jets object)
		std::vector<jets> v_jets;
		jets dummy_jet;
		//******* count the number of negative and positive muons *******//
		for (int i_jet = 0; i_jet < *rv_Jet_size; ++i_jet)
		{
			dummy_jet.setPt(ra_Jet_pT.At(i_jet));
			dummy_jet.setEta(ra_Jet_Eta.At(i_jet));
			dummy_jet.setPhi(ra_Jet_Phi.At(i_jet));
			dummy_jet.setM(ra_Jet_Mass.At(i_jet));
			dummy_jet.setBTag(ra_Jet_BTag.At(i_jet));

			v_jets.push_back(dummy_jet);
		}

		if(v_electrons.size() + v_muons.size() < 1)
		{
			bool_lepton_cut = 0; 
			n_lepton_cut++;
		}
		 
		std::vector<jets> v_light_jets;

		for (int i_jet = 0; i_jet < v_jets.size(); ++i_jet)
		{
			if(v_jets.at(i_jet).getBTag()!=1)
			{
				v_light_jets.push_back(v_jets.at(i_jet));
			}
		}

		std::vector<jets> v_b_jets;

		for (int i_jet = 0; i_jet < v_jets.size(); ++i_jet)
		{
			if(v_jets.at(i_jet).getBTag()==1)
			{
				v_b_jets.push_back(v_jets.at(i_jet));
			}
		}

		n_btag = *rv_Jet_size - v_light_jets.size();

		if(v_light_jets.size() < 2)
		{
			bool_light_jet_cut = 0;
			n_light_jet_cut++;
		}

		if(n_btag < 2)
		{
			bool_b_jet_cut = 0;
			n_btag_cut++;
		}
		
		//******* skip the event if MET < threshold *******//
		//std::cout << "MET:" << ra_MissingET_MET[0] << std::endl;
		if (ra_MissingET_MET[0] < MET_threshold)
		{
			n_MET_cut++;
			bool_MET_cut = 0;
			continue;
		}

	
		if(bool_lepton_cut==0 || bool_MET_cut==0 || bool_light_jet_cut==0 || bool_b_jet_cut==0)
		{
			continue;
		}

		std::vector<candidate> v_W_candidates;
		
		for (int i_jet = 0; i_jet < v_light_jets.size(); ++i_jet)
		{
			for (int j_jet = 0; j_jet < v_light_jets.size(); ++j_jet)
			{
				if(i_jet < j_jet)
				{
					candidate dummy_candidate(v_light_jets.at(i_jet), v_light_jets.at(j_jet));
					v_W_candidates.push_back(dummy_candidate);
				}
			}
		}

		std::vector<candidate> v_had_top_candidates;

		for (int i_w_candidate = 0; i_w_candidate < v_W_candidates.size(); ++i_w_candidate)
		{
			for (int j_bjet = 0; j_bjet < v_b_jets.size(); ++j_bjet)
			{
				candidate dummy_candidate(v_W_candidates.at(i_w_candidate), v_b_jets.at(j_bjet));
				v_had_top_candidates.push_back(dummy_candidate);
			}
		}

		std::sort(v_had_top_candidates.begin(), v_had_top_candidates.end(),[](candidate first, candidate second) {return abs(first.getM()- 173.) < abs(second.getM()- 173.);});


		histo_had_top_mass->Fill(v_had_top_candidates.at(0).getM());
		histo_had_top_eta->Fill(v_had_top_candidates.at(0).get4P().Eta());
		histo_had_top_phi->Fill(v_had_top_candidates.at(0).get4P().Phi());

		n_events_passed++;
	}

	std::cout << "n_light_jet_cut: " << n_light_jet_cut << std::endl;
	std::cout << "n_btag_cut: " << n_btag_cut << std::endl;
	std::cout << "n_lepton_cut: " << n_lepton_cut << std::endl;
	std::cout << "n_MET_cut: " << n_MET_cut << std::endl;
	std::cout << "n_events_passed/n_events: " << n_events_passed << "/" << n_events << std::endl;

	TCanvas *canvas1 = new TCanvas("canvas1", "canvas1", 600, 600); 
	histo_had_top_mass->Draw("HIST");

	TCanvas *canvas2 = new TCanvas("canvas2", "canvas2", 600, 600);
	histo_had_top_eta->Draw("HIST");

	TCanvas *canvas3 = new TCanvas("canvas3", "canvas3", 600, 600);
	histo_had_top_phi->Draw("HIST");
}
