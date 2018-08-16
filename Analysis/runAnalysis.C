
void runAnalysis() {
	TString filename[1] =  {"WJets"};
	
	for(int i = 0; i <(int)(sizeof(filename)/sizeof(TString)); i++)
	{
		gROOT->ProcessLine(".L SMCEDM_Analysis.C\n
							SMCEDM_Analysis m(\"" + filename[i] + ".root\")\n
							m.Loop(\"" + filename[i]+ +"_out.root\")");
	}
}
