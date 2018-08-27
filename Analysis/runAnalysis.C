
void runAnalysis() {
	TString filename[1] =  {"Signal"};
	
	for(int i = 0; i <(int)(sizeof(filename)/sizeof(TString)); i++)
	{
		gROOT->ProcessLine(".L SMCEDM_Analysis.C");
		gROOT->ProcessLine("SMCEDM_Analysis m(\"" + filename[i] + ".root\")");
		gROOT->ProcessLine("m.Loop(\"" + filename[i]+ +"_out.root\")");
		//gROOT->ProcessLine("TBrowser b");
	}
}
