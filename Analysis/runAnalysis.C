
void runAnalysis() {
	gROOT->ProcessLine(".L SMCEDM_Analysis.C");
	gROOT->ProcessLine("SMCEDM_Analysis m");
	gROOT->ProcessLine("m.Loop()");
	gROOT->ProcessLine("TBrowser b");
}