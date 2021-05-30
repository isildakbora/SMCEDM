import ROOT
import os
import numpy as np
from array import array

def get_mathed_xsec(delphes_filename, nevents_generated, xsec):

    tfile = ROOT.TFile(delphes_filename, "READ")
    tree = tfile.FindObjectAny("Delphes")
    nevents = tree.GetEntriesFast()

    return xsec*(nevents/nevents_generated)

def set_style():

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
    ROOT.gStyle.SetCanvasDefH(600)
    ROOT.gStyle.SetCanvasDefW(600)
    ROOT.gStyle.SetCanvasDefX(0)
    ROOT.gStyle.SetCanvasDefY(0)

    ROOT.gStyle.SetPadTopMargin(0.08)
    ROOT.gStyle.SetPadBottomMargin(0.13)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetHistLineColor(1)
    ROOT.gStyle.SetHistLineStyle(0)
    ROOT.gStyle.SetHistLineWidth(1)
    ROOT.gStyle.SetEndErrorSize(2)
    ROOT.gStyle.SetMarkerStyle(20)

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetTitleFont(42)
    ROOT.gStyle.SetTitleColor(1)
    ROOT.gStyle.SetTitleTextColor(1)
    ROOT.gStyle.SetTitleFillColor(10)
    ROOT.gStyle.SetTitleFontSize(0.05)

    ROOT.gStyle.SetTitleColor(1, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleXOffset(1.00)
    ROOT.gStyle.SetTitleYOffset(1.60)

    ROOT.gStyle.SetLabelColor(1, "XYZ")
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
    ROOT.gStyle.SetLabelSize(0.04, "XYZ")

    ROOT.gStyle.SetAxisColor(1, "XYZ")
    ROOT.gStyle.SetStripDecimals(True)
    ROOT.gStyle.SetTickLength(0.03, "XYZ")
    ROOT.gStyle.SetNdivisions(510, "XYZ")
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    ROOT.gStyle.SetPaperSize(20., 20.)
    ROOT.gStyle.SetHatchesLineWidth(5)
    ROOT.gStyle.SetHatchesSpacing(0.05)

    ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "Y")

def get_branch_names(input_file_name, tree_name):
    """
    Gets names of all branches inside the tree and puts them into a list.
    Returns the list
    """

    tfile = ROOT.TFile(input_file_name, "READ")
    tree = tfile.Get(tree_name)
    branches = [name.GetName() for name in tree.GetListOfBranches()]

    return branches

def plot_stack_histo(paths, tree_name, variables, scales, colors, labels, rebin=4):
    """
    Makes all the plots of the variables inside the variables list
    """
    tmp_file = ROOT.TFile("tmp_file.root", "RECREATE")

    for variable in variables:
        idx = 0
        print(variable, idx)
        stack = ROOT.THStack("", "")
        stack.SetTitle(variable)
        canvas = ROOT.TCanvas(variable, variable, 600, 600)
        legend = ROOT.TLegend(0.60, 0.75, 0.91, 0.88)
        legend.SetNColumns(1)

        for path in paths:
            print(path)
            #if "signal" not in path:
            #    variable = "br_"+variable
            #print(path, variable, colors[idx], scales[idx], variable + "_" + str(idx))
            tfile = ROOT.TFile(path, "READ")
            #tree = tfile.Get(tree_name)
            #tree.Draw(variable + ">>" + variable + "_" + str(idx), "", "goff")
            #histo = ROOT.gDirectory.Get(variable + "_" + str(idx))

            print(variable.replace("br_",""))
            histo = tfile.Get(variable).Clone(variable + "_" + str(idx))
            tmp_file.cd()
            histo.Write()
            histo = tmp_file.FindObjectAny(variable + "_" + str(idx))

            if rebin != 1:
                histo.Rebin(rebin)
            histo.SetFillColorAlpha(colors[idx], 1)
            histo.SetLineWidth(0)

            histo.Scale(1./histo.GetEntries())

            if "signal" in path:
                histo_signal = histo
                histo_signal.SetMarkerStyle(20)
                histo_signal.SetMarkerSize(1.2)
                histo_signal.SetMarkerColorAlpha(ROOT.kBlack, 0.9)
                legend.AddEntry(histo_signal, labels[idx], "pl")

                #histo_signal.Draw("e1")
                signal_maximum = histo_signal.GetMaximum()
            else:
                stack.Add(histo)
                legend.AddEntry(histo, labels[idx], "f")
            idx = idx + 1
        stack.Draw("hist")
        histo_signal.Draw("e1 same")
        stack.GetXaxis().SetTitle(variable.replace("br", "").replace("_", " "))
        stack.GetYaxis().SetTitle("N_{Events}")
        #histo_signal.Draw("e1 same")
        canvas.SetLogy()

        stack.SetMaximum(max(stack.GetMaximum(), signal_maximum) * 2.8)

        # Add title
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.SetTextFont(42)
        lumi = 11.467
        #latex.DrawLatex(0.6, 0.935, "{:.1f} fb^{{-1}} (2012, 8 TeV)".format(lumi * scale))
        latex.DrawLatex(
            0.16, 0.935, "#bf{MG5#oplusPYTHIA8#oplusDelphes Simulation}")
        legend.Draw()

        canvas.Update()
        canvas.SaveAs("pdf/" + variable + ".pdf")
        canvas.SaveAs("png/" + variable + ".png")
        
    os.system("rm tmp_file.root")

if __name__ == "__main__":

    set_style()

    ROOT.gROOT.SetBatch(1)

    colors = [ROOT.TColor.GetColor(155, 152, 204),
              ROOT.TColor.GetColor(222, 90, 106),
              ROOT.TColor.GetColor("#00A88F"),
              ROOT.TColor.GetColor(100, 192, 232),
              ROOT.TColor.GetColor(250, 202, 255)]

    labels = ["t#bar{t} (Signal)", "W + Jets",
              "Drell-Yan + Jets", "Single Top"]
    scales = [1.0, 1.0, 1.0, 1.0]
    paths = ["/mnt/harddisk4/scratch/signal_dtG1_flat.root",
             "/mnt/harddisk4/scratch/wjets_flat_10M.root",
             "/mnt/harddisk4/scratch/dyjets_flat_30M.root",
             "/mnt/harddisk4/scratch/single_top_merged_flat.root"]

    tree_name = "outtree"
    variables = ["scalar_ht","jet_pt_1","jet_pt_2","jet_pt_3","jet_pt_4","met","met_phi","sphericity","aplanarity","fox_wolfram_1","fox_wolfram_2","fox_wolfram_3","fox_wolfram_4"]#get_branch_names(paths[0], tree_name)
    print(variables)
    rebin = 1
    plot_stack_histo(paths, tree_name, variables, scales, colors, labels, rebin)