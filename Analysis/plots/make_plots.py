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


def branch_to_histo(input_file_name, tree_name, branch_name, bins):

    tfile = ROOT.TFile(input_file_name, "READ")
    tree = tfile.Get(tree_name)
    data = tree.AsMatrix(columns=[branch_name])
    histo = ROOT.TH1F("histo_"+branch_name, "histo_" +
                      branch_name, bins[0], bins[1], bins[2])

    for val in data:
        histo.Fill(val)

    return histo


def plot_stack_histo(paths, tree_name, variables, scales, colors, labels, lumi, rebin=1):
    """
    Makes all the plots of the variables inside the variables list
    """
    tmp_file = ROOT.TFile("tmp_file.root", "RECREATE")

    for variable in variables:
        idx = 0
        print(variable[0], idx)
        stack = ROOT.THStack("", "")
        stack.SetTitle(variable[0])
        canvas = ROOT.TCanvas(variable[0], variable[0], 600, 600)
        legend = ROOT.TLegend(0.60, 0.75, 0.91, 0.88)
        legend.SetNColumns(1)
        nbins = variable[1]
        xlow = variable[2]
        xup = variable[3]

        for path in paths:
            print(path)
            # if "signal" not in path:
            #    variable = "br_"+variable
            #print(path, variable, colors[idx], scales[idx], variable + "_" + str(idx))
            tfile = ROOT.TFile(path, "READ")
            tree = tfile.Get(tree_name)
            hist = ROOT.TH1F(variable[0] + "_" + str(idx),
                             variable[0] + "_" + str(idx), nbins, xlow, xup)
            tree.Draw(variable[0] + ">>" + variable[0] +
                      "_" + str(idx), "", "goff")
            histo = ROOT.gDirectory.Get(variable[0] + "_" + str(idx))

            tmp_file.cd()
            histo.Write()
            histo = tmp_file.FindObjectAny(variable[0] + "_" + str(idx))
            #histo.Smooth(2)

            if rebin != 1:
                histo.Rebin(rebin)
            histo.SetFillColorAlpha(colors[idx], 1)
            histo.SetLineWidth(0)

            #histo.Scale(1./histo.GetEntries())
            histo.Scale(scales[idx])

            if "signal" in path:
                histo_signal = histo
                histo_signal.SetMarkerStyle(20)
                histo_signal.SetMarkerSize(0.8)
                histo_signal.SetMarkerColorAlpha(ROOT.kBlack, 0.9)
                legend.AddEntry(histo_signal, labels[idx], "pl")

                #histo_signal.Draw("e1 same")
                signal_maximum = histo_signal.GetMaximum()
            else:
                stack.Add(histo)
                legend.AddEntry(histo, labels[idx], "f")
            idx = idx + 1
        # histo_signal.Draw("e1")
        stack.Draw("hist")
        stack.GetXaxis().SetTitle(
            variable[0].replace("br", "").replace("_", " "))
        stack.GetYaxis().SetTitle("Event Yield (" + "{:.0f}".format(lumi/100) + " fb^{-1})")
        #stack.GetYaxis().SetMoreLogLabels(1)

        histo_signal.Draw("e1 same")
        if "phi" not in variable[0]:
            canvas.SetLogy()

        stack.SetMaximum(max(stack.GetMaximum(), signal_maximum) * 1.8)

        # Add title
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.SetTextFont(42)
        #latex.DrawLatex(0.6, 0.935, "{:.1f} fb^{{-1}} (2012, 8 TeV)".format(lumi * scale))
        latex.DrawLatex(
            0.16, 0.935, "#bf{MG5#oplusPYTHIA8#oplusDelphes Simulation}")
        legend.Draw()

        canvas.Update()
        canvas.SaveAs("pdf/"+variable[0]+".pdf")
        canvas.SaveAs("png/"+variable[0]+".png")
    os.system("rm tmp_file.root")


if __name__ == "__main__":

    set_style()

    ROOT.gROOT.SetBatch(1)

    lumi = 50000.  # pb^-1
    xsecs = np.array([7.08400e+01, 15738, 3700, 46.4225])
    ngens = np.array([1.0e+06,     1.0e+07,   3.0e+07, 1.0e+07], dtype=float)

    colors = [ROOT.TColor.GetColor(155, 152, 204),
              ROOT.TColor.GetColor(209, 142, 8),
              ROOT.TColor.GetColor(111, 14, 230),
              ROOT.TColor.GetColor(222, 90, 106),
              ROOT.TColor.GetColor(250, 202, 255)]

    labels = ["t#bar{t} (Signal)", "W + Jets",
              "Drell-Yan + Jets", "Single Top"]
    scales = (lumi*xsecs)/ngens
    print(scales)
    paths = ["/mnt/harddisk4/scratch/signal_dtG1_flat.root",
             "/mnt/harddisk4/scratch/wjets_flat_10M.root",
             "/mnt/harddisk4/scratch/dyjets_flat_30M.root",
             "/mnt/harddisk4/scratch/single_top_merged_flat.root"]

    # variables = [["weight", 100, -1, 1],
    variables = [["njets", 20, 0, 10],
                 ["nbjets", 20, 0, 5],
                 ["scalar_ht", 50, 0, 2000],
                 ["jet_pt_1", 50, 0, 1000],
                 ["jet_pt_2", 25, 0, 500],
                 ["jet_pt_3", 25, 0, 400],
                 ["jet_pt_4", 25, 0, 250],
                 ["met", 50, 0, 500],
                 ["met_phi", 10, -ROOT.TMath.Pi(), ROOT.TMath.Pi()],
                 ["sphericity", 20, 0, 1],
                 ["aplanarity", 10, 0, 0.3],
                 ["fox_wolfram_1", 20, 0.5, 1],
                 ["fox_wolfram_2", 40, 0, 1],
                 ["fox_wolfram_3", 40, 0, 1],
                 ["fox_wolfram_4", 40, 0, 1]]

    tree_name = "outtree"

    rebin = 1
    plot_stack_histo(paths, tree_name, variables,
                     scales, colors, labels, lumi, rebin)
