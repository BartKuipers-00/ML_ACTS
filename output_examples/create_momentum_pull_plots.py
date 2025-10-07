#!/usr/bin/env python3

import ROOT
import sys
import os

def create_momentum_pull_plots(output_filename):
    """Create momentum pull plot: (p_gen - p_rec) / p_gen vs pT."""
    ROOT.gROOT.SetBatch(True)  # Suppress graphics
    ROOT.gStyle.SetPalette(ROOT.kBird)

    file_path = "tracksummary_ckf.root"
    
    if not os.path.exists(file_path):
        return

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        return
    
    track_tree = root_file.Get("tracksummary")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        root_file.Close()
        return

    # Create canvas with single plot
    canvas = ROOT.TCanvas("c1", "Fractional Momentum Pull vs pT", 800, 600)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()

    # Define binning - let ROOT determine the range automatically
    pt_bins = 20
    
    # Fractional momentum pull calculation: |p_truth - p_fitted| / p_truth
    h_p_pull_profile = ROOT.TProfile("h_p_pull_profile", "Fractional Momentum Pull vs p_{T};p_{T,truth} [GeV];|p_{truth} - p_{fitted}| / p_{truth}", pt_bins, 0, 0)
    
    # Calculate momentum pull: |t_p - 1/|eQOP_fit|| / t_p
    # Note: eQOP_fit = q/p, so |p_fitted| = 1/|eQOP_fit|
    track_tree.Draw("abs(t_p - 1.0/abs(eQOP_fit))/t_p:t_pT>>h_p_pull_profile", "hasFittedParams && abs(eQOP_fit) > 1e-6", "")
    
    # Get the actual x-axis range from the data, but ensure it starts at 0
    x_min = 0.0
    x_max = h_p_pull_profile.GetXaxis().GetXmax()
    h_p_pull_profile.GetXaxis().SetRangeUser(x_min, x_max)
    
    h_p_pull_profile.SetLineColor(4)
    h_p_pull_profile.SetLineWidth(3)
    h_p_pull_profile.SetMarkerColor(4)
    h_p_pull_profile.SetMarkerStyle(20)
    h_p_pull_profile.SetMarkerSize(1.2)
    h_p_pull_profile.SetMinimum(0.0)
    h_p_pull_profile.SetMaximum(0.2)
    h_p_pull_profile.Draw("PE")
    
    # Add horizontal reference lines using actual data range, starting from x=0
    line1 = ROOT.TLine(x_min, 0, x_max, 0)
    line1.SetLineColor(1)
    line1.SetLineStyle(2)
    line1.Draw()
    
    line2 = ROOT.TLine(x_min, 0.1, x_max, 0.1)
    line2.SetLineColor(1)
    line2.SetLineStyle(3)
    line2.Draw()

    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    
    root_file.Close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(1)
    
    output_filename = sys.argv[1]
    create_momentum_pull_plots(output_filename)