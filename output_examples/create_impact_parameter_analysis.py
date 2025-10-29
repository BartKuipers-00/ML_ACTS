#!/usr/bin/env python3

import ROOT
import sys
import os

def create_impact_parameter_plots(output_filename):
    """Create two line plots: average |d0| vs pT and average |z0| vs pT."""
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

    # Create canvas with 2 plots side by side
    canvas = ROOT.TCanvas("c1", "Impact Parameter Resolution vs pT", 1600, 600)
    canvas.Divide(2, 1)

    # Define binning - let ROOT determine the range automatically
    pt_bins = 20
    
    # Plot 1: Average |d0| vs pT using TProfile
    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    h_d0_profile = ROOT.TProfile("h_d0_profile", "Average |d_{0}| vs p_{T};p_{T} [GeV];Average |d_{0}| [mm]", pt_bins, 0, 0)
    track_tree.Draw("abs(eLOC0_fit):t_pT>>h_d0_profile", "hasFittedParams", "")
    h_d0_profile.SetLineColor(2)
    h_d0_profile.SetLineWidth(3)
    h_d0_profile.SetMarkerColor(2)
    h_d0_profile.SetMarkerStyle(20)
    h_d0_profile.SetMarkerSize(1.2)
    h_d0_profile.SetMinimum(0)
    # Fix y-axis scale as requested
    h_d0_profile.SetMaximum(1)
    h_d0_profile.Draw("PE")
    
    # Get the actual x-axis range from the data, but ensure it starts at 0
    x_min = 0.0
    x_max = h_d0_profile.GetXaxis().GetXmax()
    h_d0_profile.GetXaxis().SetRangeUser(x_min, x_max)

    # Plot 2: Average |z0| vs pT using TProfile
    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    h_z0_profile = ROOT.TProfile("h_z0_profile", "Average |z_{0}| vs p_{T};p_{T} [GeV];Average |z_{0}| [mm]", pt_bins, 0, 0)
    track_tree.Draw("abs(eLOC1_fit):t_pT>>h_z0_profile", "hasFittedParams", "")
    h_z0_profile.SetLineColor(4)
    h_z0_profile.SetLineWidth(3)
    h_z0_profile.SetMarkerColor(4)
    h_z0_profile.SetMarkerStyle(20)
    h_z0_profile.SetMarkerSize(1.2)
    h_z0_profile.SetMinimum(0)
    # Fix y-axis scale as requested
    h_z0_profile.SetMaximum(1)
    h_z0_profile.Draw("PE")
    
    # Ensure both plots use the same x-axis range starting from 0
    h_z0_profile.GetXaxis().SetRangeUser(x_min, x_max)

    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    
    root_file.Close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(1)
    
    output_filename = sys.argv[1]
    create_impact_parameter_plots(output_filename)