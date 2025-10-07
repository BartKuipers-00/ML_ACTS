#!/usr/bin/env python3

import ROOT
import sys
import os
import math

def create_spatial_residual_plots(output_filename):
    """Create spatial residual plots per detector layer: Δz and Δφ."""
    ROOT.gROOT.SetBatch(True)  # Suppress graphics
    ROOT.gStyle.SetPalette(ROOT.kBird)

    file_path = "trackstates_ckf.root"
    
    if not os.path.exists(file_path):
        print(f"File {file_path} not found!")
        return

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print("Could not open ROOT file!")
        return
    
    track_tree = root_file.Get("trackstates")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        print("Could not find trackstates tree!")
        root_file.Close()
        return

    print(f"Found {track_tree.GetEntries()} track states entries")

    # Create canvas with 2x2 plots (RMS vs Layer + RMS vs Momentum)
    canvas = ROOT.TCanvas("c1", "Spatial Residual RMS Analysis", 1600, 1200)
    canvas.Divide(2, 2)

    # Get first entry to check data structure
    track_tree.GetEntry(0)
    volume_ids = track_tree.volume_id
    layer_ids = track_tree.layer_id
    print(f"Number of hits per track: {len(volume_ids)}")
    print(f"Volume range: {min(volume_ids)} to {max(volume_ids)}")
    print(f"Layer range: {min(layer_ids)} to {max(layer_ids)}")

    # Get actual layer range for proper binning
    min_layer = min(layer_ids)
    max_layer = max(layer_ids)

    # Plot 1: RMS of Local0 residuals vs layer (azimuthal reconstruction precision)
    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    
    h_loc0_rms_vs_layer = ROOT.TProfile("h_loc0_rms_vs_layer", 
                                       "Local0 Residual RMS vs Layer;Layer ID;RMS eLOC0 [mm]", 
                                       max_layer - min_layer + 1, min_layer - 0.5, max_layer + 0.5, "s")
    
    n_entries = track_tree.Draw("abs(res_eLOC0_ubs):layer_id>>h_loc0_rms_vs_layer", 
                               "!TMath::IsNaN(res_eLOC0_ubs)", "goff")

    
    h_loc0_rms_vs_layer.SetLineColor(2)
    h_loc0_rms_vs_layer.SetLineWidth(3)
    h_loc0_rms_vs_layer.SetMarkerColor(2)
    h_loc0_rms_vs_layer.SetMarkerStyle(21)
    h_loc0_rms_vs_layer.SetMarkerSize(1.2)
    h_loc0_rms_vs_layer.SetMinimum(0)
    h_loc0_rms_vs_layer.Draw("PE")

    # Plot 2: RMS of Local1 residuals vs layer (longitudinal reconstruction precision)
    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    
    h_loc1_rms_vs_layer = ROOT.TProfile("h_loc1_rms_vs_layer", 
                                       "Local1 Residual RMS vs Layer;Layer ID;RMS eLOC1 [mm]", 
                                       max_layer - min_layer + 1, min_layer - 0.5, max_layer + 0.5, "s")
    
    n_entries = track_tree.Draw("abs(res_eLOC1_ubs):layer_id>>h_loc1_rms_vs_layer", 
                               "!TMath::IsNaN(res_eLOC1_ubs)", "goff")
    print(f"Plot 2: Drew {n_entries} points for Local1 RMS vs Layer")
    
    h_loc1_rms_vs_layer.SetLineColor(4)
    h_loc1_rms_vs_layer.SetLineWidth(3)
    h_loc1_rms_vs_layer.SetMarkerColor(4)
    h_loc1_rms_vs_layer.SetMarkerStyle(21)
    h_loc1_rms_vs_layer.SetMarkerSize(1.2)
    h_loc1_rms_vs_layer.SetMinimum(0)
    h_loc1_rms_vs_layer.Draw("PE")

    # Plot 3: Total RMS vs momentum (Local0 - azimuthal) averaged over all layers
    canvas.cd(3)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    
    # Use TRUTH momentum to avoid reconstruction outliers extending beyond 4 GeV
    h_loc0_rms_vs_p = ROOT.TProfile("h_loc0_rms_vs_p", 
                                   "Local0 RMS vs Predicted p_{T};p_{T,pred} [GeV];Total RMS eLOC0 [mm]", 
                                   25, 0, 0, "s")  # Auto-range (state-level predicted pT)
    

    n_entries = track_tree.Draw("abs(res_eLOC0_ubs):pT_prt>>h_loc0_rms_vs_p", 
                               "!TMath::IsNaN(res_eLOC0_ubs) && pT_prt > 0", "goff")
    print(f"Plot 3: Drew {n_entries} points for Local0 RMS vs Truth Momentum")
    
    h_loc0_rms_vs_p.SetLineColor(2)
    h_loc0_rms_vs_p.SetLineWidth(3)
    h_loc0_rms_vs_p.SetMarkerColor(2)
    h_loc0_rms_vs_p.SetMarkerStyle(21)
    h_loc0_rms_vs_p.SetMarkerSize(1.2)
    h_loc0_rms_vs_p.SetMinimum(0)
    h_loc0_rms_vs_p.Draw("PE")

    # Plot 4: Total RMS vs momentum (Local1 - longitudinal) averaged over all layers
    canvas.cd(4)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    
    h_loc1_rms_vs_p = ROOT.TProfile("h_loc1_rms_vs_p", 
                                   "Local1 RMS vs Predicted p_{T};p_{T,pred} [GeV];Total RMS eLOC1 [mm]", 
                                   25, 0, 0, "s")  # Auto-range (state-level predicted pT)
    
    n_entries = track_tree.Draw("abs(res_eLOC1_ubs):pT_prt>>h_loc1_rms_vs_p", 
                               "!TMath::IsNaN(res_eLOC1_ubs) && pT_prt > 0", "goff")
    print(f"Plot 4: Drew {n_entries} points for Local1 RMS vs Truth Momentum")
    
    h_loc1_rms_vs_p.SetLineColor(4)
    h_loc1_rms_vs_p.SetLineWidth(3)
    h_loc1_rms_vs_p.SetMarkerColor(4)
    h_loc1_rms_vs_p.SetMarkerStyle(21)
    h_loc1_rms_vs_p.SetMarkerSize(1.2)
    h_loc1_rms_vs_p.SetMinimum(0)
    h_loc1_rms_vs_p.Draw("PE")

    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    
    # Print some statistics
    print(f"\nSpatial residual RMS analysis saved as {output_filename}.png")
    print("This shows the reconstruction precision analysis:")
    print("- Top plots: RMS per detector layer (which layers perform best/worst)")
    print("- Bottom plots: Total RMS vs predicted per-state p_{T} (auto-range)")
    print("- Local0: Azimuthal (φ) reconstruction precision")
    print("- Local1: Longitudinal (z) reconstruction precision")
    print("- Lower RMS = Better reconstruction precision")

    root_file.Close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_spatial_residual_analysis.py <output_filename>")
        sys.exit(1)
    
    output_filename = sys.argv[1]
    create_spatial_residual_plots(output_filename)