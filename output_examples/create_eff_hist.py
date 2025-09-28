import ROOT
import sys
import os
import math

def create_combined_efficiency_and_pt_plot(output_filename):
    ROOT.gStyle.SetPalette(ROOT.kBird)

    file_path = "tracksummary_ckf.root"
    
    if not os.path.exists(file_path):
        print(f"Error: Required file '{file_path}' was not found.")
        return

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"Error: Could not open the file '{file_path}'.")
        return
    
    track_tree = root_file.Get("tracksummary")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        print("Error: Could not find TTree 'tracksummary' in the file.")
        root_file.Close()
        return


    canvas = ROOT.TCanvas("c1", "Tracking Performance Summary", 2400, 1350)
    canvas.Divide(4, 3)

    pt_bins = 100
    pt_max = 4.0
    eta_bins = 100
    eta_min = -3.0
    eta_max = 3.0

    # --- First row: plots vs. pT ---
    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    h_den_pt = ROOT.TH1F("h_den_pt", "Total Tracks vs. pT", pt_bins, 0, pt_max)
    h_num_pt = ROOT.TH1F("h_num_pt", "Good Tracks vs. pT", pt_bins, 0, pt_max)
    track_tree.Draw("t_pT>>h_den_pt", "", "goff")
    track_tree.Draw("t_pT>>h_num_pt", "trackClassification == 1", "goff")
    efficiency_plot_pt = ROOT.TEfficiency(h_num_pt, h_den_pt)
    efficiency_plot_pt.SetTitle("Tracking Efficiency vs. pT;p_{T} [GeV];Efficiency")
    efficiency_plot_pt.Draw("AP")

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.15)
    h_clone_pt = ROOT.TH1F("h_clone_pt", "Clone Tracks vs. pT", pt_bins, 0, pt_max)
    track_tree.Draw("t_pT>>h_clone_pt", "trackClassification == 2", "goff")
    clone_rate_plot_pt = ROOT.TEfficiency(h_clone_pt, h_den_pt)
    clone_rate_plot_pt.SetTitle("Clone Rate vs. pT;p_{T} [GeV];Clone Rate")
    clone_rate_plot_pt.Draw("AP")

    canvas.cd(3)
    ROOT.gPad.SetLeftMargin(0.15)
    h_fake_pt = ROOT.TH1F("h_fake_pt", "Fake Tracks vs. pT", pt_bins, 0, pt_max)
    track_tree.Draw("t_pT>>h_fake_pt", "trackClassification == 3", "goff")
    fake_rate_plot_pt = ROOT.TEfficiency(h_fake_pt, h_den_pt)
    fake_rate_plot_pt.SetTitle("Fake Rate vs. pT;p_{T} [GeV];Fake Rate")
    fake_rate_plot_pt.Draw("AP")

    canvas.cd(4)
    ROOT.gPad.SetLeftMargin(0.15)
    pt_hist = ROOT.TH1F("pt_dist", "Reconstructed Track pT Distribution;p_{T} [GeV];Number of Tracks", pt_bins, 0, pt_max)
    track_tree.Draw("t_pT>>pt_dist")
    pt_hist.SetLineColor(4)
    pt_hist.SetLineWidth(2)
    pt_hist.Draw()

    # --- Second row: plots vs. eta ---
    canvas.cd(5)
    ROOT.gPad.SetLeftMargin(0.15)
    h_den_eta = ROOT.TH1F("h_den_eta", "Total Tracks vs. eta", eta_bins, eta_min, eta_max)
    h_num_eta = ROOT.TH1F("h_num_eta", "Good Tracks vs. eta", eta_bins, eta_min, eta_max)
    track_tree.Draw("t_eta>>h_den_eta", "", "goff")
    track_tree.Draw("t_eta>>h_num_eta", "trackClassification == 1", "goff")
    efficiency_plot_eta = ROOT.TEfficiency(h_num_eta, h_den_eta)
    efficiency_plot_eta.SetTitle("Tracking Efficiency vs. eta;#eta;Efficiency")
    efficiency_plot_eta.Draw("AP")

    canvas.cd(6)
    ROOT.gPad.SetLeftMargin(0.15)
    h_clone_eta = ROOT.TH1F("h_clone_eta", "Clone Tracks vs. eta", eta_bins, eta_min, eta_max)
    track_tree.Draw("t_eta>>h_clone_eta", "trackClassification == 2", "goff")
    clone_rate_plot_eta = ROOT.TEfficiency(h_clone_eta, h_den_eta)
    clone_rate_plot_eta.SetTitle("Clone Rate vs. eta;#eta;Clone Rate")
    clone_rate_plot_eta.Draw("AP")

    canvas.cd(7)
    ROOT.gPad.SetLeftMargin(0.15)
    h_fake_eta = ROOT.TH1F("h_fake_eta", "Fake Tracks vs. eta", eta_bins, eta_min, eta_max)
    track_tree.Draw("t_eta>>h_fake_eta", "trackClassification == 3", "goff")
    fake_rate_plot_eta = ROOT.TEfficiency(h_fake_eta, h_den_eta)
    fake_rate_plot_eta.SetTitle("Fake Rate vs. eta;#eta;Fake Rate")
    fake_rate_plot_eta.Draw("AP")

    canvas.cd(8)
    ROOT.gPad.SetLeftMargin(0.15)
    eta_hist = ROOT.TH1F("eta_dist", "Reconstructed Track #eta Distribution;#eta;Number of Tracks", eta_bins, eta_min, eta_max)
    track_tree.Draw("t_eta>>eta_dist")
    eta_hist.SetLineColor(4)
    eta_hist.SetLineWidth(2)
    eta_hist.Draw()

    # --- Third row: Matching fraction plots ---
    # Pad 9: Matching Fraction vs. pT (2D histogram, log color scale)
    canvas.cd(9)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetLogz()
    h_matching_pt = ROOT.TH2F("h_matching_pt", "Matching Fraction vs. pT (log color);p_{T} [GeV];Matching Fraction", pt_bins, 0, pt_max, 100, 0, 1.05)
    track_tree.Draw("matchingFraction:t_pT>>h_matching_pt", "", "goff")
    h_matching_pt.SetTitle("Matching Fraction vs. pT; p_{T} [GeV]; Matching Fraction")
    h_matching_pt.SetMinimum(1)
    h_matching_pt.Draw("COLZ")

    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    print(f"Combined plot saved as '{output_filename}.png'")
    root_file.Close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_combined_efficiency_and_pt_plot.py <output_filename_without_extension>")
        sys.exit(1)
    
    output_filename = sys.argv[1]
    create_combined_efficiency_and_pt_plot(output_filename)