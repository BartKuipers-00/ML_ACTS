import ROOT
import sys
import os
import math

def create_combined_efficiency_and_pt_plot(output_filename):
    ROOT.gStyle.SetPalette(ROOT.kBird)

    file_path = "tracksummary_ckf.root"
    particles_file_path = "particles.root"
    
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

    # Try to open particles.root for simulated tracks
    particles_file = None
    particles_tree = None
    if os.path.exists(particles_file_path):
        particles_file = ROOT.TFile.Open(particles_file_path)
        if particles_file and not particles_file.IsZombie():
            particles_tree = particles_file.Get("particles")
            if particles_tree and isinstance(particles_tree, ROOT.TTree):
                print(f"Found particles tree with {particles_tree.GetEntries()} entries")
            else:
                print("Warning: Could not find 'particles' tree in particles.root")
        else:
            print(f"Warning: Could not open {particles_file_path}")
    else:
        print(f"Warning: {particles_file_path} not found")


    canvas = ROOT.TCanvas("c1", "Tracking Performance Summary", 2400, 1800)
    canvas.Divide(4, 4)

    pt_bins = 100
    pt_max = 4 # Auto-range - let ROOT determine from data
    p_bins = 100
    p_max = 4 # Total momentum range
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

    # --- Third row: plots vs. total momentum (p) ---
    canvas.cd(9)
    ROOT.gPad.SetLeftMargin(0.15)
    h_den_p = ROOT.TH1F("h_den_p", "Total Tracks vs. p", p_bins, 0, p_max)
    h_num_p = ROOT.TH1F("h_num_p", "Good Tracks vs. p", p_bins, 0, p_max)
    track_tree.Draw("t_p>>h_den_p", "", "goff")
    track_tree.Draw("t_p>>h_num_p", "trackClassification == 1", "goff")
    efficiency_plot_p = ROOT.TEfficiency(h_num_p, h_den_p)
    efficiency_plot_p.SetTitle("Tracking Efficiency vs. p;p [GeV];Efficiency")
    efficiency_plot_p.Draw("AP")

    canvas.cd(10)
    ROOT.gPad.SetLeftMargin(0.15)
    h_clone_p = ROOT.TH1F("h_clone_p", "Clone Tracks vs. p", p_bins, 0, p_max)
    track_tree.Draw("t_p>>h_clone_p", "trackClassification == 2", "goff")
    clone_rate_plot_p = ROOT.TEfficiency(h_clone_p, h_den_p)
    clone_rate_plot_p.SetTitle("Clone Rate vs. p;p [GeV];Clone Rate")
    clone_rate_plot_p.Draw("AP")

    canvas.cd(11)
    ROOT.gPad.SetLeftMargin(0.15)
    h_fake_p = ROOT.TH1F("h_fake_p", "Fake Tracks vs. p", p_bins, 0, p_max)
    track_tree.Draw("t_p>>h_fake_p", "trackClassification == 3", "goff")
    fake_rate_plot_p = ROOT.TEfficiency(h_fake_p, h_den_p)
    fake_rate_plot_p.SetTitle("Fake Rate vs. p;p [GeV];Fake Rate")
    fake_rate_plot_p.Draw("AP")

    canvas.cd(12)
    ROOT.gPad.SetLeftMargin(0.15)
    p_hist = ROOT.TH1F("p_dist", "Reconstructed Track p Distribution;p [GeV];Number of Tracks", p_bins, 0, p_max)
    track_tree.Draw("t_p>>p_dist")
    p_hist.SetLineColor(4)
    p_hist.SetLineWidth(2)
    p_hist.Draw()

    # --- Fourth row: Reconstructable tracks and Simulated tracks plots ---
    # Pad 13: Reconstructable Tracks vs. pT
    canvas.cd(13)
    ROOT.gPad.SetLeftMargin(0.15)
    h_reconstructable_pt = h_den_pt.Clone("h_reconstructable_pt")
    h_reconstructable_pt.SetTitle("Reconstructable Tracks vs. pT;p_{T} [GeV];Number of Reconstructable Tracks")
    h_reconstructable_pt.SetLineColor(2)  # Red color
    h_reconstructable_pt.SetLineWidth(2)
    h_reconstructable_pt.Draw()

    # Pad 14: Reconstructable Tracks vs. eta
    canvas.cd(14)
    ROOT.gPad.SetLeftMargin(0.15)
    h_reconstructable_eta = h_den_eta.Clone("h_reconstructable_eta")
    h_reconstructable_eta.SetTitle("Reconstructable Tracks vs. #eta;#eta;Number of Reconstructable Tracks")
    h_reconstructable_eta.SetLineColor(2)  # Red color
    h_reconstructable_eta.SetLineWidth(2)
    h_reconstructable_eta.Draw()

    # Pad 15: Reconstructable Tracks vs. total momentum (p)
    canvas.cd(15)
    ROOT.gPad.SetLeftMargin(0.15)
    h_reconstructable_p = h_den_p.Clone("h_reconstructable_p")
    h_reconstructable_p.SetTitle("Reconstructable Tracks vs. p;p [GeV];Number of Reconstructable Tracks")
    h_reconstructable_p.SetLineColor(2)  # Red color
    h_reconstructable_p.SetLineWidth(2)
    h_reconstructable_p.Draw()

    # Pad 16: Simulated Tracks vs. eta (from particles.root)
    canvas.cd(16)
    ROOT.gPad.SetLeftMargin(0.15)
    if particles_tree:
        h_simulated_eta = ROOT.TH1F("h_simulated_eta", "Simulated Tracks vs. #eta;#eta;Number of Simulated Tracks", eta_bins, eta_min, eta_max)
        particles_tree.Draw("eta>>h_simulated_eta", "", "goff")
        h_simulated_eta.SetLineColor(8)  # Green color
        h_simulated_eta.SetLineWidth(2)
        h_simulated_eta.Draw()
    else:
        # Create empty histogram with message if particles.root not available
        h_empty = ROOT.TH1F("h_empty", "Simulated Tracks (particles.root not found);#eta;Number of Tracks", eta_bins, eta_min, eta_max)
        h_empty.Draw()

    # Pad 12: Could add another plot here if needed (currently empty)


    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    print(f"Combined plot saved as '{output_filename}.png'")
    root_file.Close()
    if particles_file:
        particles_file.Close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python create_combined_efficiency_and_pt_plot.py <output_filename_without_extension>")
        sys.exit(1)
    
    output_filename = sys.argv[1]
    create_combined_efficiency_and_pt_plot(output_filename)