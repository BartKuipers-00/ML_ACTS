import ROOT
import sys
import os

def get_impact_profile(directory, color, label, pt_bins=20, pt_min=0, pt_max=4):
    file_path = os.path.join(directory, "tracksummary_ckf.root")
    if not os.path.exists(file_path):
        print(f"Error: Required file '{file_path}' was not found in {directory}.")
        return None, None

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"Error: Could not open the file '{file_path}'.")
        return None, None

    track_tree = root_file.Get("tracksummary")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        print("Error: Could not find TTree 'tracksummary' in the file.")
        root_file.Close()
        return None, None

    h_d0_profile = ROOT.TProfile(f"h_d0_profile_{label}", f"{label}: Average |d0| vs pT", pt_bins, pt_min, pt_max)
    h_z0_profile = ROOT.TProfile(f"h_z0_profile_{label}", f"{label}: Average |z0| vs pT", pt_bins, pt_min, pt_max)
    n_d0 = track_tree.Draw(f"abs(eLOC0_fit):t_pT>>h_d0_profile_{label}", "hasFittedParams", "goff")
    n_z0 = track_tree.Draw(f"abs(eLOC1_fit):t_pT>>h_z0_profile_{label}", "hasFittedParams", "goff")
    if n_d0 == 0 or n_z0 == 0:
        print(f"Warning: No entries for {label}")
        root_file.Close()
        return None, None

    h_d0_profile.SetLineColor(color)
    h_d0_profile.SetMarkerColor(color)
    h_d0_profile.SetLineWidth(3)
    h_d0_profile.SetMarkerStyle(20)
    h_d0_profile.SetMarkerSize(1.2)

    h_z0_profile.SetLineColor(color)
    h_z0_profile.SetMarkerColor(color)
    h_z0_profile.SetLineWidth(3)
    h_z0_profile.SetMarkerStyle(20)
    h_z0_profile.SetMarkerSize(1.2)

    root_file.Close()
    return h_d0_profile, h_z0_profile

def main():
    if len(sys.argv) < 3:
        print("Usage: python create_impact_parameter_combined.py <output_filename_without_extension> <folder1> [<folder2> ...]")
        sys.exit(1)

    output_filename = sys.argv[1]
    folders = sys.argv[2:]

    ROOT.gStyle.SetPalette(ROOT.kBird)
    canvas = ROOT.TCanvas("c1", "Impact Parameter Resolution vs pT (Multi-Species)", 1600, 600)
    canvas.Divide(2, 1)

    colors = [2, 4, 8, 6, 7, 9, 46, 38]
    legend_d0 = ROOT.TLegend(0.15, 0.75, 0.45, 0.90)
    legend_z0 = ROOT.TLegend(0.15, 0.75, 0.45, 0.90)
    legend_d0.SetBorderSize(0)
    legend_d0.SetFillStyle(0)
    legend_z0.SetBorderSize(0)
    legend_z0.SetFillStyle(0)

    d0_profiles = []
    z0_profiles = []
    labels = []
    for idx, folder in enumerate(folders):
        label = os.path.basename(folder).replace("_files", "")
        color = colors[idx % len(colors)]
        h_d0, h_z0 = get_impact_profile(folder, color, label)
        if h_d0 and h_z0:
            d0_profiles.append(h_d0)
            z0_profiles.append(h_z0)
            labels.append(label)

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    if d0_profiles:
        d0_profiles[0].SetTitle("Average |d_{0}| vs p_{T};p_{T} [GeV];Average |d_{0}| [mm]")
        d0_profiles[0].SetMinimum(0)
        d0_profiles[0].SetMaximum(0.3)
        d0_profiles[0].Draw("PE")
        legend_d0.AddEntry(d0_profiles[0], labels[0], "lp")
        for i in range(1, len(d0_profiles)):
            d0_profiles[i].Draw("PE SAME")
            legend_d0.AddEntry(d0_profiles[i], labels[i], "lp")
        legend_d0.Draw()

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    if z0_profiles:
        z0_profiles[0].SetTitle("Average |z_{0}| vs p_{T};p_{T} [GeV];Average |z_{0}| [mm]")
        z0_profiles[0].SetMinimum(0)
        z0_profiles[0].SetMaximum(0.3)
        z0_profiles[0].Draw("PE")
        legend_z0.AddEntry(z0_profiles[0], labels[0], "lp")
        for i in range(1, len(z0_profiles)):
            z0_profiles[i].Draw("PE SAME")
            legend_z0.AddEntry(z0_profiles[i], labels[i], "lp")
        legend_z0.Draw()

    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    print(f"Combined impact parameter plot saved as '{output_filename}.png'")

if __name__ == "__main__":
    main()