
import ROOT
import sys
import os

def get_impact_profiles(directory, color, label, pt_bins=20):

    file_path = os.path.join(directory, "tracksummary_ckf.root")
    if not os.path.exists(file_path):
        print(f"[ERROR] Required file '{file_path}' was not found in {directory}.")
        return None, None
    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"[ERROR] Could not open the file '{file_path}'.")
        return None, None
    track_tree = root_file.Get("tracksummary")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        print(f"[ERROR] Could not find TTree 'tracksummary' in the file '{file_path}'.")
        root_file.Close()
        return None, None
    h_d0_profile = ROOT.TProfile(f"h_d0_profile_{label}", f"Average |d_{{0}}| vs p_{{T}} ({label});p_{{T}} [GeV];Average |d_{{0}}| [mm]", pt_bins, 0, 0)
    h_z0_profile = ROOT.TProfile(f"h_z0_profile_{label}", f"Average |z_{{0}}| vs p_{{T}} ({label});p_{{T}} [GeV];Average |z_{{0}}| [mm]", pt_bins, 0, 0)
    track_tree.Draw(f"abs(eLOC0_fit):t_pT>>h_d0_profile_{label}", "hasFittedParams", "goff")
    track_tree.Draw(f"abs(eLOC1_fit):t_pT>>h_z0_profile_{label}", "hasFittedParams", "goff")

    print(f"[INFO] {label} d0 entries: {h_d0_profile.GetEntries()}")
    print(f"[INFO] {label} z0 entries: {h_z0_profile.GetEntries()}")

    h_d0_profile.SetLineColor(color)
    h_d0_profile.SetLineWidth(3)
    h_d0_profile.SetMarkerColor(color)
    h_d0_profile.SetMarkerStyle(20)
    h_d0_profile.SetMarkerSize(1.2)
    h_d0_profile.SetMinimum(0)
    h_d0_profile.SetMaximum(1)

    h_z0_profile.SetLineColor(color)
    h_z0_profile.SetLineWidth(3)
    h_z0_profile.SetMarkerColor(color)
    h_z0_profile.SetMarkerStyle(20)
    h_z0_profile.SetMarkerSize(1.2)
    h_z0_profile.SetMinimum(0)
    h_z0_profile.SetMaximum(1)

    x_min = 0.0
    x_max = h_d0_profile.GetXaxis().GetXmax()
    h_d0_profile.GetXaxis().SetRangeUser(x_min, x_max)
    h_z0_profile.GetXaxis().SetRangeUser(x_min, x_max)

    h_d0_profile.SetDirectory(0)
    h_z0_profile.SetDirectory(0)

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
    legend_d0 = ROOT.TLegend(0.65, 0.75, 0.90, 0.90)
    legend_z0 = ROOT.TLegend(0.65, 0.75, 0.90, 0.90)
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
        print(f"[INFO] Processing folder: {folder}")
        h_d0, h_z0 = get_impact_profiles(folder, color, label)
        valid = True
        if h_d0 and h_z0:
            if h_d0.GetEntries() > 0:
                d0_profiles.append(h_d0)
            else:
                print(f"[WARNING] {label} d0 profile is empty.")
                valid = False
            if h_z0.GetEntries() > 0:
                z0_profiles.append(h_z0)
            else:
                print(f"[WARNING] {label} z0 profile is empty.")
                valid = False
            if valid:
                labels.append(label)
        else:
            print(f"[WARNING] Skipping {folder} due to missing or invalid data.")

    canvas.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    if d0_profiles:
        for i, h_d0 in enumerate(d0_profiles):
            h_d0.SetStats(0)
            draw_opt = "PE" if i == 0 else "PE SAME"
            h_d0.Draw(draw_opt)
            legend_d0.AddEntry(h_d0, labels[i], "lp")
        legend_d0.Draw()
    else:
        print("[ERROR] No valid d0 profiles to plot.")

    canvas.cd(2)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()
    if z0_profiles:
        for i, h_z0 in enumerate(z0_profiles):
            h_z0.SetStats(0)
            draw_opt = "PE" if i == 0 else "PE SAME"
            h_z0.Draw(draw_opt)
            legend_z0.AddEntry(h_z0, labels[i], "lp")
        legend_z0.Draw()
    else:
        print("[ERROR] No valid z0 profiles to plot.")

    if d0_profiles or z0_profiles:
        canvas.Update()
        canvas.SaveAs(f"{output_filename}.png")
        print(f"Combined impact parameter plot saved as '{output_filename}.png'")
    else:
        print("[ERROR] No valid profiles found in any directory. No plot was created.")

if __name__ == "__main__":
    main()