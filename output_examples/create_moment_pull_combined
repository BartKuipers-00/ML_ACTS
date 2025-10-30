import ROOT
import sys
import os

def get_momentum_pull_profile(directory, color, label, pt_bins=20):
    file_path = os.path.join(directory, "tracksummary_ckf.root")
    if not os.path.exists(file_path):
        print(f"[ERROR] Required file '{file_path}' was not found in {directory}.")
        return None

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"[ERROR] Could not open the file '{file_path}'.")
        return None

    track_tree = root_file.Get("tracksummary")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        print(f"[ERROR] Could not find TTree 'tracksummary' in the file '{file_path}'.")
        root_file.Close()
        return None

    h_p_pull_profile = ROOT.TProfile(f"h_p_pull_profile_{label}",
                                     f"Fractional Momentum Pull vs p_{{T}} ({label});p_{{T,truth}} [GeV];|p_{{truth}} - p_{{fitted}}| / p_{{truth}}",
                                     pt_bins, 0, 0)
    draw_expr = "abs(t_p - 1.0/abs(eQOP_fit))/t_p:t_pT>>h_p_pull_profile_" + label
    track_tree.Draw(draw_expr, "hasFittedParams && abs(eQOP_fit) > 1e-6", "goff")

    h_p_pull_profile.SetLineColor(color)
    h_p_pull_profile.SetLineWidth(3)
    h_p_pull_profile.SetMarkerColor(color)
    h_p_pull_profile.SetMarkerStyle(20)
    h_p_pull_profile.SetMarkerSize(1.2)
    h_p_pull_profile.SetMinimum(0.0)
    h_p_pull_profile.SetMaximum(0.2)
    h_p_pull_profile.SetDirectory(0)

    root_file.Close()
    return h_p_pull_profile

def main():
    if len(sys.argv) < 3:
        print("Usage: python create_momentum_pull_combined.py <output_filename_without_extension> <folder1> [<folder2> ...]")
        sys.exit(1)

    output_filename = sys.argv[1]
    folders = sys.argv[2:]

    ROOT.gStyle.SetPalette(ROOT.kBird)
    canvas = ROOT.TCanvas("c1", "Fractional Momentum Pull vs pT (Multi-Species)", 1200, 900)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gPad.SetGrid()

    colors = [2, 4, 8, 6, 7, 9, 46, 38]
    legend = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    profiles = []
    labels = []
    for idx, folder in enumerate(folders):
        label = os.path.basename(folder).replace("_files", "")
        color = colors[idx % len(colors)]
        prof = get_momentum_pull_profile(folder, color, label)
        if prof:
            profiles.append(prof)
            labels.append(label)

    if profiles:
        # Disable statistics box for all profiles
        for prof in profiles:
            prof.SetStats(0)
        profiles[0].Draw("PE")
        legend.AddEntry(profiles[0], labels[0], "lp")
        for i in range(1, len(profiles)):
            profiles[i].Draw("PE SAME")
            legend.AddEntry(profiles[i], labels[i], "lp")
        legend.Draw()

        # Add horizontal reference lines
        x_min = profiles[0].GetXaxis().GetXmin()
        x_max = profiles[0].GetXaxis().GetXmax()
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
        print(f"Combined momentum pull plot saved as '{output_filename}.png'")
    else:
        print("[ERROR] No valid momentum pull profiles to plot.")

if __name__ == "__main__":
    main()