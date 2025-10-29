import ROOT
import sys
import os

def get_efficiency_histogram(directory, color, label, pt_bins=100, pt_max=4):
    file_path = os.path.join(directory, "tracksummary_ckf.root")
    if not os.path.exists(file_path):
        print(f"Error: Required file '{file_path}' was not found in {directory}.")
        return None

    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"Error: Could not open the file '{file_path}'.")
        return None

    track_tree = root_file.Get("tracksummary")
    if not track_tree or not isinstance(track_tree, ROOT.TTree):
        print("Error: Could not find TTree 'tracksummary' in the file.")
        root_file.Close()
        return None

    h_den_pt = ROOT.TH1F(f"h_den_pt_{label}", f"Total Tracks vs. pT ({label})", pt_bins, 0, pt_max)
    h_num_pt = ROOT.TH1F(f"h_num_pt_{label}", f"Good Tracks vs. pT ({label})", pt_bins, 0, pt_max)
    track_tree.Draw("t_pT>>h_den_pt_" + label, "", "goff")
    track_tree.Draw("t_pT>>h_num_pt_" + label, "trackClassification == 1", "goff")
    efficiency_plot_pt = ROOT.TEfficiency(h_num_pt, h_den_pt)
    efficiency_plot_pt.SetLineColor(color)
    efficiency_plot_pt.SetMarkerColor(color)
    efficiency_plot_pt.SetTitle("Tracking Efficiency vs. pT;p_{T} [GeV];Efficiency")
    root_file.Close()
    return efficiency_plot_pt

def main():
    if len(sys.argv) < 3:
        print("Usage: python create_eff_hist_multi_species.py <output_filename_without_extension> <folder1> [<folder2> ...]")
        sys.exit(1)

    output_filename = sys.argv[1]
    folders = sys.argv[2:]

    ROOT.gStyle.SetPalette(ROOT.kBird)
    canvas = ROOT.TCanvas("c1", "Tracking Efficiency vs. pT (Multi-Species)", 1200, 900)
    colors = [2, 4, 8, 6, 7, 9, 46, 38]  # Red, Blue, Green, etc.

    legend = ROOT.TLegend(0.6, 0.15, 0.88, 0.35)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    effs = []
    labels = []
    for idx, folder in enumerate(folders):
        label = os.path.basename(folder).replace("_files", "")
        color = colors[idx % len(colors)]
        eff = get_efficiency_histogram(folder, color, label)
        if eff:
            effs.append(eff)
            labels.append(label)
    if effs:
        effs[0].Draw("AP")
        ROOT.gPad.Update()
        ROOT.gPad.SetTopMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.12)
        ROOT.gPad.SetLeftMargin(0.12)
        ROOT.gPad.SetRightMargin(0.05)
        # Set y-axis range for efficiency
        effs[0].GetPaintedGraph().GetYaxis().SetRangeUser(0.4, 1.0)
        legend.AddEntry(effs[0], labels[0], "lp")
        for i in range(1, len(effs)):
            effs[i].Draw("P SAME")
            legend.AddEntry(effs[i], labels[i], "lp")
        legend.Draw()
    else:
        print("No efficiency histograms to plot.")

    canvas.Update()
    canvas.SaveAs(f"{output_filename}.png")
    print(f"Combined efficiency plot saved as '{output_filename}.png'")


if __name__ == "__main__":
    main()