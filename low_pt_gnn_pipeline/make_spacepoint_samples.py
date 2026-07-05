#!/usr/bin/env python
"""Generate the two-panel spacepoint sample plot (cf. hit_density.ipynb cell 10)
for one event from each of data_40, data_100, data_200."""
from pathlib import Path
import torch
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent

DATASETS  = ['data_10', 'data_40', 'data_100', 'data_200']
SUBSET    = 'testset'
EVENT_IDX = 4   # which event file (sorted) to sample from each dataset


def make_plot(graph, out_path):
    hx = graph.hit_x.numpy()
    hy = graph.hit_y.numpy()
    z  = graph.hit_z.numpy()

    fontsize = 15
    label_fs = fontsize + 8
    tick_fs  = fontsize + 1
    lim      = 1050

    BASE_FS  = 14.0
    fs_scale = fontsize / BASE_FS

    plot_h  = 7.0
    plot_w1 = 7.0
    plot_w2 = plot_h * (1200 / (2 * lim))

    left_m  = 1.0 * fs_scale
    right_m = 0.3
    gap     = 0.5
    top_m   = 0.6 * fs_scale
    bot_m   = 0.8 * fs_scale

    fig_w = left_m + plot_w1 + gap + plot_w2 + right_m
    fig_h = top_m + plot_h + bot_m

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=300)

    def to_frac(x, total): return x / total

    ax1 = fig.add_axes([to_frac(left_m,             fig_w), to_frac(bot_m, fig_h),
                        to_frac(plot_w1,             fig_w), to_frac(plot_h, fig_h)])
    ax2 = fig.add_axes([to_frac(left_m+plot_w1+gap, fig_w), to_frac(bot_m, fig_h),
                        to_frac(plot_w2,             fig_w), to_frac(plot_h, fig_h)])
    ax2.sharey(ax1)

    # Left: down-beampipe (x-y)
    ax1.scatter(hx, hy, s=3, alpha=1, linewidths=0.7, color='steelblue', rasterized=True)
    ax1.set_xlim(-lim, lim)
    ax1.set_ylim(-lim, lim)
    ax1.set_xlabel('x [mm]', fontsize=label_fs)
    ax1.set_ylabel('y [mm]', fontsize=label_fs)
    ax1.set_title('Down-beampipe', fontsize=fontsize+10)
    ax1.tick_params(labelsize=tick_fs)
    ax1.xaxis.set_major_locator(plt.MultipleLocator(200))
    ax1.yaxis.set_major_locator(plt.MultipleLocator(200))

    # Right: longitudinal (z-y)
    ax2.scatter(z, hy, s=3, alpha=1, linewidths=0.7, color='steelblue', rasterized=True)
    ax2.set_xlim(-600, 600)
    ax2.set_xlabel('z [mm]', fontsize=label_fs)
    ax2.set_title('Longitudinal', fontsize=fontsize+10)
    ax2.tick_params(labelsize=tick_fs, labelleft=False)
    ax2.xaxis.set_major_locator(plt.MultipleLocator(200))

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved: {out_path}  ({hx.shape[0]} hits)')


for ds in DATASETS:
    files = sorted((SCRIPT_DIR / ds / 'feature_store' / SUBSET).glob('*.pyg'))
    if not files:
        print(f'[skip] no .pyg files in {ds}/feature_store/{SUBSET}')
        continue
    idx = min(EVENT_IDX, len(files) - 1)
    graph = torch.load(files[idx], weights_only=False)
    eid = getattr(graph, 'event_id', files[idx].stem)
    out_path = SCRIPT_DIR / ds / 'visuals' / f'spacepoints_pyg_{eid}.pdf'
    make_plot(graph, out_path)
