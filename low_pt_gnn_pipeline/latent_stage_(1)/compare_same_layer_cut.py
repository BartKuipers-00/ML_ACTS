#!/usr/bin/env python3
"""Compare same-layer cut variants on identical embeddings: none / |dr|<cut / same-(vol,layer).
Reports precision & recall so we can see the gain from a geometry-exact same-layer cut.
"""
import sys, glob, os
from pathlib import Path
import numpy as np, pandas as pd, torch
from tqdm import tqdm
from torch_geometric.loader import DataLoader

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(PIPELINE_ROOT.parent / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from low_pt_custom_utils.graph_utils import build_edges, edge_truth_labels
from acorn.stages.graph_construction.graph_construction_stage import EventDataset
from test_my_latent_model import load_model

CKPT = sys.argv[1]
if not Path(CKPT).is_absolute() and not Path(CKPT).exists():
    CKPT = str(PIPELINE_ROOT / CKPT)
NUM = int(sys.argv[2]) if len(sys.argv) > 2 else 50
KNN, RMAX, RGEO, DRCUT = 180, 0.14, 600.0, 10.0
TRUTH_DIR = PIPELINE_ROOT / 'data' / 'feature_store' / 'testset'

model, hparams = load_model(CKPT)
node_scales = hparams['node_scales']
input_dir = PIPELINE_ROOT / 'data' / 'feature_store'
testset = EventDataset(input_dir=str(input_dir), data_name='testset', num_events=NUM, hparams=hparams)
loader = DataLoader(testset, batch_size=1, num_workers=0, shuffle=False)

def vol_lay_for_batch(batch):
    """Per-node (volume_id, layer_id) packed into a single int key, via truth.csv join on hit_id."""
    eid = batch.event_id
    eid = eid[0] if hasattr(eid, '__getitem__') and not isinstance(eid, str) else eid
    try: eid = int(eid)
    except Exception: pass
    tcsv = TRUTH_DIR / (f"event{eid:09d}-truth.csv" if isinstance(eid, int) else f"{eid}-truth.csv")
    df = pd.read_csv(tcsv)
    key = (df.volume_id.astype(np.int64) * 1000 + df.layer_id.astype(np.int64)).values
    lut = dict(zip(df.hit_id.values, key))
    hid = batch.hit_id.cpu().numpy()
    return np.array([lut.get(h, -1) for h in hid], dtype=np.int64)

def metrics_for(pred_edges, true_edges, n_nodes, n_unique_true):
    pe, ey = edge_truth_labels(pred_edges, true_edges, num_nodes=n_nodes,
                               undirected=hparams.get("undirected", False))
    ey = ey.cpu().numpy().astype(bool)
    tpe = pe[:, ey]
    if tpe.shape[1] > 0:
        can = torch.stack([torch.min(tpe, 0)[0], torch.max(tpe, 0)[0]])
        tp = torch.unique(can, dim=1).shape[1]
    else:
        tp = 0
    fp = int((~ey).sum()); fn = n_unique_true - tp
    return tp, fp, fn

agg = {k: [0,0,0] for k in ('none','dr','geom')}
with torch.no_grad():
    for batch in tqdm(loader, desc="events"):
        emb = model.apply_embedding(batch)
        pred = build_edges(query=emb, database=emb, indices=None, r_max=RMAX, k_max=KNN, backend="FRNN")
        src, dst = pred
        hx, hy = batch.hit_x.float(), batch.hit_y.float()
        hz = batch.hit_z.float() * node_scales[2]
        d3 = torch.sqrt((hx[src]-hx[dst])**2 + (hy[src]-hy[dst])**2 + (hz[src]-hz[dst])**2)
        pred = pred[:, d3 <= RGEO]
        src, dst = pred
        # variant masks
        hr = batch.hit_r.float() * node_scales[0]
        dr = torch.abs(hr[src] - hr[dst])
        m_dr = dr >= DRCUT
        vl = torch.from_numpy(vol_lay_for_batch(batch))
        m_geom = vl[src] != vl[dst]
        true_edges = batch.track_edges
        n_unique_true = torch.unique(true_edges, dim=1).shape[1]
        n_nodes = emb.shape[0]
        for name, mask in (('none', None), ('dr', m_dr), ('geom', m_geom)):
            pe = pred if mask is None else pred[:, mask]
            tp, fp, fn = metrics_for(pe, true_edges, n_nodes, n_unique_true)
            agg[name][0]+=tp; agg[name][1]+=fp; agg[name][2]+=fn

print(f"\n{'variant':22s} {'TP':>9s} {'FP':>12s} {'FN':>6s} {'precision':>11s} {'recall':>9s}")
labels = {'none':'no same-layer cut','dr':'|dr|<10 (current)','geom':'same-(vol,layer) EXACT'}
for k in ('none','dr','geom'):
    tp,fp,fn = agg[k]
    p = tp/(tp+fp) if tp+fp else 0; r = tp/(tp+fn) if tp+fn else 0
    print(f"{labels[k]:22s} {tp:9d} {fp:12d} {fn:6d} {p*100:10.4f}% {r*100:8.4f}%")
