#!/usr/bin/env python3
"""Top-3 #3: mutation-aware ESM2-650M representations (robustness control).

For each variant we build the mutant sequence from the WT and extract, from ESM2-650M
final-layer per-residue hidden states:
  esm2_global       : mean over residues of the MUTANT sequence (代际对照 for ESM-1v)
  esm2_globaldelta  : mean(mutant) - mean(WT)          (encode change relative to WT)
  esm2_sitedelta    : h_mut[p] - h_wt[p] at the mutated residue(s), averaged  (focused local semantics)
  esm2_window16     : mean over p +/- 16 of (h_mut - h_wt)                     (local context)
Long proteins (JAK1, 1154 aa > ESM2's 1022-residue limit) are cropped to a 1022-residue
window centred on the mutation (global reps then span the crop; documented). Multi-substitution
variants (e.g. GATA1 C204R,V205A) apply all substitutions and average their site deltas.
Outputs 4 feature npz keyed {GENE}__{variant} (+ {GENE}__WT), ready for the 6-head grid.
"""
import json, re, os, numpy as np, torch, esm, sys
sys.path.insert(0, "/data/boom/NUS/VCCompass/unified")
import harness as H

BASE = "/data/boom/NUS/VCCompass"
OUT = f"{BASE}/esm2_control"; os.makedirs(OUT, exist_ok=True)
MAXLEN = 1022; WIN = 16
WT = json.load(open(f"{BASE}/wt_seqs.json"))
MUT_RE = re.compile(r'^([A-Z])(\d+)([A-Z])$')

model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
model = model.eval().cuda()
bc = alphabet.get_batch_converter()
LAYER = model.num_layers


@torch.no_grad()
def per_residue(seq):
    _, _, toks = bc([("x", seq)])
    out = model(toks.cuda(), repr_layers=[LAYER])
    rep = out["representations"][LAYER][0]          # (len+2, D) incl BOS/EOS
    return rep[1:len(seq) + 1].float().cpu().numpy()  # (len, D)


def parse(variant):
    subs = []
    for part in variant.split(','):
        m = MUT_RE.match(part.strip())
        if not m:
            return None
        subs.append((m.group(1), int(m.group(2)), m.group(3)))
    return subs


reps = {r: {} for r in ['esm2_global', 'esm2_globaldelta', 'esm2_sitedelta', 'esm2_window16']}
D = model.embed_dim
for g in H.GENES:
    wt_seq = WT[g]
    X, lab = H.load_gene(g); lab = np.asarray(lab)
    variants = [v for v in np.unique(lab) if v not in ('WT', 'wt', 'WT_control')]
    # WT reference (crop-independent global uses the full or window-per-variant; for WT keys store full-seq global)
    ok = skip = 0
    for v in variants:
        subs = parse(v)
        if subs is None:
            skip += 1; continue
        # validate + build mutant
        mut = list(wt_seq); positions = []; bad = False
        for wa, p, ma in subs:
            if p < 1 or p > len(wt_seq) or wt_seq[p - 1] != wa:
                bad = True; break
            mut[p - 1] = ma; positions.append(p)
        if bad:
            skip += 1; continue
        mut = ''.join(mut)
        # crop long proteins around the mutation centre
        if len(wt_seq) > MAXLEN:
            c = int(np.mean(positions)); half = MAXLEN // 2
            s = max(0, min(c - half, len(wt_seq) - MAXLEN)); e = s + MAXLEN
        else:
            s, e = 0, len(wt_seq)
        wt_c, mut_c = wt_seq[s:e], mut[s:e]
        pos_c = [p - 1 - s for p in positions]   # 0-indexed in crop
        hw = per_residue(wt_c); hm = per_residue(mut_c)
        gd = hm.mean(0) - hw.mean(0)
        site = np.mean([hm[pc] - hw[pc] for pc in pos_c], axis=0)
        win = []
        for pc in pos_c:
            lo, hi = max(0, pc - WIN), min(len(mut_c), pc + WIN + 1)
            win.append((hm[lo:hi] - hw[lo:hi]).mean(0))
        win = np.mean(win, axis=0)
        k = f"{g}__{v}"
        reps['esm2_global'][k] = hm.mean(0).astype(np.float32)
        reps['esm2_globaldelta'][k] = gd.astype(np.float32)
        reps['esm2_sitedelta'][k] = site.astype(np.float32)
        reps['esm2_window16'][k] = win.astype(np.float32)
        ok += 1
    # WT keys: global = WT global mean (cropped to first MAXLEN if long); deltas = 0
    wt_c = wt_seq[:MAXLEN]
    hw = per_residue(wt_c)
    reps['esm2_global'][f"{g}__WT"] = hw.mean(0).astype(np.float32)
    for r in ['esm2_globaldelta', 'esm2_sitedelta', 'esm2_window16']:
        reps[r][f"{g}__WT"] = np.zeros(D, np.float32)
    print(f"{g}: {ok} variants embedded, {skip} skipped (unparseable/validation)", flush=True)

for r, d in reps.items():
    np.savez_compressed(f"{OUT}/{r}.npz", **d)
    print(f"saved {OUT}/{r}.npz  ({len(d)} keys, dim {D})")
