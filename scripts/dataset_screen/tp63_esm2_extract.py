#!/usr/bin/env python3
"""Mutation-aware ESM2-650M features for the 18 arrayed TP63 alleles of GSE311877.

Implements section 3 of docs/PREREG_GSE311877_FEATURE_TO_RESIDUAL_v1.md. Three
representations per variant, all differences against the wild-type forward pass:

    esm2_sitedelta    h_mut[p] - h_wt[p]                 at the mutated residue
    esm2_window16     mean over p +/- 16 of h_mut - h_wt  local context
    esm2_globaldelta  mean(h_mut) - mean(h_wt)            whole sequence

Numbering: the GSE311877 labels use neither UniProt isoform's coordinates. One
constant offset reconciles them and is asserted here against every wild-type
residue, so a wrong offset aborts instead of silently producing features for the
wrong site:

    label N = TAp63-alpha (Q9H3D4) residue N+39 = dNp63-alpha (Q9H3D4-2) residue N-55

Primary background is dNp63-alpha, the keratinocyte isoform driving the assay.
TAp63-alpha is available via --isoform as the declared sensitivity check; the two
are identical from TAp63-alpha residue 115 on, and every mutated site sits at
318-595, so only the global mean representation can differ.

Runs on the remote GPU box (conda env `Agent`): torch + fair-esm, ESM2-650M
weights already cached under ~/.cache/torch/hub/checkpoints. Both sequences are
under the 1022-residue limit, so no cropping is needed.

    python tp63_esm2_extract.py --out tp63_esm2_dnp63a.npz
    python tp63_esm2_extract.py --isoform TAp63a --out tp63_esm2_tap63a.npz
"""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

# UniProt Q9H3D4 (TAp63-alpha, 680 aa) and Q9H3D4-2 (dNp63-alpha, 586 aa),
# retrieved 2026-08-03. md5 asserted below so a silently edited constant fails.
SEQS = {
    "TAp63a": (
        "MNFETSRCATLQYCPDPYIQRFVETPAHFSWKESYYRSTMSQSTQTNEFLSPEVFQHIWD"
        "FLEQPICSVQPIDLNFVDEPSEDGATNKIEISMDCIRMQDSDLSDPMWPQYTNLGLLNSM"
        "DQQIQNGSSSTSPYNTDHAQNSVTAPSPYAQPSSTFDALSPSPAIPSNTDYPGPHSFDVS"
        "FQQSSTAKSATWTYSTELKKLYCQIAKTCPIQIKVMTPPPQGAVIRAMPVYKKAEHVTEV"
        "VKRCPNHELSREFNEGQIAPPSHLIRVEGNSHAQYVEDPITGRQSVLVPYEPPQVGTEFT"
        "TVLYNFMCNSSCVGGMNRRPILIIVTLETRDGQVLGRRCFEARICACPGRDRKADEDSIR"
        "KQQVSDSTKNGDGTKRPFRQNTHGIQMTSIKKRRSPDDELLYLPVRGRETYEMLLKIKES"
        "LELMQYLPQHTIETYRQQQQQQHQHLLQKQTSIQSPSSYGNSSPPLNKMNSMNKLPSVSQ"
        "LINPQQRNALTPTTIPDGMGANIPMMGTHMPMAGDMNGLSPTQALPPPLSMPSTSHCTPP"
        "PPYPTDCSIVSFLARLGCSSCLDYFTTQGLTTIYQIEHYSMDDLASLKIPEQFRHAIWKG"
        "ILDHRQLHEFSSPSHLLRTPSSASTVSVGSSETRGERVIDAVRFTLRQTISFPPRDEWND"
        "FNFDMDARRNKQQRIKEEGE"
    ),
    "dNp63a": (
        "MLYLENNAQTQFSEPQYTNLGLLNSMDQQIQNGSSSTSPYNTDHAQNSVTAPSPYAQPSS"
        "TFDALSPSPAIPSNTDYPGPHSFDVSFQQSSTAKSATWTYSTELKKLYCQIAKTCPIQIK"
        "VMTPPPQGAVIRAMPVYKKAEHVTEVVKRCPNHELSREFNEGQIAPPSHLIRVEGNSHAQ"
        "YVEDPITGRQSVLVPYEPPQVGTEFTTVLYNFMCNSSCVGGMNRRPILIIVTLETRDGQV"
        "LGRRCFEARICACPGRDRKADEDSIRKQQVSDSTKNGDGTKRPFRQNTHGIQMTSIKKRR"
        "SPDDELLYLPVRGRETYEMLLKIKESLELMQYLPQHTIETYRQQQQQQHQHLLQKQTSIQ"
        "SPSSYGNSSPPLNKMNSMNKLPSVSQLINPQQRNALTPTTIPDGMGANIPMMGTHMPMAG"
        "DMNGLSPTQALPPPLSMPSTSHCTPPPPYPTDCSIVSFLARLGCSSCLDYFTTQGLTTIY"
        "QIEHYSMDDLASLKIPEQFRHAIWKGILDHRQLHEFSSPSHLLRTPSSASTVSVGSSETR"
        "GERVIDAVRFTLRQTISFPPRDEWNDFNFDMDARRNKQQRIKEEGE"
    ),
}
MD5 = {"TAp63a": "7c3a66b8403223e802e1a24563190ea7",
       "dNp63a": "5d9b54fec18acb11ffa5ac6e03fcf79e"}
# label position -> sequence index offset, per the pre-registration
OFFSET = {"TAp63a": 39, "dNp63a": -55}

VARIANTS = ["R279Q", "R279S", "G293H", "G297L", "R304Q", "R304T", "G310Q", "R311I",
            "L514D", "L514F", "C522D", "C522G", "T527I", "L531E", "L531R",
            "L548N", "I550D", "H556A"]
WINDOW = 16


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--isoform", choices=sorted(SEQS), default="dNp63a",
                    help="sequence background (default: dNp63a, the keratinocyte isoform)")
    ap.add_argument("--out", required=True, help="output .npz")
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    wt = SEQS[args.isoform]
    got = hashlib.md5(wt.encode()).hexdigest()
    if got != MD5[args.isoform]:
        raise SystemExit(f"{args.isoform} sequence md5 {got} != expected {MD5[args.isoform]}")

    off = OFFSET[args.isoform]
    sites = {}
    for v in VARIANTS:
        aa_wt, pos_label, aa_mut = v[0], int(v[1:-1]), v[-1]
        idx = pos_label + off                       # 1-based index into `wt`
        if not 1 <= idx <= len(wt):
            raise SystemExit(f"{v}: mapped index {idx} outside 1..{len(wt)}")
        if wt[idx - 1] != aa_wt:
            raise SystemExit(f"{v}: expected {aa_wt} at {args.isoform} index {idx}, "
                             f"found {wt[idx - 1]}. Offset is wrong; aborting.")
        sites[v] = idx
    print(f"{args.isoform}: {len(wt)} aa, all {len(VARIANTS)} wild-type residues verified")

    import torch
    import esm

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(args.device)
    bc = alphabet.get_batch_converter()
    layer = model.num_layers

    @torch.no_grad()
    def reps(seq: str) -> np.ndarray:
        _, _, toks = bc([("x", seq)])
        out = model(toks.to(args.device), repr_layers=[layer])["representations"][layer]
        # drop BOS/EOS so row i is residue i+1
        return out[0, 1:len(seq) + 1].float().cpu().numpy()

    h_wt = reps(wt)
    site, win, glob = {}, {}, {}
    for v in VARIANTS:
        idx = sites[v]
        mut = wt[:idx - 1] + v[-1] + wt[idx:]
        assert len(mut) == len(wt) and mut[idx - 1] == v[-1]
        h_mut = reps(mut)
        d = h_mut - h_wt
        lo, hi = max(0, idx - 1 - WINDOW), min(len(wt), idx + WINDOW)
        site[v] = d[idx - 1]
        win[v] = d[lo:hi].mean(axis=0)
        glob[v] = h_mut.mean(axis=0) - h_wt.mean(axis=0)
        print(f"  {v:<7} idx={idx:<4} |site|={np.linalg.norm(site[v]):8.3f} "
              f"|win|={np.linalg.norm(win[v]):8.3f} |glob|={np.linalg.norm(glob[v]):8.3f}")

    order = np.array(VARIANTS)
    np.savez_compressed(
        args.out,
        variants=order,
        esm2_sitedelta=np.stack([site[v] for v in VARIANTS]),
        esm2_window16=np.stack([win[v] for v in VARIANTS]),
        esm2_globaldelta=np.stack([glob[v] for v in VARIANTS]),
        meta=json.dumps({"isoform": args.isoform, "length": len(wt), "md5": got,
                         "offset_label_to_index": off, "window": WINDOW,
                         "model": "esm2_t33_650M_UR50D", "layer": layer,
                         "sites": sites}),
    )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
