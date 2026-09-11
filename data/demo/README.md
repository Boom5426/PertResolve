# Self-contained resolution demo

`pertresolve_resolution_demo.npz` is a small, deterministic demonstration of
the public `pertresolve.resolution` API:

- 160 profiles: 32 controls and 32 profiles for each of four TP53 variants;
- 939 response features;
- four disjoint groups of eight profiles per condition, matching the demo depth;
- generated with a fixed seed from four public TP53 response vectors in the
  [PertResolve-Bench Hugging Face release](https://huggingface.co/datasets/Boom5426/PertResolve_Bench).

The file is intentionally small and runnable on a laptop. It is a derived
demonstration matrix, not a replacement for the paper's raw single-cell data
and not a paper-reproduction input. The output reports detection,
identification, an empirical split-half reproducibility reference and
model-ranking resolution as independent axes; it does not collapse them into a
single verdict. Run it with:

```bash
python examples/run_resolution_demo.py
```

The full benchmark metadata and larger processed arrays are linked from the
root README.
