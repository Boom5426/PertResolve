# Fig 1 redesign — pipeline-first (design spec; render pending)

Goal (reviewer points 1–4): make the **input → output → evaluation → power-aware
verdict** pipeline the entry point of the figure and the Results, and show the
**single-cell vs mean** relationship honestly. Text (main .tex caption + Results
opening) is already updated to this spec; only the rendered panels lag.

## Panel plan
- **a (NEW): the AllelePerturb pipeline** — replaces the old "gene-vs-variant
  concept" panel and subsumes the old **g** workflow panel. Full spec below.
- b, c, d, e, f: **unchanged** (coverage / depth / θ features / eval axes / splits).
- **g: removed** (its content is folded into the new a).

## Panel a — pipeline layout (left → right, single row with a parallel lower track)

```
 [1 MEASUREMENT]              [2 MODEL]                 [3 EVALUATION]              [5 VERDICT]
 held-out variant v                                     ┌ Direction (Pearson-δ)     power-aware
 its single cells      variant features θ / ESM         │   δ̂ points the right way   triage:
 (+ WT cells)   ──────►  (NOT identity / NOT cells) ──►  │                            • measurable?
   │ cell cloud                    │                     └ Discrimination (PDS)       • benchmarkable?
   │ mean_v − mean_WT              ▼                        δ̂ nearest its own δ        • worth modeling?
   ▼                          model → δ̂_v  ─────────────►  among sibling variants     ──► expand /
 δ_v  = pseudobulk TARGET                                                                redesign /
 (ground truth, a G-vector)                                                              exclude
   │
   └──[4 RESOLUTION DIAGNOSTIC — single-cell, parallel track]──────────────────────►
        SAME single cells (no averaging):
          detection      = v vs WT distribution      (energy distance / cell classifier)
          identification = v vs sibling distributions (energy distance / cell classifier)
        → "can the ground truth reward an allele model at all?"  feeds the verdict (5)
```

## Key requirements (must be visually unambiguous)
1. **Input to the model = variant features (θ / ESM), not single cells.** Draw the
   feature vector entering the model box; explicitly *not* the cells.
2. **Prediction target = the pseudobulk mean δ_v** built from single cells
   (mean_v − mean_WT). Show the cells collapsing to one vector.
3. **The single-cell distributions feed the resolution diagnostic (track 4), not the
   model.** This is the legitimate "single-cell" contribution — draw it as a parallel
   lower track operating on the cell clouds directly (no averaging).
4. **Two evaluation axes** (direction Pearson-δ; discrimination PDS = retrieval of the
   variant among its siblings). PDS: show δ̂ compared against a small set of candidate
   sibling δ's, nearest-match = correct.
5. **Verdict box (5)** is the paper's significance: measurable / benchmarkable / worth
   modeling → expand / redesign / exclude. This is the power-aware output.
6. Embed the gene-vs-variant motivation lightly (siblings occupying distinct regions of
   response space) inside the discrimination axis, so the old concept panel is not lost.

## Rendering notes
- Style: match existing panels (#5185C0 palette; vector PDF; nm_style). Panel a is a
  schematic (like the old Fig1a / `fig1a_prompt.md`), not a data plot.
- Update `fig1_assemble.tex`: put the new `Fig1a.pdf` (pipeline) in the a slot, **drop
  the `Fig1g.pdf` node**, keep b–f. Re-emit `fig1_composite.pdf` → `../latex/figures/fig1.pdf`.
- The main .tex caption (panel a) and the Results opening already describe exactly this
  pipeline; render to match that wording.
