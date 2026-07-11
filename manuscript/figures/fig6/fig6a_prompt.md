# Figure 6a — AI schematic prompt

**Panel role:** define pilot-based rankability prediction. It predicts NOT the expression
response and NOT model performance, but whether a perturbation's ground truth will be
measurable enough to enter model comparison. No results.

**Target tool:** general raster image model. **Aspect ratio:** ~4:3 (left→right, 3 steps).

**Palette:** white background; ink `#1A1A1A`; pilot/variant blue `#5185C0`; grey `#9AA0A6`;
rankable green `#55966B`; below-floor muted red `#D98C8C`.

---

## Prompt (paste into the image model)

> A clean flat 2D vector scientific schematic, Nature-Methods style, pure white
> background, thin 0.5 pt strokes, no 3D, no shadows, no gradients. Left-to-right, three steps.
>
> STEP 1 "Small pilot": a tiny blue cluster of cells labelled "pilot: 25-50 cells per
> perturbation" with a small balanced grey "WT / control" cluster beside it.
>
> STEP 2 "Estimate measurement properties": two small readouts, "variant-to-WT effect
> size" and "within-variant split-half noise", combining into a chip "pilot signal-to-noise".
>
> STEP 3 "Predict rankability": a small classifier box outputs two labelled outcomes — a
> green tag "likely rankable -> benchmark it" and a red tag "likely below floor -> triage
> out / redesign".
>
> A single caption strip along the bottom: "can a small pilot identify perturbations whose
> ground truth will be measurable at higher depth?". Minimalist, tidy, small neat
> sans-serif labels, colour-blind-safe, no numbers, no charts.

---

**Exact text to place/verify in Illustrator:**
- Step titles: "Small pilot" / "Estimate measurement properties" / "Predict rankability"
- "pilot: 25-50 cells per perturbation", "WT / control", "variant-to-WT effect size",
  "within-variant split-half noise", "pilot signal-to-noise"
- outcome tags "likely rankable -> benchmark it" / "likely below floor -> triage out / redesign"
- caption strip (verbatim above)

**Do NOT include:** any AUROC value, cell-number cutoff, or model name. Use "predicted
rankability" / "triage" wording, never "guaranteed" or "required cell number".
