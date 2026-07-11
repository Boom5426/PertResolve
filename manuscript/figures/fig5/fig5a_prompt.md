# Figure 5a — AI schematic prompt

**Panel role:** define the oracle discrimination ceiling as an *independent second
measurement* of the same variant (NOT a model trained on labels). It estimates how
much allele discrimination the measurement itself permits at the available depth.
No results.

**Target tool:** general raster image model. **Aspect ratio:** ~4:3 (left→right, 3 steps).

**Palette:** white background; ink `#1A1A1A`; variant blue `#5185C0`; WT/candidates grey
`#9AA0A6`; accent (oracle path) `#55966B`.

---

## Prompt (paste into the image model)

> A clean flat 2D vector scientific schematic, Nature-Methods style, pure white
> background, thin 0.5 pt strokes, no 3D, no shadows, no gradients. Left-to-right, three
> numbered steps.
>
> STEP 1 "Split the variant's cells": a single blue cluster of cells labelled
> "variant v" is divided by a dashed line into two equal halves, "half 1" and "half 2".
>
> STEP 2 "Two independent pseudobulks": an arrow from half 1 leads to a small profile
> labelled "measured target δ(v)"; an arrow from half 2 leads to an identical-looking
> profile labelled "oracle prediction (2nd measurement)", drawn in green to mark it as
> the prediction path.
>
> STEP 3 "Score in the same harness": both feed a box labelled "PDS harness" alongside a
> row of small grey candidate profiles "other held-out variants"; an output points to a
> gauge/dial labelled "oracle discrimination ceiling".
>
> A single caption strip along the bottom reads: "the oracle prediction is a genuine
> second measurement of the same variant, so its score is the discrimination attainable
> at this depth". Minimalist, tidy, small neat sans-serif labels, colour-blind-safe, no
> numbers, no charts.

---

**Exact text to place/verify in Illustrator:**
- Step titles: "Split the variant's cells" / "Two independent pseudobulks" / "Score in the same harness"
- "variant v", "half 1", "half 2", "measured target δ(v)", "oracle prediction (2nd measurement)"
- "PDS harness", "other held-out variants", "oracle discrimination ceiling"
- caption strip (verbatim above)

**Do NOT include:** any ceiling value, any model name, any claim about which genes are
measurement-limited. This panel only *defines* the oracle; the per-gene values are 5b.
