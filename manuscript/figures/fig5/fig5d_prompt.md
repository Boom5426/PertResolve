# Figure 5d — AI schematic prompt

**Panel role:** explain how synthetic predictors of *known, graded* quality are built by
interpolating between the gene-mean profile and the variant-specific profile, so that
the true model ordering is known by construction (the basis for the 5f resolution test).
No results.

**Target tool:** general raster image model. **Aspect ratio:** ~16:9 (left→right).

**Palette:** white background; ink `#1A1A1A`; gene-mean grey `#9AA0A6`; variant blue
`#5185C0`; a 5-step gradient from grey to blue for the α series.

---

## Prompt (paste into the image model)

> A clean flat 2D vector scientific schematic, Nature-Methods style, pure white
> background, thin strokes, no 3D, no shadows. Left-to-right.
>
> LEFT: two reference profiles stacked — a grey one labelled "gene-mean δ (direction
> only, α = 0)" and a blue one labelled "variant-specific δ (α = 1)".
>
> MIDDLE: a formula chip "δ̂(α) = (1−α)·gene-mean + α·variant", and below it a row of
> five small predictor profiles that morph gradually from grey to blue, tagged
> "α = 0", "0.25", "0.5", "0.75", "1", getting closer to the blue variant profile as α
> increases.
>
> RIGHT: a short vertical "true quality" ladder with the five predictors ordered
> P(α=0) at the bottom up to P(α=1) at the top, and a small label "true order known by
> construction".
>
> Minimalist, tidy, small neat sans-serif labels, colour-blind-safe grey-to-blue only,
> no numbers besides the α tags, no charts.

---

**Exact text to place/verify in Illustrator:**
- "gene-mean δ (direction only, α = 0)", "variant-specific δ (α = 1)"
- formula "δ̂(α) = (1−α)·gene-mean + α·variant"
- α tags "0 / 0.25 / 0.5 / 0.75 / 1"
- "true quality", "true order known by construction"

**Do NOT include:** any PDS value, any real dataset, any claim about recovery success.
This panel only defines the synthetic predictor family; recovery is 5f.
