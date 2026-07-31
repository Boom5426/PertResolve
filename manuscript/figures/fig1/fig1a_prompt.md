# Figure 1a — AI schematic prompt (pipeline-first; replaces the old concept panel)

**Panel role:** the entry panel of the whole paper. Show the AllelePerturb pipeline
end to end, so a reader enters Figure 2 already knowing (i) what is measured, (ii) what
the model does and does *not* see, (iii) how a prediction is scored, and (iv) the
parallel single-cell resolution diagnostic that decides whether the ground truth can
reward an allele model at all. **No results.** This panel replaces the old
gene-vs-variant concept panel and **subsumes the old workflow panel (old g)**; the
gene-vs-variant motivation survives, lightly, inside the discrimination axis (stage 3).

**Target tool:** general raster image model (GPT-Image / Nano Banana).
**Aspect ratio:** ~3:2 landscape (a main left-to-right row with a thin parallel lower
track). It will be placed in an enlarged panel-a slot at assembly (row-1 left, widened
into the space freed by removing old g).

**Palette (use exactly, matches the sibling panels):** background pure white `#FFFFFF`;
ink `#1A1A1A`; wild-type / control neutral grey `#7A7A7A`; the variant / TP53 blue
`#0072B2`, with two sibling shades `#00A0C6` and `#5AB4E6` (same blue family = "same
gene, different alleles"); held-out target and the verdict accent in orange `#E69F00`;
model and evaluation step boxes in a single muted blue-grey. No 3D, no gradients, no
shadows.

---

## Prompt (paste into the image model)

> A clean, flat 2D vector scientific schematic on a pure white background,
> Nature-Methods editorial style, thin 0.5 pt strokes, no 3D, no gradients, no drop
> shadows, generous white space. A single left-to-right pipeline of three main stages
> ending in a verdict box, with one thin parallel lower track running beneath the whole
> row and feeding up into the verdict.
>
> STAGE 1, header "Measurement": on the left, a small blue cluster of dots labelled
> "held-out variant, single cells" sits above a grey cluster labelled "wild-type cells".
> A short downward bracket collapses the two clusters into ONE small horizontal vector
> bar labelled "delta-v = pseudobulk target (ground truth)"; a tiny sub-label reads
> "mean(variant) - mean(WT)". This makes clear the target is the averaged profile built
> from the cells.
>
> STAGE 2, header "Model": from STAGE 1, an arrow carrying ONLY a small feature chip
> labelled "variant features (theta / ESM)" enters a rounded box labelled
> "feature -> response model". A short red-free grey note beside the arrow reads
> "not the variant identity, not its cells". The model box outputs a second small vector
> bar labelled "delta-hat-v (predicted profile)". Draw the feature chip, NOT cells,
> entering the model.
>
> STAGE 3, header "Evaluation": the predicted bar delta-hat-v is scored on two stacked
> mini-axes. Top mini-axis "Direction (Pearson-delta)": two short overlaid line profiles
> pointing the same way, caption "points the right way". Bottom mini-axis
> "Discrimination (PDS)": delta-hat-v on the left with an arrow to a small set of three
> candidate sibling profiles drawn as three separate blue-family clusters occupying
> distinct regions; the variant's own cluster is ringed in orange as the nearest match,
> the two siblings stay grey-blue, caption "nearest its own allele among siblings".
>
> STAGE 5, header "Verdict" (rightmost, orange-accented rounded box): a short stacked
> checklist "measurable?", "benchmarkable?", "worth modelling?" with a small branching
> arrow out to three outcomes "expand", "redesign", "exclude". This is the power-aware
> output of the pipeline.
>
> LOWER TRACK, header "Resolution diagnostic (single cell)": a thin full-width band
> beneath stages 1 to 3 that operates on the SAME cell clusters with no averaging. Two
> small side-by-side motifs: "detection = variant vs wild-type distribution" and
> "identification = variant vs sibling distributions", each drawn as two cell clouds with
> a small distance caliper between them. A thin arrow rises from this band into the
> Verdict box, with a caption "can the ground truth reward an allele model?".
>
> Keep stages evenly sized, thin connecting arrows left to right, small neat sans-serif
> labels, colour-blind-safe blues, grey and one orange accent only, minimalist, no charts
> with axes ticks, no numbers.

---

**Exact text to place / verify in Illustrator (raster text is unreliable):**
- Stage headers: "1 Measurement" / "2 Model" / "3 Evaluation" / "4 Resolution diagnostic (single cell)" / "5 Verdict"
- Stage 1: "held-out variant, single cells", "wild-type cells", "delta_v = pseudobulk target (ground truth)", sub-label "mean(variant) - mean(WT)"
- Stage 2: "variant features (theta / ESM)", grey note "not the variant identity, not its cells", "feature -> response model", "delta-hat_v (predicted profile)"
- Stage 3 top: "Direction (Pearson-delta)", "points the right way"
- Stage 3 bottom: "Discrimination (PDS)", "nearest its own allele among siblings"; three sibling clusters, own ringed in orange as "held-out target"
- Stage 4 (lower track): "detection = variant vs wild type", "identification = variant vs siblings", "same single cells, no averaging", "can the ground truth reward an allele model?"
- Stage 5: "measurable?", "benchmarkable?", "worth modelling?", outcomes "expand", "redesign", "exclude"

Use the Greek delta glyph (δ) and the hat (δ̂) in the final Illustrator text; the prompt
spells them out only because raster models mangle glyphs.

**Do NOT include (results-leakage guard):** no score value, no PDS/Pearson number, no
statement that any model or the benchmark succeeds or fails, no mention of "chance",
"measurement floor", "0.50", or any per-gene outcome. This is a task-and-pipeline
definition panel only; every result begins in Figure 2 and later.

---

## Assembly notes (step 2, after the raster is made)

- Render, overlay the exact text above in Illustrator, export vector-preserving
  `Fig1a.pdf` into this directory at >= 300 dpi effective size.
- Edit `fig1_assemble.tex`: place the new `Fig1a.pdf` in the panel-a slot, **delete the
  `Fig1g.pdf` node and its `g` panel letter**, and reflow so panel a takes part of the
  width freed by removing g (b-f keep their content). Panel letters become **a-f**.
- Recompile `lualatex fig1_assemble.tex`, then copy `fig1_composite.pdf` over
  `../../latex/figures/fig1.pdf` so the rendered figure matches the a-f caption already
  in the manuscript.
- This also resolves the current mismatch: the embedded `fig1.pdf` is still the old
  seven-panel a-g composite, while the caption already describes a-f.
