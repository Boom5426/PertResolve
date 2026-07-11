# Figure 1a — AI schematic prompt

**Panel role:** the single most important definition panel. Show, with no results,
that gene-level perturbation collapses different variants of one gene into one
averaged response, whereas AllelePerturb keeps each protein-coding variant as an
independent, held-out prediction target.

**Target tool:** general raster image model (GPT-Image / Nano Banana).
**Aspect ratio:** ~4:3 (portrait-ish), two stacked rows.

**Palette (use exactly):** background pure white `#FFFFFF`; ink `#1A1A1A`;
wild-type / control neutral grey `#7A7A7A`; TP53 blue `#0072B2`; the three TP53
variants in three distinguishable hues — R175H `#0072B2`, R273C `#00A0C6`,
R248Q `#5AB4E6` (all in the blue family so they read as "same gene, different
alleles"); held-out target highlighted with orange accent `#E69F00`.

---

## Prompt (paste into the image model)

> A clean, flat 2D vector scientific schematic on a pure white background, Nature-
> Methods editorial style, thin 0.5 pt strokes, no 3D, no gradients, no drop
> shadows, no photorealism, generous white space. Two stacked rows separated by a
> thin horizontal divider.
>
> TOP ROW, labelled "Gene-level perturbation": on the left a small grey cluster of
> dots labelled "wild-type cells"; a thin arrow points right to a single blue cluster
> labelled "TP53 perturbation"; a second arrow points to ONE merged blue response
> cloud on the right. Inside that merged cloud, three faint overlapping ghosted
> labels "R175H", "R273C", "R248Q" are blended together to show the alleles are
> averaged away. A short caption under the row reads "one gene → one averaged
> response".
>
> BOTTOM ROW, labelled "Protein-coding variant-level perturbation": on the left the
> same grey "wild-type cells" cluster; three separate arrows fan out to three
> distinct small clusters in three shades of blue, labelled "TP53 R175H",
> "TP53 R273C", "TP53 R248Q"; each maps to its own SEPARATE response region on the
> right, three clouds that do NOT overlap. One of them, "TP53 R175H", is outlined in
> an orange dashed ring and tagged "held-out prediction target". A short caption
> under the row reads "each allele → a distinct response".
>
> Minimalist, lots of white space, small neat sans-serif labels, colour-blind-safe
> blues and grey only. No title, no numbers, no charts.

---

**Exact text to place/verify in Illustrator (raster text is unreliable):**
- Row titles: "Gene-level perturbation" / "Protein-coding variant-level perturbation"
- "wild-type cells" (both rows)
- Top: "TP53 perturbation"; ghost labels "R175H", "R273C", "R248Q"; caption "one gene → one averaged response"
- Bottom: "TP53 R175H", "TP53 R273C", "TP53 R248Q"; badge "held-out prediction target"; caption "each allele → a distinct response"

**Do NOT include (results leakage guard):** no mention of models, accuracy, failure,
PDS, measurement floor, or any performance claim. This is a task-definition panel only.
