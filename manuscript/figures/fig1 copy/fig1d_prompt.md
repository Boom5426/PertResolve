# Figure 1d — AI schematic prompt (feature concept; pairs with fig1d_theta_pca.pdf)

**Panel role:** explain why held-out variant prediction is possible — the model input
is a protein-level feature vector, not a one-hot variant ID. The data-direct
"shared feature space" scatter is generated separately (`fig1d_theta_pca.pdf`) and
placed as the right-hand sub-part; this prompt covers the left/middle schematic.

**Target tool:** general raster image model. **Aspect ratio:** ~16:9 (wide, left→right).

**Palette:** white background; ink `#1A1A1A`; TP53 blue `#0072B2`; feature icons in
muted grey `#7A7A7A` with a single orange accent `#E69F00` for the hotspot feature.

---

## Prompt (paste into the image model)

> A clean flat 2D vector scientific schematic, Nature-Methods style, pure white
> background, thin strokes, no 3D, no shadows, no gradients. Left-to-right flow in
> three stages.
>
> LEFT: a single worked example. A short protein ribbon segment in blue labelled
> "TP53 DNA-binding domain", with one residue highlighted showing an amino-acid
> substitution "R → H" (arginine to histidine) tagged "R175H".
>
> MIDDLE: an arrow leads to a vertical stack of six small clean icon cards, one per
> feature, each a simple pictogram with a short label: "Δ hydrophobicity",
> "Δ side-chain volume", "Δ charge", "fold-core location", "functional-switch
> residue", "hotspot / pathogenic" (this last card outlined in orange). A small
> bracket labels the six cards together as "θ — 6-dimensional biophysical vector".
>
> RIGHT: an arrow leads to an empty rounded panel labelled "shared variant feature
> space" (leave this area blank — a real scatter plot will be inserted here).
>
> Small grey footnote text at the bottom: "alternative feature spaces tested later:
> ESM, ESM + θ". Minimalist, tidy, colour-blind-safe, lots of white space, no charts,
> no numbers besides the labels given.

---

**Exact text to place/verify in Illustrator:**
- "TP53 DNA-binding domain", "R → H", "R175H"
- Six feature labels exactly as above; bracket "θ — 6-dimensional biophysical vector"
- "shared variant feature space" (the PCA scatter `fig1d_theta_pca.pdf` goes here)
- footnote "alternative feature spaces tested later: ESM, ESM + θ"

**Do NOT include:** ESM architecture detail, model performance, or any result. Keep θ
the star; ESM is only a small footnote so Fig 4e keeps its punch.
