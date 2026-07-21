# Figure 1e — AI schematic prompt

**Panel role:** show that evaluation is split into three axes, and that the central
requirement is allele *discrimination*, not merely getting the direction right.
No results — this defines the metrics conceptually.

**Target tool:** general raster image model. **Aspect ratio:** ~16:9 (three columns).

**Palette:** white background; ink `#1A1A1A`; direction column blue `#0072B2`;
discrimination column orange accent `#E69F00`; DE column green `#009E73`;
neutral grey `#7A7A7A` for candidates.

---

## Prompt (paste into the image model)

> A clean flat 2D vector scientific schematic, Nature-Methods style, pure white
> background, thin strokes, no 3D, no shadows. Three equal columns separated by thin
> vertical rules, and one horizontal banner across the bottom.
>
> LEFT column, header "Direction recovery": two overlaid line profiles pointing the
> same way (a predicted profile and an observed profile roughly parallel), small
> label "Pearson-δ" and a faint sub-label "δ-cosine". Sub-caption "does the
> prediction point the right way?".
>
> MIDDLE column, header "Allele discrimination": one predicted profile on the left
> with an arrow to a small ranked list of candidate variant profiles on the right;
> the correct/own variant is highlighted in orange as the nearest match, the others
> are grey. Label "PDS" and a small inequality "d(pred, own) < d(pred, others)".
> Sub-caption "does it identify its own allele among held-out variants?".
>
> RIGHT column, header "DE fidelity": three small stacked metric chips labelled
> "DE overlap", "DE-LFC rank corr.", "direction agreement", next to a tiny cartoon of
> a predicted vs measured gene-change list. Sub-caption "do the changed genes agree?".
>
> BOTTOM banner (full width, subtle grey fill): a single centred sentence
> "Central question: can models distinguish alleles, not merely recover the shared
> direction of a perturbed gene?".
>
> Minimalist, tidy, colour-blind-safe, small neat sans-serif labels, no real numbers.

---

**Exact text to place/verify in Illustrator:**
- Headers: "Direction recovery" / "Allele discrimination" / "DE fidelity"
- Metric names: "Pearson-δ", "δ-cosine", "PDS", "d(pred, own) < d(pred, others)",
  "DE overlap", "DE-LFC rank corr.", "direction agreement"
- Bottom banner sentence verbatim (above)

**Do NOT include:** any score value, any statement that models succeed or fail. This is
the metric-definition panel; the dissociation result belongs to Figure 2.
