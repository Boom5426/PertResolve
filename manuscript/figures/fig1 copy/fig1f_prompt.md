# Figure 1f — AI schematic prompt

**Panel role:** show the six generalization splits at a glance. Pure design, no data.

**Target tool:** general raster image model. **Aspect ratio:** ~3:1 wide strip, or a
2 × 3 grid of small cards.

**Palette:** white background; ink `#1A1A1A`; six muted card tints from a single
colour-blind-safe family; gene colours only if a card needs them
(TP53 `#0072B2`, KRAS `#D55E00`, GATA1 `#009E73`, JAK1 `#CC79A7`).

---

## Prompt (paste into the image model)

> A clean flat 2D vector schematic, Nature-Methods style, pure white background, thin
> strokes, no 3D, no shadows. Six small equal rounded cards in a tidy 2×3 grid, each
> with a tiny minimalist pictogram at top and two short text lines below (a name and
> a one-line description). A small shared icon motif shows a protein bar with some
> residues marked as "train" (grey) and some marked as "held-out" (orange), varied per
> card to convey what is held out.
>
> Card 1 "Random" — "random held-out variants".
> Card 2 "Positional" — "hold out one protein region (C-terminal half)".
> Card 3 "Mechanistic" — "hold out hotspot / functional-switch residues".
> Card 4 "Cross-gene" — "hold out a whole gene, train on the others".
> Card 5 "Low-depth" — "hold out shallow-sampling variants".
> Card 6 "Compatibility" — "match a prior evaluation protocol".
>
> Minimalist, evenly spaced, small neat sans-serif labels, muted colour-blind-safe
> palette, no charts, no numbers.

---

**Exact text to place/verify in Illustrator (card title / subtitle):**
- Random / "random held-out variants"
- Positional / "hold out one protein region (C-terminal half)"
- Mechanistic / "hold out hotspot or functional-switch residues"
- Cross-gene / "hold out a whole gene, train on the others"
- Low-depth / "hold out shallow-sampling variants"
- Compatibility / "match a prior evaluation protocol"

**Do NOT include:** any per-split result or PDS value. Names and one-liners only.
