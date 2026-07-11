# Figure 5e — AI schematic prompt

**Panel role:** make the abstract "benchmark validity" concrete: the SAME true model
order is scrambled by a low-resolution (near-floor) benchmark but recovered by a
higher-resolution one. Illustrative concept (uses the known true order from 5d), not a
data plot. No real numbers.

**Target tool:** general raster image model. **Aspect ratio:** ~16:9 (two side-by-side
mini-leaderboards).

**Palette:** white background; ink `#1A1A1A`; true-order accent blue `#5185C0`; observed
grey `#9AA0A6`; mismatch/crossing lines muted red `#D98C8C`; recovered lines green
`#55966B`.

---

## Prompt (paste into the image model)

> A clean flat 2D vector schematic, Nature-Methods style, pure white background, thin
> strokes, no 3D, no shadows. Two side-by-side mini-panels separated by a thin vertical
> rule.
>
> LEFT mini-panel, header "Ceiling near chance (low resolution)": on its left a vertical
> "true quality" list of five predictors P1..P5 (best at top, blue); on its right an
> "observed benchmark score" list of the same five in a SCRAMBLED order with overlapping
> error bars; thin muted-red lines connect each predictor between the two lists, crossing
> each other to show the ranking is scrambled.
>
> RIGHT mini-panel, header "Ceiling above the floor (high resolution)": the same true
> quality list on the left and the observed list on the right, but now in the SAME order
> with separated error bars; thin green connecting lines run straight across (no
> crossings), showing the ranking is recovered.
>
> A single caption strip along the bottom: "the truly better model need not rank higher
> when the measurement ceiling is near chance". Minimalist, tidy, small neat sans-serif
> labels, colour-blind-safe, no real numbers.

---

**Exact text to place/verify in Illustrator:**
- Headers: "Ceiling near chance (low resolution)" / "Ceiling above the floor (high resolution)"
- "true quality", "observed benchmark score", predictor labels P1–P5
- caption strip (verbatim above)

**Do NOT include:** real PDS values, real dataset names, oracle numbers. This is an
illustrative concept panel; the quantitative curve is 5f.
