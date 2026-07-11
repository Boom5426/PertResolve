# Figure 6g — AI schematic prompt (hero: power-aware workflow)

**Panel role:** the paper's final deliverable — an executable, decision-based protocol
for variant-level perturbation benchmarking. It must read as a real workflow with a
decision node and two branches, not a flat 4-step arrow. Positioned as a **triage tool,
not a universal cell-number calculator**. No results.

**Target tool:** general raster image model (or Illustrator by hand — text accuracy
matters, so overlay the exact labels rather than trusting generated text).
**Aspect ratio:** full width, ~4:1 (left→right with one downward branch).

**Palette:** white background; ink `#1A1A1A`; process boxes muted blue-grey; decision
diamond amber `#E8C468`; "yes" path green `#55966B`; "no" path muted red `#D98C8C`.

---

## Prompt (paste into the image model)

> A clean flat 2D vector process/decision diagram, Nature-Methods style, pure white
> background, thin strokes, no 3D, no shadows. A left-to-right pipeline of rounded boxes
> with one central decision diamond that branches into two labelled paths.
>
> BOX 1 "Pilot screen": 25-50 cells per perturbation, balanced WT/control, initial variant panel.
> BOX 2 "Estimate measurement properties": variant-to-WT effect size, within-variant split-half
> noise, pilot signal-to-noise.
> BOX 3 "Predict rankability": rankable probability, expected measurement regime.
> DECISION DIAMOND (amber) "Predicted to clear the measurement floor?" with two arrows:
>   - a green "Yes" arrow to BOX 4a "Proceed to full-depth benchmarking": size experiment to
>     target depth; evaluate direction, PDS and DE fidelity; compare models.
>   - a red "No" arrow to BOX 4b "Redesign or report as measurement-limited": increase depth
>     where plausible; enrich stronger phenotypes; change assay / stimulation / cell state;
>     improve variant assignment; do not interpret model ranking.
> Both paths converge into a final full-width box "Report benchmark by measurement regime":
> rankable / measurable perturbations vs below-floor / measurement-limited perturbations.
>
> A bottom principle strip reads: "benchmark models only where the ground truth is
> sufficiently resolved". Minimalist, evenly spaced boxes, thin connectors, small neat
> sans-serif labels, colour-blind-safe, no numbers besides "25-50 cells".

---

**Exact text to place/verify in Illustrator (this is text-heavy — overlay, don't trust AI text):**
- Box titles: "Pilot screen" / "Estimate measurement properties" / "Predict rankability" /
  "Proceed to full-depth benchmarking" / "Redesign or report as measurement-limited" /
  "Report benchmark by measurement regime"
- Decision: "Predicted to clear the measurement floor?"  ·  path labels "Yes" / "No"
- The bullet contents of each box, verbatim from the prompt above
- principle strip: "benchmark models only where the ground truth is sufficiently resolved"

**Do NOT include:** exact required cell numbers, a guaranteed-rankability claim, or a
universal sufficient depth. Wording must stay "predicted rankability / relative triage /
effect-size-conditioned design / expected measurement regime".
