# Figure 1g — AI schematic prompt

**Panel role:** tie the definitions together into one end-to-end benchmark loop, so the
reader enters Figure 2 already knowing how a prediction is made and scored. No results.

**Target tool:** general raster image model. **Aspect ratio:** ~4:1 wide horizontal flow.

**Palette:** white background; ink `#1A1A1A`; steps in a single muted blue-grey family;
orange accent `#E69F00` only on the "held-out" element; gene colours optional.

---

## Prompt (paste into the image model)

> A clean flat 2D vector process diagram, Nature-Methods style, pure white background,
> thin strokes, no 3D, no shadows. A single left-to-right pipeline of four boxes
> connected by thin arrows.
>
> Box 1 "Variant + features": a protein icon with a small θ feature chip attached.
> Box 2 "Train on observed variants": a few example variants marked as training
> (grey), feeding a small generic "feature → response" predictor block.
> Box 3 "Predict held-out response": the trained predictor outputs a pseudobulk
> response profile for a NEW variant marked "held-out" in orange.
> Box 4 "Evaluate": the predicted profile is compared to the measured one, with three
> small tags "direction", "discrimination", "DE fidelity".
>
> Keep the four boxes evenly sized and spaced, thin connecting arrows, small neat
> sans-serif labels, minimalist and colour-blind-safe, no charts, no numbers.

---

**Exact text to place/verify in Illustrator:**
- Box titles: "Variant + features" / "Train on observed variants" /
  "Predict held-out response" / "Evaluate"
- Inside labels: "θ features", "feature → response predictor",
  "pseudobulk perturbation profile", "held-out variant"
- Evaluation tags: "direction", "discrimination", "DE fidelity"

**Do NOT include:** any outcome, score, ranking, or claim that the pipeline succeeds or
fails. This is a workflow schematic; results start in Figure 2.
