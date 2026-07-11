# Figure 2a — AI schematic prompt

**Panel role:** small orientation schematic reminding the reader what Fig 2 compares: a
held-out variant's predicted response is judged two ways — direction recovery (is the
predicted profile pointing the right way?) and allele discrimination (does it identify the
correct held-out variant?). No result numbers, no chance line, no measurement-floor content.

**Target tool:** general raster image model. **Aspect ratio:** ~3:4 or ~1:1 (compact, left→right→branch).

**Palette:** white background; ink `#1A1A1A`; TP53 blue `#5185C0`; direction branch blue;
discrimination branch orange `#E99D4E`; candidate profiles grey `#9AA0A6`.

---

## Prompt (paste into the image model)

> A clean, flat 2D vector scientific schematic, Nature-Methods style, pure white
> background, thin 0.5 pt strokes, no 3D, no shadows, no gradients. Left-to-right with a
> two-way branch on the right.
>
> LEFT "Training": a short list of observed variants "TP53 R248Q", "TP53 R273C",
> "TP53 G245S". Below/beside, a held-out variant "TP53 R175H" tagged "held-out".
>
> MIDDLE: an arrow from the held-out variant's features through a small box
> "feature → response predictor" to a predicted profile labelled "predicted response".
>
> RIGHT, TWO branches from the predicted profile:
>   - UPPER branch "Direction recovery": the predicted profile overlaid with the measured
>     profile of the same variant, roughly parallel, labelled "Pearson-delta". Small
>     question "Is the predicted direction correct?".
>   - LOWER branch "Allele discrimination": the predicted profile compared to a small set of
>     candidate measured profiles (grey) for the other held-out variants, one highlighted as
>     the nearest match, labelled "PDS". Small question "Does it identify the correct allele?".
>
> Minimalist, tidy, small neat sans-serif labels, colour-blind-safe, no numbers, no charts.

---

**Exact text to place/verify in Illustrator:**
- "Training", "TP53 R248Q / R273C / G245S", "TP53 R175H", "held-out"
- "feature → response predictor", "predicted response"
- "Direction recovery", "Pearson-delta", "Is the predicted direction correct?"
- "Allele discrimination", "PDS", "Does it identify the correct allele?"

**Do NOT include:** any PDS/Pearson value, chance line, per-gene result, or model name from
the benchmark. This is an orientation schematic only; results are the data panels 2b-2g.
