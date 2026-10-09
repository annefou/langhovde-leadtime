# 05 — FORRT Replication Outcome

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Verify the actual numerical results first** by reading `results/` and `notebooks/03_analysis.py`. Don't quote numbers from memory. See `docs/verify-before-drafting.md`.

> Template fields (live 2026-10-09): outcome, label, study (filled by the wizard), repo, date, validationStatus, conclusion, evidence, confidenceLevel, limitations (optional).
>
> `repo` is the v1.0.0 version DOI (10.5281/zenodo.23266427, from `CITATION.cff`).
>
> Confidence: the Science Live vocabulary is phrased as agreement with an original. "low — limited evidence" is the closest fit to a one-event, inconclusive result. Anne to confirm.
>
> Numbers from `results/`: `lead_labels.csv`, `step6_*`, `step7*`, `step8_*`, `step9_*`, `kinematic_*`, `nested_models*`, `robust_delay_check.csv`.

## Field-by-field draft

<!-- field: outcome -->
### Short URI suffix for outcome ID (text input, required)

Slug. Use kebab-case.

```
langhovde-basal-lead-time-2022-outcome
```

<!-- field: label -->
### Plain-text label for the outcome (text input, required)

Descriptive title.

```
Basal pressure vs acceleration onsets at Langhovde Glacier, January 2022: inconclusive (one event)
```

<!-- field: study -->
### Choose study (search/select, required)

URI of the Replication Study published in step 04. Pull from `nanopubs/PUBLISHED.md`.

```

```

<!-- field: repo -->
### Repository URL (text input, required)

Use the Zenodo **version DOI** URL for the release the results came from — not a
bare branch URL, and not the concept DOI.

> **Why not the bare repo URL.** `https://github.com/ORG/REPO` names a *moving
> branch*. This Outcome asserts "this code produced this number", in a signed,
> immutable record. A branch URL means that assertion points at whatever `main`
> happens to be years from now — code that may never have produced the number
> above. A concept DOI has the same flaw: it resolves to the latest version.
> The version DOI pins the exact release. `docs/chain-decision-tree.md` § Anchor
> ranks the options: SWHID > Zenodo DOI > repo URL > Wayback.
>
> Both DOIs and the SWHID are in `CITATION.cff` under `identifiers:`, recorded
> automatically at release by `.github/workflows/release-identifiers.yml`. Take
> the one described as *"Version DOI"*.

```
https://doi.org/10.5281/zenodo.23266427
```

<!-- field: date -->
### Choose completion date (text input, required)

```
2026-10-09
```

<!-- field: validationStatus -->
### Choose validation status (dropdown, required)


This dropdown maps to the CiTO intention in step 06: Validated → `confirms`, PartiallySupported → `qualifies`, Contradicted → `disputes`.

- [ ] contradicted
- [x] inconclusive
- [ ] not tested
- [ ] partially supported
- [ ] validated

<!-- field: confidenceLevel -->
### Choose confidence level (dropdown, required)

_Vocabulary not yet captured._

```

```

- [ ] high - Strong evidence, mostly agrees with original
- [x] low - Limited evidence, significant disagreement
- [ ] moderate - Adequate evidence, partial agreement
- [ ] very high - Extensive evidence, high agreement with original
- [ ] very low - Minimal evidence, major disagreement

<!-- field: conclusion -->
### Describe the overall conclusion about the original claim (textarea, required)

Substantive interpretation. Headline comparison: replication's number vs the paper's number, sign + significance.

```
This record cannot tell whether basal water pressure rises before grounded-ice acceleration. The claim is therefore neither supported nor contradicted.
What is clear:
- In the one instrumented event (2–6 January 2022), rain and melt came first: rain began at 01:50 UTC on 2 January, and the basal pressure rise followed about 19 h later in every processing.
What is not clear:
- The order of pressure and acceleration depends on the processing.
  - Centred smoothing of the published 15 min GNSS puts the speed rise about 8 h before the pressure rise.
  - The same smoothing of the reprocessed 1 s GNSS puts pressure about 6 h before speed.
  - Every past-only (causal) estimate is indeterminate, with spreads of more than ±10 h.
Why sharper data did not resolve it:
- The grounded and floating trunk stations share a ~12 h variability in speed that is larger than the event's speed-up. It is mostly not GNSS multipath and not a simple tidal response.
- The borehole tilt sensor was still settling.
The proposed held-out comparison of three nested models needs several events. In a single-event hold-out, adding basal pressure and uplift did not improve on weather and speed history.
```

<!-- field: evidence -->
### Describe the evidence that supports your conclusion (textarea, required)

Numerical results, test statistics, model coefficients. Read directly from `results/`.

```
Basal pressure and forcing:
- BH2201 basal water level peaked at 51.0 m a.s.l. (97% of flotation) at 23:07 UTC on 3 January 2022; 50.95 m after the measured air-pressure correction.
- Rain led the pressure rise by a median of 18.6–18.9 h (causal) and 20.8–23.4 h (centred), in 30–32 of 32 rule × bandwidth × reference combinations.
Pressure to GNSS1 (grounded) speed, onset leads (positive = pressure first):
- 15 min deposit:
  - centred median −8.0 h ("level lags", 20 of 30 combinations);
  - robust causal median +0.6 h, 13 earlier and 13 later ("indeterminate").
- 1 s kinematic GNSS (100% fixed; 5.7 mm and 3.7 mm horizontal SD against the published 15 min positions):
  - centred median +5.75 h ("level leads", 16 of 32 combinations defined);
  - robust causal, corrected for detection delay, median +3.6 h, range −7.6 to +36.1 h, 13 of 32 defined ("indeterminate");
  - after sidereal filtering and tidal correction, still "indeterminate" (median +4.9 to +6.8 h).
- Detection delay of the causal speed estimators on synthetic onsets: 3.0–8.9 h (15 min data); 3.4 h (1 s data, h = 3 h).
Sub-daily variability at GNSS1:
- 12–13 h period, SD 0.047 m/d at h = 1 h, against an event speed-up of about 0.03 m/d (12.5%).
- Sidereal filtering removed 14% (GNSS1) and 32% (GNSS2) of the 10–15 h variance.
- A tidal model fitted before the event explained R² = 0.001 and 0.07 afterwards.
- GNSS1 and GNSS2 band-passed (6–18 h) residuals correlate at r = 0.97 (east), 0.72 (north) and 0.84 (up). GNSS1 and receiver LG1 correlate at 0.17, 0.23 and 0.00.
Bed tilt and nested models:
- Bed tilt in BH2201 drifted about 0.9 mrad/h, with 13 steps above 5 mrad (largest 77 mrad).
- Nested models, one-event hold-out at GNSS1 (h = 3 h): skill of the pressure-and-uplift model relative to weather-and-speed was −0.24, −0.05 and −0.01 at 1, 3 and 6 h horizons.
- Rolling-origin forecasts began only after the event onset.
```

<!-- field: limitations -->
### Describe what limits the conclusions of the study (textarea, optional)

Honest caveats. If the result is partial or contradicted, say so plainly. Don't overclaim.

```
- One event, at one site, in one season. No conclusion about lead time in general, or at other glaciers, can be drawn.
- The proposed predictive test was not run. The lead labels describe this event only.
- Basal pressure comes from one borehole, about 190 m from the grounded GNSS station.
- Onset times depend on the smoothing and the onset rule. All combinations are reported, not a chosen one.
- The origin of the ~12 h speed variability is not established: real trunk motion is likely, but differential troposphere cannot be excluded. What the two extra GNSS receivers stand on is not documented in the archive.
- The bed-tilt sensor was settling (drift and steps), so its onsets were withdrawn. The tilt-rate change during the pressure rise and the stick-slip-like steps after the peak are post hoc observations, not tested.
- Several steps were amended after results were seen (dated in ANALYSIS_PLAN.md).
```

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 05.
