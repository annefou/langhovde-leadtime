# 04 — FORRT Replication Study

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Verify code first:** read the actual reproduction script in `notebooks/03_analysis.py` before writing the methodology field. See `docs/verify-before-drafting.md`.

> Template fields (live 2026-10-09): study, label, type, claim (filled by the wizard), scope, methodology, deviation (optional), keyword (optional, repeatable, Wikidata), discipline (optional, Wikidata).
>
> Study type: no Science Live option fits a first test of a new question exactly. "Replication Study (different methodology or conditions)" is the closest; the deviation field says so.
>
> Methodology checked against `notebooks/03_analysis.py`, `06_kinematic.py`, `07_step6_analysis.py`, `08_step7.py`, `09_step8.py` and `10_step9.py`.

## Field-by-field draft

<!-- field: study -->
### Short URI suffix for study ID (text input, required)

Slug. Use kebab-case.

```
langhovde-basal-lead-time-2022
```

<!-- field: label -->
### Label/name of replication study (text input, required)

Human-readable title.

```
Lead time of basal water pressure before glacier acceleration, Langhovde Glacier, January 2022 (exploratory, one event)
```

<!-- field: type -->
### Choose the study type (dropdown, required)

- [x] Replication Study - replication with different methodology or conditions
- [ ] Reproduction/Replication Study - study that is both, reproduction and replication
- [ ] Reproduction Study - direct reproduction: same methodology, same tools

<!-- field: claim -->
### Choose FORRT claim (search/select, required)

URI of the Claim published in step 03. Pull from `nanopubs/PUBLISHED.md`.

```

```

<!-- field: scope -->
### Describe what part of the claim is reproduced/replicated. (textarea, required)

The **scope** of the claim being tested. Which aspect, what's in/out of scope. NOT methodology. NOT results. See `docs/pico-study-outcome-levels.md`.

```
Whether basal water pressure (and bed uplift) changes before the acceleration of grounded ice, tested on the one meltwater- and rain-driven speed-up at Langhovde Glacier, East Antarctica, that was observed with basal pressure, GNSS and weather together: 2–6 January 2022 (summer 2021/22).
- In scope: the order and separation in time of the onsets of rain and melt, basal water pressure, bed uplift, bed tilt, and the speed and acceleration of a grounded and a floating GNSS station; the lead–lag structure between these series; and a descriptive version of the proposed comparison of three nested models (speed history; plus weather and melt; plus basal pressure and uplift).
- Out of scope: the proposed held-out prediction across events. It needs several events with basal pressure, and this record has one, so it could not be run. Also out of scope: other glaciers or seasons, and the physical mechanism linking pressure and sliding.
```

<!-- field: methodology -->
### Describe how the claim is reproduced/replicated. (textarea, required)

The **method** in plain prose. Read `notebooks/03_analysis.py` and any config files first. NOT exact numerical results.

```
Exploratory analysis of one event, pre-registered step by step before each computation, with dated amendments.
Data:
- the Sugiyama et al. (2026) deposit (Mendeley Data, doi:10.17632/8wvtxg53ry.1): 15 min GNSS positions, 1 min borehole pressure (BH2201), 10 min weather;
- raw field data at the NIPR Arctic Data archive System (CC BY 4.0): 1 s GNSS from five receivers (A20220506-004), on-glacier air pressure (A20220506-001), and a 100 Hz accelerometer in borehole BH2201 (A20220506-002).
The deposit is parsed with the code of the earlier reproduction (doi:10.5281/zenodo.23257925).
Processing:
- The raw GNSS is processed kinematically with RTKLIB 2.5.1 against a reference receiver on rock: forward filter only, so that no future data are used; GPS L1+L2; integer ambiguities fixed; checked against the published 15 min positions.
- Basal pressure is corrected with the measured air pressure.
- Every series is smoothed three ways: centred, causal (past data only) and robust causal Gaussian local-linear regression, at several bandwidths. Speed comes from the smoothed positions; acceleration from the speed.
Analysis:
- Onsets are defined by declared rules (threshold exceedance at 2, 3 and 5 standard deviations of a pre-event day, persisting for 1 h; and a hinge change-point fit).
- Leads are computed for every rule, bandwidth and reference window. A pre-registered majority rule with a 1 h tolerance labels each pair: basal leads, simultaneous, basal lags or indeterminate.
- The detection delay of each causal speed estimator is measured on synthetic onsets with the observed noise and subtracted.
- Further steps: lagged cross-correlation; lagged regressions for the three nested models, with rolling-origin and event hold-out evaluation; tests of the sub-daily speed variability (sidereal filtering for GNSS multipath, a tidal regression fitted before the event, a common-mode comparison with two further GNSS receivers); and bed-tilt series from the borehole accelerometer.
```

<!-- field: deviation -->
### Describe any deviations from original methodology. (textarea, optional)

What's different from the original method. Verify against the actual code, don't guess.

```
There is no original methodology to reproduce: the question, and the three-model held-out test, were proposed by [name to be added with the researcher's agreement], who had not run it.
Deviations from that proposed test:
- one event only, so no held-out prediction across events (descriptive onset timing, lead–lag and a one-event hold-out instead);
- speed and uplift come from GNSS positions smoothed with the processing of the earlier reproduction, extended with past-only (causal) and robust variants;
- post hoc amendments, each dated and recorded:
  - the same smoother for all series;
  - a robust causal fit, after single-epoch GNSS outliers dominated the causal speed;
  - one base position and the base receiver's broadcast ephemeris in the kinematic processing;
  - bed-tilt onset labels withdrawn because the sensor was still settling.
```

<!-- field: keyword -->
### Search keywords (Wikidata) (search/select, optional)

Provide labels (not QIDs) — the Wikidata search picks up labels.

Checked with wikidata_lookup on 2026-10-09: Q6486120, Q3962812, Q179435.

- Langhovde Glacier
- Basal sliding
- global navigation satellite system

<!-- field: discipline -->
### Search discipline (Wikidata) (search/select, optional)

Provide labels.

Q52120; Wikidata types it as a field of study, not an academic discipline.

- glaciology

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 04.
