# langhovde-leadtime

> **Does the ice–bed boundary give an earlier signal of glacier acceleration than ice speed?**
> An exploratory test at Langhovde Glacier, East Antarctica (one event, 2–6 January 2022).
> Follow-up to the [FORRT reproduction](https://w3id.org/sciencelive/np/RAQn6_6v0w8OARZmzh0BBV8maDh_dqNHBlYmgwzNMwGyo) of Sugiyama et al. (2026), [doi:10.1038/s41467-026-72724-x](https://doi.org/10.1038/s41467-026-72724-x).

**Question.** Do basal water pressure and bed uplift, added to weather and ice-speed
history, give an *earlier and better* signal of glacier acceleration than the speed
record itself? The proposed test is held-out prediction from three nested models (speed
history; plus weather and melt; plus basal pressure and uplift). The key result is
whether boundary observations change *before* acceleration (lead time), rather than at
the same moment.

**Answer from this record: inconclusive.**
- Only one acceleration event was observed with basal pressure, so the held-out test
  cannot be run.
- For that event, rain led the basal pressure rise by about 19 h.
- The order of pressure and acceleration flips with the processing. Every past-only
  estimate is indeterminate, by ± 10 h or more.
- Reprocessing the raw 1 s GNSS (100% fixed solutions) moves the limit from GNSS noise
  to a ~12 h trunk motion that is larger than the event's speed-up.
- The borehole tilt sensor was still settling.

Details: [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md) (pre-registered steps, dated amendments),
`results/`, `figures/main_result.png`, and the drafted FORRT chain in `nanopubs/drafts/`.

**Credit.**
- The research question was proposed by an external researcher; their name will be
  added with their agreement.
- The code and analysis are by Anne Fouilloux, developed with an AI coding agent
  (Claude, Anthropic).

**Data.**
- Sugiyama et al. (2026) deposit, Mendeley Data
  [doi:10.17632/8wvtxg53ry.1](https://doi.org/10.17632/8wvtxg53ry.1).
- Raw field data from the NIPR Arctic Data archive System (Sugiyama, Minowa, Kondo &
  Aoki 2022, CC BY 4.0), all released 2024-05-05:
  - GNSS, 1 s: [A20220506-004](https://ads.nipr.ac.jp/dataset/A20220506-004);
  - weather station: [A20220506-001](https://ads.nipr.ac.jp/dataset/A20220506-001);
  - borehole pressure and accelerometer: [A20220506-002](https://ads.nipr.ac.jp/dataset/A20220506-002).
- Code reused from [doi:10.5281/zenodo.23257925](https://doi.org/10.5281/zenodo.23257925).

**Pipeline.**
- `pixi run snakemake --cores 1` runs Steps 0–4 (deposit only).
- Steps 6–9 need about 7 GB of NIPR downloads and RTKLIB (built by
  `scripts/build_rtklib.sh`). They are notebooks `05`–`10`, run in order.

## Quick start

```bash
git clone https://github.com/annefou/langhovde-leadtime.git
cd langhovde-leadtime
pixi install
pixi run snakemake --cores 1
```

Or with Docker:

```bash
docker run --rm ghcr.io/annefou/langhovde-leadtime:latest
```

## Structure

- `paper/` — the source paper PDF (drop yours in there).
- `notebooks/` — jupytext `.py` notebooks that drive the pipeline.
- `data/` — downloaded by `notebooks/01_data_download.py`, never committed.
- `nanopubs/` — drafts of the FORRT chain field-by-field, plus the published-URI registry.
- `docs/` — operating manuals (FORRT form fields, chain decision tree, claim-type vocabulary).
- `figures/` — curated figures used in the Jupyter Book.

## Nanopublication chain

The published chain is listed in [`nanopubs/PUBLISHED.md`](nanopubs/PUBLISHED.md). Each step links to its viewer URL on the Science Live platform.

## Citation

If you use this work, please cite both:

- This software: [`CITATION.cff`](CITATION.cff) → DOI [{{ZENODO_DOI}}]({{ZENODO_DOI}}).
- The original paper: [10.1038/s41467-026-72724-x](https://doi.org/10.1038/s41467-026-72724-x).
