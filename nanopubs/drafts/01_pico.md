# 01 — PICO Research Question (question-rooted chains, comparative)

> Use this draft instead of `01_quote.md` if your chain is question-rooted with a clear comparator (X vs Y). For descriptive/scoping question-rooted chains, use `01_pcc.md`. See `docs/chain-decision-tree.md`.
>
> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.
>
> **Chain shape:** question-rooted (PICO), chosen by Anne on 2026-10-09. The question was
> proposed by an external researcher; their credit is the unresolved token
> `{{RESEARCHER_CREDIT}}`, which blocks publication until they have agreed on how to
> be credited.
>
> Template (live, checked 2026-10-09): https://w3id.org/np/RA5e5XeXy_-aNK5giB7kBAEQslTLVydHeM4YYEzhmEE2w
> Fields: pico, label, description, type, populationDescription,
> interventionGroupDescription, comparatorGroupDescription, outcomeGroupDescription.

**Form heading:** *"PICO Research Question — Define a research question using the PICO framework (Population, Intervention, Comparator, Outcome)"*

## Field-by-field draft

<!-- field: pico -->
### Short ID used as URI suffix (text input, required)

Slug becomes part of the nanopub URI. Use kebab-case.

```
basal-signals-lead-glacier-acceleration
```

<!-- field: label -->
### Label for the research question (text input, required)

10-200 characters. Length-bounded.

```
Do basal water pressure and uplift give an earlier and better signal of glacier acceleration than ice speed?
```

<!-- field: description -->
### Description of the research question (textarea, required)

One coherent sentence/paragraph that names P, I, C, O inline.

```
For grounded outlet glaciers whose flow responds to surface meltwater or rain reaching the bed, do observations of the ice–bed boundary (basal water pressure and bed uplift), added to weather and melt observations and the glacier-speed history, predict glacier acceleration earlier and more skilfully than the speed history alone or the speed history with weather and melt? The key question is whether boundary observations change before acceleration, providing genuine lead time, rather than at the same moment as speed. Question proposed by {{RESEARCHER_CREDIT}}.
```

<!-- field: type -->
### Question Type (dropdown, required)

- [ ] causation research question - (Does factor X cause outcome Y?)
- [ ] descriptive research question - (What are the characteristics of X?)
- [ ] effectiveness research question - (Does approach X work better than Y?)
- [ ] experience research question - (How do people experience phenomenon X?)
- [x] prediction research question - (What outcomes can we expect from X?)

<!-- field: populationDescription -->
### Description of the population (textarea, required)

Who/what is being studied. Discipline-level concept — not implementation. See `docs/pico-study-outcome-levels.md`.

```
Grounded outlet glaciers, at the contact between grounded ice and its bed, during acceleration events driven by surface meltwater or rain reaching the bed.
```

<!-- field: interventionGroupDescription -->
### Description of the intervention group (textarea, required)

The intervention or exposure being examined. Discipline-level concept.

```
Observations of the ice–bed boundary — basal water pressure and bed uplift — used together with weather and melt observations and the glacier-speed history.
```

<!-- field: comparatorGroupDescription -->
### Description of the comparator group (textarea, required)

The comparison or control condition. Discipline-level concept.

```
The glacier-speed history alone, and the glacier-speed history together with weather and melt observations, without observations of the ice–bed boundary.
```

<!-- field: outcomeGroupDescription -->
### Description of the outcome group (textarea, required)

What outcomes are being measured. The kind of measurement, not the value.

```
Glacier acceleration: skill in predicting acceleration on events not used to build the prediction, and the lead time between the onset of change at the ice–bed boundary and the onset of acceleration.
```

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 01.
