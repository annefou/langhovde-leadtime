# 02 — AIDA Sentence

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Form heading:** *"AIDA Sentence — Make structured scientific claims following the AIDA model"*

> Template fields (live 2026-10-09): aida, topic (optional, left empty), project (filled by the wizard from the PICO URI), dataset (optional, skipped), publication (optional).
>
> The AIDA states the hypothesis the PICO asks about, as a declarative answer (`docs/chain-decision-tree.md` § Question-rooted). It is atomic: basal **pressure** only. Uplift would be a second AIDA. The Outcome reports whether this event bears it out.

## Field-by-field draft

<!-- field: aida -->
### AIDA sentence (text input, required)

Atomic, Independent, Declarative, Absolute. One empirical finding. Must end with a full stop.

> _If your draft AIDA contains "and" linking two distinct findings, split into two AIDA nanopubs._

```
Basal water pressure at the bed of a grounded glacier rises before the glacier accelerates during a meltwater-driven speed-up.
```

<!-- field: topic -->
### Select related topics/tags (search/select, optional)

Predefined topic vocabulary — list the labels you intend to pick from the dropdown.

```

```

<!-- field: project -->
### Relates to this nanopublication (search/select, required)

URI of the nanopub the AIDA derives from.

- For paper-rooted chains: the Quote-with-comment URI (from step 01).
- For question-rooted chains: the PICO or PCC URI (from step 01).

Pull the URI from `nanopubs/PUBLISHED.md`.

```

```

<!-- field: dataset -->
### Supported by datasets (text input, optional)

DOIs/URLs of datasets that ground the AIDA claim.

- _(skip — optional; the datasets test the claim rather than support it)_

<!-- field: publication -->
### Supported by other publications (text input, optional)

DOIs/URLs of publications that support the AIDA claim — e.g. peer-reviewed methods papers, or the original paper if not already cited via the Quote.

Sugiyama et al. (2026): meltwater reaching the bed raised basal pressure and accelerated grounded ice at Langhovde.

- https://doi.org/10.1038/s41467-026-72724-x

> **Known platform bug (2026-04-26):** if both *Supported by datasets* AND *Supported by other publications* are populated and publishing fails, fall back to publishing this AIDA via Nanodash. The URI namespace becomes `https://w3id.org/np/...` (still valid and citable).

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 02.
