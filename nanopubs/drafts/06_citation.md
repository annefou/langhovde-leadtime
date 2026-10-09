# 06 — CiTO Citation

> Run the pre-flight checklist in `docs/forrt-form-fields.md` § Pre-flight checklist before drafting.

**Description:** *"Declare citations between papers or other works, using Citation Typing Ontology"*

> Template fields (live 2026-10-09): work (filled by the wizard with the Outcome URI), st02 citations [{cites, cited}] ≥ 1.
>
> The Outcome is `Inconclusive`, which has no canonical CiTO intention, so citation 1 uses the neutral `discusses` for Sugiyama et al. (2026). The other citations:
> - `extends`: the prior chain's Outcome (reproduction of the same paper; this work uses its code and data and adds the lead-time analysis);
> - `credits`: the prior repository, whose v1.0.2 GNSS code is reused;
> - `citesAsDataSource`: the deposit and the three NIPR ADS records;
> - `usesMethodIn`: Agnew & Larson (GPS repeat period for the sidereal filter);
> - `credits` with `{{RESEARCHER_URI}}`: the researcher who proposed the question. This blocks publication until resolved.
>
> All DOIs were checked with `resolve_doi` on 2026-10-09.

## Field-by-field draft

<!-- field: work -->
### Identifier for the citing creative work (text input, required)

URI of the Outcome published in step 05. Pull from `nanopubs/PUBLISHED.md`.

```

```

### List citations (repeatable group, required ≥1)

#### Citation 1 — back to the original paper

##### Citation Type (dropdown)

Choose based on the Outcome's validation status:

- Validated → `confirms`
- PartiallySupported → `qualifies`
- Contradicted → `disputes`

For question-rooted chains where there is no original paper to confirm/dispute, use `usesMethodIn` or `citesAsAuthority` for the methodology paper(s).

Write the chosen type in the block below (a vocabulary label such as `cites as authority`, or `citesAsAuthority`). `build-chain-draft` uses it as written; leave the block empty to have the type derived from the Outcome's validation status, which is right for paper-rooted chains only.

> **Note:** `replicates` is NOT in the Science Live dropdown (despite existing in upstream CiTO). When citing a notebook/tutorial that was directly reused, use **`credits`** instead.

```
discusses
```

##### DOI or other URL of the cited work (text input)

```
https://doi.org/10.1038/s41467-026-72724-x
```

#### Additional citations (optional)

If the Outcome cites methods papers, related replications, or upstream tools, add them here.

One line per further citation, in this exact form (each becomes a pre-filled row):

- Type: extends → URL: https://w3id.org/sciencelive/np/RAQn6_6v0w8OARZmzh0BBV8maDh_dqNHBlYmgwzNMwGyo
- Type: credits → URL: https://doi.org/10.5281/zenodo.23257925
- Type: citesAsDataSource → URL: https://doi.org/10.17632/8wvtxg53ry.1
- Type: citesAsDataSource → URL: https://ads.nipr.ac.jp/dataset/A20220506-004
- Type: citesAsDataSource → URL: https://ads.nipr.ac.jp/dataset/A20220506-001
- Type: citesAsDataSource → URL: https://ads.nipr.ac.jp/dataset/A20220506-002
- Type: usesMethodIn → URL: https://doi.org/10.1007/s10291-006-0038-4
- Type: credits → URL: {{RESEARCHER_URI}}

## Publication note

After publishing, paste the resulting URI into `nanopubs/PUBLISHED.md` step 06.

This completes the six-step FORRT chain. Optional next layers:

- **Research Software** (`drafts/07_research_software.md`) — if the repo *produces* a reusable software artefact.
- **Research Synthesis** (`drafts/08_synthesis.md`) — if this chain is one of several testing facets of a shared property.
