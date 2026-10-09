# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # 05 — NIPR ADS downloads (Step 6 inputs)
#
# Langhovde 2021/22 field data from the NIPR Arctic/Antarctic Data archive System
# (ADS), all by Sugiyama, S., M. Minowa, K. Kondo, S. Aoki (2022), CC BY 4.0. The
# embargo ended on 2024-05-05 (`open_date` in the ADS catalogue; the "Embargo" label
# was not updated). Access and licence were checked by Anne on 2026-10-09.
#
# | ADS ID | Content |
# |---|---|
# | A20220506-001 | On-glacier AWS at LG4, 10 min, incl. air pressure |
# | A20220506-002 | Borehole 3 (field name) water pressure + accelerometer |
# | A20220506-004 | Raw GNSS, 1 s: LG1, LG4, LG5, LGFX (reference), LGLK |
#
# Usage: `python 05_nipr_download.py [ID ...] [--prefix LG5 ...]`. Files are skipped
# when already present with the listed size; a manifest with sizes is written.

# %%
import json
import sys
from pathlib import Path

import requests

API = "https://ads.nipr.ac.jp/api/v1/metadata"
OUT = Path("../data/raw/nipr")
DATASETS = ["A20220506-001", "A20220506-002", "A20220506-004"]


def listing(ds: str, path: str | None = None) -> list[dict]:
    params = {"sub_directory": "true"}
    if path:
        params["path"] = path
    r = requests.get(f"{API}/{ds}/0.00/directory/DATA", params=params, timeout=120)
    r.raise_for_status()
    return r.json()


def files(ds: str, path: str | None = None) -> list[dict]:
    """Every file under `path`, descending into folders whose children are not listed."""
    out = []
    for e in listing(ds, path):
        stack = [e]
        while stack:
            n = stack.pop()
            if not n["directory"]:
                out.append(n)
            elif n["children"]:
                stack.extend(n["children"])
            elif int(n["number_of_files"]) > 0:
                out.extend(files(ds, n["path"]))
    return out


def fetch(ds: str, f: dict, kind: str = "DATA") -> Path:
    dest = OUT / ds / f["path"]
    if dest.exists() and dest.stat().st_size == int(f["size"]):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    url = f"{API}/{ds}/0.00/data/{kind}"
    with requests.get(url, params={"path": f["path"]}, stream=True, timeout=600) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 22):
                fh.write(chunk)
    if tmp.stat().st_size != int(f["size"]):
        raise RuntimeError(f"{ds}/{f['path']}: size {tmp.stat().st_size} != listed {f['size']}")
    tmp.rename(dest)
    return dest


# %%
if __name__ == "__main__":
    args = sys.argv[1:]
    prefixes = args[args.index("--prefix") + 1:] if "--prefix" in args else []
    ids = [a for a in (args[: args.index("--prefix")] if "--prefix" in args else args)] or DATASETS
    for ds in ids:
        fl = files(ds)
        if prefixes:
            fl = [f for f in fl if any(f["path"].startswith(p) for p in prefixes)]
        total = sum(int(f["size"]) for f in fl)
        print(f"{ds}: {len(fl)} files, {total / 1e9:.2f} GB", flush=True)
        for i, f in enumerate(fl, 1):
            fetch(ds, f)
            if i % 10 == 0 or i == len(fl):
                print(f"  {i}/{len(fl)}", flush=True)
        manifest = OUT / ds / "manifest.json"
        old = json.loads(manifest.read_text()) if manifest.exists() else {}
        old.update({f["path"]: int(f["size"]) for f in fl})
        manifest.write_text(json.dumps(old, indent=1, sort_keys=True))
        pdf = OUT / ds / f"{ds}_v000.pdf"
        if not pdf.exists():
            r = requests.get(f"{API}/{ds}/0.00/data/PDF", params={"path": pdf.name}, timeout=120)
            if r.ok:
                pdf.write_bytes(r.content)
