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
# # 01 — Data download (Step 0)
#
# Two inputs, both verified on download (`ANALYSIS_PLAN.md` § 0):
#
# - **Prior code**, `annefou/langhovde-meltwater-replication` v1.0.2, Zenodo
#   [doi:10.5281/zenodo.23257925](https://doi.org/10.5281/zenodo.23257925). Its
#   `notebooks/gnss.py` is imported, unmodified, by `03_analysis.py`. Checked against
#   the MD5 that Zenodo publishes for the archive.
# - **Field data**, Sugiyama et al. (2026), Mendeley Data
#   [doi:10.17632/8wvtxg53ry.1](https://doi.org/10.17632/8wvtxg53ry.1), CC BY 4.0. Each
#   extracted file is checked against the SHA-256 published by the Mendeley API.
#
# No credentials are needed.

# %%
import hashlib
import json
import zipfile
from pathlib import Path

import requests

RAW_DIR = Path("../data/raw")
RAW_DIR.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# ## Prior code (Zenodo v1.0.2)

# %%
ZENODO_RECORD = "23257925"
ZENODO_API = f"https://zenodo.org/api/records/{ZENODO_RECORD}"
PRIOR_DIR = RAW_DIR / "prior_v1.0.2"


def file_hash(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists():
        return
    with requests.get(url, stream=True, timeout=600) as r:
        r.raise_for_status()
        tmp = path.with_suffix(".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
        tmp.rename(path)


record = requests.get(ZENODO_API, timeout=60).json()
assert record["metadata"]["version"] == "v1.0.2", record["metadata"]["version"]
(zfile,) = record["files"]
zip_path = RAW_DIR / "prior_v1.0.2.zip"
download(zfile["links"]["self"], zip_path)
algo, want = zfile["checksum"].split(":")
got = file_hash(zip_path, algo)
if got != want:
    raise RuntimeError(f"Zenodo archive checksum mismatch: {got} != {want}")
with zipfile.ZipFile(zip_path) as z:
    (top,) = {n.split("/")[0] for n in z.namelist()}
    z.extractall(RAW_DIR)
if PRIOR_DIR.exists():
    import shutil
    shutil.rmtree(PRIOR_DIR)
(RAW_DIR / top).rename(PRIOR_DIR)
print(f"prior code v1.0.2 ({zfile['key']}, {algo} verified) -> {PRIOR_DIR}")

# %% [markdown]
# ## Field data (Mendeley deposit)
#
# Same download and verification as v1.0.2 `01_data_download.py`.

# %%
DATASET_ID, VERSION = "8wvtxg53ry", 1
MENDELEY_API = "https://data.mendeley.com/public-api"
MENDELEY_DIR = RAW_DIR / "mendeley"
mzip = RAW_DIR / f"mendeley_{DATASET_ID}_v{VERSION}.zip"


def expected_checksums() -> dict[str, str]:
    def files_in(folder_id: str) -> list[dict]:
        r = requests.get(f"{MENDELEY_API}/datasets/{DATASET_ID}/files",
                         params={"folder_id": folder_id, "version": VERSION}, timeout=60)
        r.raise_for_status()
        return r.json()

    out = {f["filename"]: f["content_details"]["sha256_hash"] for f in files_in("root")}
    folders = requests.get(f"{MENDELEY_API}/datasets/{DATASET_ID}/folders/{VERSION}", timeout=60)
    folders.raise_for_status()
    for folder in folders.json():
        for f in files_in(folder["id"]):
            out[f"{folder['name']}/{f['filename']}"] = f["content_details"]["sha256_hash"]
    return out


expected = expected_checksums()
download(f"{MENDELEY_API}/zip/{DATASET_ID}/download/{VERSION}", mzip)
with zipfile.ZipFile(mzip) as z:
    z.extractall(MENDELEY_DIR)

extracted = {p.relative_to(MENDELEY_DIR).as_posix(): p for p in MENDELEY_DIR.rglob("*") if p.is_file()}
problems = []
for rel, want in sorted(expected.items()):
    matches = [p for k, p in extracted.items() if k == rel or k.endswith("/" + rel)]
    if len(matches) != 1:
        problems.append(f"missing or ambiguous: {rel}")
    elif file_hash(matches[0], "sha256") != want:
        problems.append(f"checksum mismatch: {rel}")
if problems:
    raise RuntimeError("Deposit verification failed:\n" + "\n".join(problems))
print(f"Mendeley deposit: {len(expected)} files verified")

# %%
sources = [
    {"name": "langhovde-meltwater-replication v1.0.2 (prior chain code)",
     "doi": "10.5281/zenodo.23257925", "license": "MIT", "checksum": zfile["checksum"]},
    {"name": "Sugiyama et al. 2026 data deposit (Mendeley Data)", "doi": "10.17632/8wvtxg53ry.1",
     "license": "CC-BY-4.0", "checksum": "per-file SHA-256 from the Mendeley API, all verified"},
]
(RAW_DIR / "sources.json").write_text(json.dumps(sources, indent=2))
