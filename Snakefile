# Snakefile — orchestrates the replication pipeline end-to-end.
#
# Replace the placeholder rules with your actual replication steps. The
# canonical pattern is one rule per pipeline stage, and each rule wraps a
# notebook executed via jupytext (so the notebook stays the source of truth
# and the Snakefile just sequences them).
#
# Usage:
#   snakemake --cores 1                  # run everything
#   snakemake --cores 1 -n               # dry run

NOTEBOOKS = "notebooks"
DATA = "data"
RESULTS = "results"
FIGURES = "figures"


rule all:
    input:
        f"{FIGURES}/main_result.png",
        f"{RESULTS}/lead_labels.csv",


rule data_download:
    output:
        f"{DATA}/raw/sources.json",
        f"{DATA}/raw/prior_v1.0.2/notebooks/gnss.py",
    log:
        f"{RESULTS}/logs/01_data_download.log",
    shell:
        f"cd {{NOTEBOOKS}} && jupytext --to notebook --execute 01_data_download.py 2>&1 | tee ../{{log}}"


rule data_clean:
    input:
        f"{DATA}/raw/sources.json",
    output:
        expand(f"{DATA}/clean/{{f}}.nc", f=["gnss_GNSS1", "gnss_GNSS2", "pressure_BH2201", "aws"]),
    log:
        f"{RESULTS}/logs/02_data_clean.log",
    shell:
        f"cd {{NOTEBOOKS}} && jupytext --to notebook --execute 02_data_clean.py 2>&1 | tee ../{{log}}"


rule analysis:
    input:
        expand(f"{DATA}/clean/{{f}}.nc", f=["gnss_GNSS1", "gnss_GNSS2", "pressure_BH2201", "aws"]),
        f"{NOTEBOOKS}/leadtime.py",
    output:
        f"{RESULTS}/series_15min.nc",
        f"{RESULTS}/onsets.csv",
        f"{RESULTS}/leads.csv",
        f"{RESULTS}/lead_labels.csv",
        f"{RESULTS}/xcorr.csv",
        f"{RESULTS}/nested_models.csv",
        f"{RESULTS}/robust_delay_check.csv",
    log:
        f"{RESULTS}/logs/03_analysis.log",
    shell:
        f"cd {{NOTEBOOKS}} && jupytext --to notebook --execute 03_analysis.py 2>&1 | tee ../{{log}}"


rule figures:
    input:
        f"{RESULTS}/series_15min.nc",
        f"{RESULTS}/onsets.csv",
        f"{RESULTS}/leads.csv",
        f"{RESULTS}/robust_delay_check.csv",
    output:
        f"{FIGURES}/main_result.png",
    log:
        f"{RESULTS}/logs/04_figures.log",
    shell:
        f"cd {{NOTEBOOKS}} && jupytext --to notebook --execute 04_figures.py 2>&1 | tee ../{{log}}"
