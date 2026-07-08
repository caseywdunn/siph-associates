# siph_associates

Cross-study mining of non-host sequence — **bacteria, viruses, eukaryotic
parasites, and prey** — from whole-genome shotgun data of siphonophores
(Cnidaria: Hydrozoa). This is the public analysis repository for the paper; it
pools three WGS datasets and analyzes them under one symmetric pipeline in which
each study is a first-class covariate, not a hierarchy.

| Study | Data | Role |
|-------|------|------|
| Church et al. 2025, *Curr. Biol.* ([10.1016/j.cub.2025.05.066](https://doi.org/10.1016/j.cub.2025.05.066)) | 151 *Physalia* libraries (5 spp.) + chromosome-scale *P. physalis* genome | Physalia sampling; a host reference |
| Ahuja et al. 2024, *GBE* ([10.1093/gbe/evae048](https://doi.org/10.1093/gbe/evae048)) | 32-species genome skim across the siphonophore phylogeny | phylogenetic breadth |
| Ahuja et al. 2026, *PLoS One* ([10.1371/journal.pone.0351247](https://doi.org/10.1371/journal.pone.0351247)) | *Nanomia* population set + chromosome-scale *N. septata* genome | Nanomia sampling; a host reference |

## Repository layout

```
README.md              this file
PLAN.md                the analysis plan — symmetric pipeline (S0–S9) + locked parameters
manifest.csv           unified specimen/library manifest (206 libraries; the S0 input)
scripts/
  build_manifest.py    reproducibly builds manifest.csv from data/sources/
data/sources/          source supplements from the three studies (provenance)
data/results/          analysis outputs (not tracked)
```

## The manifest

[`manifest.csv`](manifest.csv) is the single source of truth: **one row per
sequencing library** (the analysis unit), with `specimen_id` (YPM voucher)
grouping libraries of the same animal and `also_in_studies` flagging specimens
shared across studies. It carries identity/provenance, taxonomy, collection data,
and sequencing/data columns (SRA run, BioProject, McCleary path). Cells still
needing external gathering are marked `TODO:<what>` (SRA runs for Church/Ahuja
2024, Church raw-read paths, congener CO1 species, lat/long normalization).

Rebuild it with:

```bash
python3 scripts/build_manifest.py
```

## Status

In preparation. The pipeline (PLAN.md) is being built from a completed
exploratory pilot; the manifest is the first committed artifact. The manuscript
(LaTeX) lives in a separate repository and cites this one.

## Citation / license

TBD (add Zenodo DOI on release).
