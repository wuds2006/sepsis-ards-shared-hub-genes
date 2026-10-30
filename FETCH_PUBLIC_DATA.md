# Retrieval Runbook for the Public Inputs

Every input listed in [`DATA_INVENTORY.md`](DATA_INVENTORY.md) is public and third-party hosted.
**No raw public data are redistributed in this repository.** This runbook reproduces the input
tree from scratch; the commands are the ones actually used for this project (see
`analysis/download_data.sh`, `analysis/download_data2.sh`, `analysis/extract_MR_data.sh`).

Set one variable first:

```bash
export RAW=/path/to/your/00_rawdata
mkdir -p "$RAW/geo" "$RAW/validation" "$RAW/gwas" "$RAW/eqtlgen"
```

---

## 1. GEO series matrices

Each series matrix is fetched from the NCBI FTP tree (`GSEnnnnnn` is zero-padded; the directory
is the first six digits with the last three replaced by `nnn`):

```bash
# GSE76293 (ARDS discovery, GPL570)
curl -sL -C - -o "$RAW/geo/GSE76293_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76293/matrix/GSE76293_series_matrix.txt.gz"

# GSE69528 (sepsis discovery, GPL10558)
curl -sL -C - -o "$RAW/geo/GSE69528_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE69nnn/GSE69528/matrix/GSE69528_series_matrix.txt.gz"

# GSE100159 (sepsis validation, GPL6884)
curl -sL -C - -o "$RAW/geo/GSE100159_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100159/matrix/GSE100159_series_matrix.txt.gz"

# GSE65682 (prognosis, GPL13667) — also fetch the family SOFT file for platform annotation
curl -sL -C - -o "$RAW/geo/GSE65682_series_matrix.txt.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/GSE65682_series_matrix.txt.gz"
curl -sL -C - -o "$RAW/geo/GSE65682_family.soft.gz" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/soft/GSE65682_family.soft.gz"
```

The two single-cell series are supplied as per-sample processed matrices rather than a single
series matrix; download the `suppl/` archives and unpack them:

```bash
curl -sL -C - -o "$RAW/geo/GSE151263_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151263/suppl/GSE151263_RAW.tar"
curl -sL -C - -o "$RAW/geo/GSE167363_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167363/suppl/GSE167363_RAW.tar"
```

GSE243066 (ARDS validation, PBMC RNA-seq) is downloaded as a processed-count archive:

```bash
curl -sL -C - -o "$RAW/validation/GSE243066_RAW.tar" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243066/suppl/GSE243066_RAW.tar"
```

---

## 2. GWAS summary statistics

```bash
# Sepsis, IEU OpenGWAS ID ieu-b-4980 (GWAS Catalog GCST90281174)
curl -sL -o "$RAW/gwas/ieu-b-4980.vcf.gz" \
  "https://gwas-api.mrcieu.ac.uk/files/ieu-b-4980/ieu-b-4980.vcf.gz"

# ARDS, FinnGen R11 endpoint J10_ARDS
curl -sL -C - -o "$RAW/gwas/finngen_R11_J10_ARDS.gz" \
  "https://storage.googleapis.com/finngen-public-data-r11/summary_stats/finngen_R11_J10_ARDS.gz"
```

Both are harmonised and clumped inside the pipeline; see `analysis/convert_gwas_ma.py`,
`analysis/extract_MR_4980.py`, and `analysis/match_ards_4980.py`.

---

## 3. eQTLGen cis-eQTL (SMR input)

```bash
curl -sL -C - -m 7200 \
  -o "$RAW/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" \
  "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
```

≈2.5 GB; `-C -` resumes an interrupted transfer. Converted to BESD format with SMR v1.03
(`analysis/step117_run_smr_heidi.sh`).

---

## 4. Reference panels, protein and pharmacology resources

- **1000 Genomes Phase 3 (EUR)** — used as the LD reference for the HEIDI test and LDSC.
  Local PLINK binaries are built with `--make-bld` (±2,000 kb windows) as described in
  `analysis/step117_run_smr_heidi.sh`.
- **Human Protein Atlas, Blood Atlas v24** — https://www.proteinatlas.org/ (accessed 2025-11-06).
- **PDB 3EUF (UPP1) and 7SLV (VNN1)** — https://www.rcsb.org/structure/3EUF and
  https://www.rcsb.org/structure/7SLV.
- **AlphaFold DB** — https://alphafold.ebi.ac.uk/.
- **ChEMBL release 37** — https://www.ebi.ac.uk/chembl/ (accessed 2026-09-17).
- **Open Targets, September 2026 release** — https://platform.opentargets.org/.
- **LINCS / CMap** — https://clue.io/ (requires a free account for API access).
- **MSigDB Hallmark** — https://www.gsea-msigdb.org/ (seven metabolic gene sets).

---

## 5. Verify before you rerun

```bash
bash analysis/run_all.sh --list          # 17 stages, mapped to Supplementary Methods S1-S14
bash analysis/run_all.sh --stage 1       # smoke-test the first stage
```

Raw inputs are expected under `00_rawdata/`, which is git-ignored. Third-party resources are
governed by their own terms of use; this runbook downloads them for research reuse only.
