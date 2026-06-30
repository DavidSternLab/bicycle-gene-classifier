# bicycle-gene-classifier

A logistic-regression classifier that identifies **bicycle genes** in eukaryotic
genome annotations using only gene-structure features (exon counts, exon
lengths, intron phase distribution). Trained on *Hormaphis cornu* (Hcor) and
shown to recover bicycle homologs in other species.

Background and the published model: Stern DL & Han Y, 2022.
*Genome Biology and Evolution* **14**:evac069. [PMID 35660862](https://pmc.ncbi.nlm.nih.gov/articles/PMC9168663/) · [DOI 10.1093/gbe/evac069](https://doi.org/10.1093/gbe/evac069).

---

## What you need

- **R ≥ 4.0** (tested with 4.2.3)
- The R packages installed by `install.R`:
  CRAN — `dplyr`, `tidyr`, `ggplot2`, `optparse`
  Bioconductor — `rtracklayer`
- A trained model file (`Hcor.glm.full_v5.5.6`, ~1.3 MB). See
  [Downloading the model](#downloading-the-model).
- An input annotation in GFF3 (or GTF) with **CDS** features that include
  **`phase`** plus a **`Parent`**, `transcript_id`, or `gene_id` attribute.
  Genes with fewer than 3 CDS exons are dropped (the classifier needs first,
  last, and internal exon features).

## Install

```bash
git clone https://github.com/DavidSternLab/bicycle-gene-classifier.git
cd bicycle-gene-classifier
Rscript install.R                       # one-time: install R dependencies
chmod +x bin/bicycle_classifier scripts/get_model.sh tests/test_example.sh
```

## Downloading the model

The trained Hcor GLM model is distributed separately from the code so that
the repository stays small and the model can be cited with a DOI.

> **Status (2026):** the model file has not yet been uploaded to Zenodo.
> Until it is, get the file directly from the lab and either pass it with
> `-m /path/to/model` or set `BICYCLE_MODEL=/path/to/model`.

Once the Zenodo deposit exists, you'll be able to fetch it with one command:

```bash
bin/bicycle_classifier --download-model     # → $HOME/.bicycle-classifier/models/
```

The classifier looks for a model in this order:
1. `-m / --model` flag (explicit path)
2. `$BICYCLE_MODEL` environment variable
3. `$HOME/.bicycle-classifier/models/Hcor.glm.full_v5.5.6` (download target)
4. `<repo>/models/Hcor.glm.full_v5.5.6` (if present locally)

## Usage

```bash
bin/bicycle_classifier -g my_genes.gff3 -o my_species -c 0.72 -d results/
```

| Flag | Default | Description |
|---|---|---|
| `-g, --gff` | (required) | Input GFF3/GTF |
| `-m, --model` | from env / cache | Path to trained `.rda` model |
| `-c, --cutoff` | `0.72` | Classification threshold (0–1) |
| `-o, --output` | `bicycle_output` | Output file prefix |
| `-d, --outdir` | `bicycle_results` | Output directory |
| `--download-model` | — | Fetch the default model and exit |
| `-h, --help` | — | Show help |

Three output files land in `--outdir`:

- `<prefix>_classifier_all_transcripts_response.txt` — every transcript with its predicted probability
- `<prefix>_classifier_bicycle_gene_names.txt` — gene names with probability ≥ cutoff (one per line)
- `<prefix>_classifier_response_histogram.pdf` — distribution plot with cutoff overlay

## Quick smoke test

A tiny synthetic GFF3 is bundled at `data/example.gff3` for confirming the
install is wired correctly. **It is not biologically meaningful** — only
useful for "does the pipeline run end-to-end."

```bash
# After install.R and a model is available
tests/test_example.sh
```

The test exits 0 on success, 77 if it had to skip (no model resolvable),
non-zero on real failure.

## Filtering a GTF by gene list

`bin/filter_gtf` is a tiny awk helper for keeping only the genes called as
bicycle by the classifier:

```bash
bin/filter_gtf <gene_list.txt> <input.gtf> <output.gtf>
```

## Repository layout

```
bicycle-gene-classifier/
├── README.md
├── CITATION.cff
├── install.R                       # one-shot R dependency installer
├── bin/
│   ├── bicycle_classifier          # bash wrapper (entry point)
│   └── filter_gtf                  # helper: subset GTF by gene list
├── R/
│   └── bicycle_classifier.R        # the actual classifier
├── scripts/
│   └── get_model.sh                # downloads model from Zenodo (URL TBD)
├── data/
│   └── example.gff3                # synthetic input for smoke testing
├── tests/
│   └── test_example.sh             # end-to-end smoke test
└── models/                         # (gitignored) downloaded models live here
```

## Cite

```
Stern DL & Han Y. Genetic Innovations in Aphids' Salivary Gland Effectors via
Convergent Evolution Identified by Gene-Structure-Based Search.
Genome Biol Evol. 2022;14(6):evac069. doi:10.1093/gbe/evac069. PMID:35660862.
```

A `CITATION.cff` is included so GitHub's "Cite this repository" button picks
it up automatically.
