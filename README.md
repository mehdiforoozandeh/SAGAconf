<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/mehdiforoozandeh/SAGAconf/main/docs/graphical_abstract_dark.svg">
  <img alt="SAGAconf compares two replicate chromatin state annotations, gives each genomic bin an r-value, and keeps the reproducible calls." src="https://raw.githubusercontent.com/mehdiforoozandeh/SAGAconf/main/docs/graphical_abstract_light.svg">
</picture>

# SAGAconf

[![PyPI](https://img.shields.io/pypi/v/sagaconf)](https://pypi.org/project/sagaconf/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Paper](https://img.shields.io/badge/Genome%20Research-10.1101%2Fgr.278343.123-informational)](https://doi.org/10.1101/gr.278343.123)

SAGAconf gives a calibrated confidence score to each call in a chromatin state annotation.
Segmentation and genome annotation (SAGA) methods, such as ChromHMM and Segway, label every
genomic bin with a chromatin state. Many of these labels do not reproduce across replicates.
SAGAconf compares a base annotation with a verification annotation and gives each bin an
**r-value**: a score for how well its state call reproduces. You can then keep only the
reproducible calls for downstream analysis. SAGAconf works with any SAGA method, because it
reads only the posterior probability matrix.

- Paper: Foroozandeh Shahraki, Farahbod and Libbrecht, *Robust chromatin state annotation*,
  [Genome Research 34(3):469–483, 2024](https://doi.org/10.1101/gr.278343.123).
- Blog post: [Reproducibility as a measure of confidence](https://medium.com/@mehdiforoozandehsh/reproducibility-as-a-measure-of-confidence-7454b72f984e).

## Install

```bash
pip install sagaconf
```

SAGAconf needs Python 3.9 or later. The `active-regions` mode also needs
[bedtools](https://bedtools.readthedocs.io/) on your `PATH`.

## Quick start

Parse the posteriors of each annotation, then compare the two parsed files:

```bash
sagaconf parse --saga chmm base_posteriors/ 200 base/
sagaconf parse --saga chmm verif_posteriors/ 200 verif/
sagaconf run base/parsed_posterior.bed verif/parsed_posterior.bed results/
```

`results/r_values.bed` holds one r-value per bin. `results/confident_segments.bed` holds the
reproducible subset of the base annotation (this file needs mnemonics; see
[Known issues](#known-issues)).

To run the full example on ChromHMM's sample data, run [`example/run.sh`](example/run.sh).

## Input

SAGAconf compares two annotations of the same genome:

- The **base** annotation is the one you want to score.
- The **verification** annotation comes from a replicate. It can come from different data and
  a different model (setting S1 in the paper), the same model trained on both replicates (S2),
  or the same data with a different random initialization (S3).

`sagaconf run` reads one **parsed posterior** file per annotation. The file has one row per
genomic bin and `K + 3` columns: `chr`, `start`, `end`, and one posterior probability for
each of the `K` states. The file can be BED (tab-separated) or CSV. `sagaconf parse` writes
this file from the output of ChromHMM or Segway:

- **ChromHMM:** run `LearnModel` with `-printposterior`. Give `sagaconf parse` the directory
  of `*_posterior.txt` files (one file per cell type and chromosome).
- **Segway:** give `sagaconf parse` the directory of `posterior<i>.bedGraph` files (one file
  per state).

## Commands

### `sagaconf parse`

```
sagaconf parse --saga {chmm,segway} [--out-format {bed,csv}] POSTERIORDIR RESOLUTION OUTDIR
```

| Argument | Meaning |
|---|---|
| `POSTERIORDIR` | Directory with the SAGA model's posterior files. |
| `RESOLUTION` | Bin size of the SAGA model, in bp. |
| `OUTDIR` | Directory for `parsed_posterior.bed` (or `.csv`). SAGAconf creates it if it does not exist. |
| `--saga` | The SAGA model that wrote the posteriors: `chmm` or `segway`. |
| `--out-format` | `bed` (default) or `csv`. |

### `sagaconf run`

```
sagaconf run BASE VERIF OUTDIR [--mode MODE] [options]
```

| Option | Default | Meaning |
|---|---|---|
| `--mode` | `full` | The analysis to run. See the table below. |
| `-bm`, `--base-mnemonics FILE` | none | State names for the base annotation. |
| `-vm`, `--verif-mnemonics FILE` | none | State names for the verification annotation. |
| `-s`, `--chr21-only` | off | Analyse chr21 only, for a quick check. |
| `-w`, `--window-size BP` | 1000 | Window around each bin, in bp. |
| `-to`, `--iou-threshold` | 0.75 | IoU overlap that makes two states correspond. |
| `-tr`, `--repr-threshold` | 0.8 | r-value that makes a bin reproduced (α in the paper). |
| `-k`, `--merge-k K` | none | Merge the base states down to `K` states. |
| `--ccre-file FILE` | none | cCRE BED file. The `active-regions` mode needs it. |
| `--meuleman-file FILE` | none | Meuleman et al. DHS index for the `active-regions` mode. If you omit it, SAGAconf skips that step. |
| `-v`, `--verbose` | off | Report the analysis steps that fail. |

| Mode | What it writes |
|---|---|
| `full` | All analyses: r-values, confident segments, overlap, calibration, granularity and misalignment reports. |
| `quick` | The essential reports only. |
| `rvalues` | `r_values.bed` only. |
| `celltype` | The per-annotation analyses (overlap, calibration, granularity, misalignment). |
| `seglength` | Reproducibility as a function of position in a segment. |
| `merge` | The reports after merging the base states down to `-k` states. |
| `active-regions` | r-values genome-wide and in cCREs (and in Meuleman DHSs if you give the file). |

`python -m sagaconf` runs the same command as `sagaconf`.

### Mnemonics file

A mnemonics file gives each state a biological name. It is a tab-separated file with the
header `old` and `new`. The `old` column holds the state number, in the column order of the
parsed posterior file. The numbers can start at 0 or at 1. The `new` column holds the name.
SAGAconf shortens each name to its first 4 characters (and to 3 characters after an
underscore), so `Enhancer_low` becomes `Enha_low`.
[`example/example_mnemonics.txt`](example/example_mnemonics.txt) is an example:

```
old	new
1	Enhancer_low
2	Enhancer
3	Promoter_flanking
```

## Output

`sagaconf run` writes these files in `OUTDIR` in `full` mode:

| File | Content |
|---|---|
| `r_values.bed` | The r-value of each bin: `chr`, `start`, `end`, `MAP` (the most probable state), `r_value`. |
| `r_values_UCSC_GenomeBrowser.bed` | The r-values as a UCSC Genome Browser track. |
| `r_values_report.txt` | The average r-value, genome-wide and per state. |
| `rval_hist` | Histograms of r-values per state. |
| `confident_segments.bed` | The reproducible subset of the base annotation, one row per bin. |
| `confident_segments_dense.bed` | The same subset, with adjacent bins of one state merged into segments. |
| `ratio_robust.txt` | The fraction of bins that are reproduced, genome-wide and per state. |
| `overlap_ratio.txt` | The naive overlap per state. |
| `NMI.txt` | Mutual information, with and without the posteriors. |
| `coverages1.txt`, `coverages2.txt` | The genome coverage of each state in the base and verification annotations. |
| `heatmap`, `heatmap_w` | IoU overlap between the states of the two annotations, with `w = 0` and with `w > 0`. |
| `binned_posterior_heatmap` | IoU overlap against the binned posteriors of the base annotation. |
| `granularity`, `barplot`, `AUC_mAUC.txt` | The state-merging curve and the area under it (auSMC) per state. |
| `len_bound`, `len_bound_overall` | Overlap as a function of `w`, per state and genome-wide. |
| `calib/` | Posterior calibration curves per state. |
| `Dist_vs_Corresp/`, `Dist_vs_Corresp_3/` | Overlap and correspondence as a function of `w`. |

Each plot is written as PDF and SVG. Each `.txt` file next to a plot holds the plotted values.

## Known issues

These issues come from the version of SAGAconf that the paper used. This release keeps them,
so that its results match that version exactly. Fixes will come in later releases.

- In `full` mode, `-k` does not merge states.
- In `full` mode without mnemonics, SAGAconf does not write the genome-wide results
  (including `confident_segments.bed` and the UCSC track). Give `-bm` and `-vm` to get them.
- In `merge` mode, SAGAconf can stop with a `KeyError` after it writes the merge reports.

## Legacy interface

The original scripts still work, with the same flags and the same output:

```bash
python SAGAconf_parser.py --saga chmm POSTERIORDIR 200 OUTDIR
python SAGAconf.py [-q | --r_only | --ct_only | ...] BASE VERIF OUTDIR
```

`sagaconf run` also accepts the original flag names (`--r_only`, `--windowsize`,
`--base_mnemonics`, and so on). In the original scripts, the `--active_regions` mode reads
`src/biointerpret/GRCh38-cCREs.bed` and `src/biointerpret/Meuleman.tsv` relative to the
working directory.

The code as the paper used it is on the [`legacy`](https://github.com/mehdiforoozandeh/SAGAconf/tree/legacy)
branch (tag `v0-legacy`). That branch also keeps the scripts that produced the paper's figures.

## Tests

[`tests/golden/`](tests/golden) checks that this version gives byte-identical output to the
legacy version. It runs both versions on the same input and compares every file they write,
and everything they print. [`tests/golden/gate.py`](tests/golden/gate.py) describes the method.
The check passes on synthetic data and on the paper's GM12878 ChromHMM and Segway annotations.

## Citation

```bibtex
@article{foroozandeh2024robust,
  title   = {Robust chromatin state annotation},
  author  = {Foroozandeh Shahraki, Mehdi and Farahbod, Marjan and Libbrecht, Maxwell W.},
  journal = {Genome Research},
  volume  = {34},
  number  = {3},
  pages   = {469--483},
  year    = {2024},
  doi     = {10.1101/gr.278343.123}
}
```

![Graphical abstract of the paper](https://raw.githubusercontent.com/mehdiforoozandeh/SAGAconf/main/docs/paper_graphical_abstract.png)

## License

SAGAconf is available under the [MIT License](LICENSE).
