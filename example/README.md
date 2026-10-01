# Example

`run.sh` runs SAGAconf from start to end on the sample data that comes with ChromHMM. It
downloads ChromHMM, trains a 10-state model, parses the posteriors of GM12878 (base) and
K562 (verification) on chr11, and writes a full report to `example/sagaconf_base/`. The two
cell types are not replicates; the example only shows the steps.

```bash
pip install sagaconf
bash example/run.sh
```

The script needs Java, `curl` and `unzip`. `example_mnemonics.txt` is an example mnemonics
file.
