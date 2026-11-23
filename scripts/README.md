# README - `./scripts`

`./scripts` is where all code scripts, regardless of language (e.g. R, python), are to be stored.

## md5 checksums

Before running the workflow scripts, `scripts/00_md5sum_checks.sh` should be ran using:

```bash
bash scripts/00_md5sum_checks.sh                                          ─╯
md5sum -c data/metadata/data-checksums.txt
```

Which will recursively find files and generate a checksum for the `data/` and `results/` files.

## Workflow scripts

Workflow scripts depend on functions sourced from `scripts/functions.R`.

For each output, there will be a numbered set of scripts to be run sequentially, as indicated by the number in the name.

For the **best candidate outputs** the order of running is:

1. `scripts/10-create-output-best-candidate.R`
2. `scripts/11-create-output-best-candidate.R`

