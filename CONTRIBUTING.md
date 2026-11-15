# Contributing

Bug reports, questions and pull requests are welcome.

## Reporting a problem

Open an issue and include:
- the command you ran and the profile you used;
- your Nextflow version (`nextflow -version`);
- the error message, and `.nextflow.log` if possible.

## Development setup

```bash
git clone https://github.com/AKULA-PRASAD/luad-rnaseq-nextflow.git
cd luad-rnaseq-nextflow
python3 tests/make_test_data.py          # regenerates tests/data*/ (deterministic)
```

## Before opening a pull request

Run the same checks as CI:

```bash
nextflow lint . -exclude nf-test.config -exclude .nf-test -exclude work -exclude .git
python3 -m unittest discover -s tests -p 'test_*.py' -v
nextflow run main.nf -profile test,docker
python3 tests/check_outputs.py results_test
nextflow run main.nf -profile test_paired,docker
python3 tests/check_outputs.py results_test_paired tests/data_paired/samplesheet.csv
nf-test test                              # module and pipeline tests (needs nf-test >= 0.9)
```

If you change an output on purpose, update its snapshot with `nf-test test --update-snapshot`, and explain the change in the pull request.

On Apple Silicon, add `-process.containerOptions --platform=linux/amd64` to the `nextflow run` commands, and run nf-test as `NFT_PLATFORM=linux/amd64 nf-test test`.

## Conventions

- **Modules:** one process per file in `modules/`. Each module pins its tool version both as a container and as a conda package, and writes `versions.yml`.
- **Helper scripts:** scripts go in `bin/` and must be executable. Python helpers use only the standard library.
- **New parameters:** add each to `nextflow.config` with a comment, and to the parameter table in `README.md`.
- **Tests:** each module has a test in `tests/modules/`. A new module needs one too.
- **Changelog:** record user-visible changes in `CHANGELOG.md`.
- **Test data:** change it only through `tests/make_test_data.py`, and commit the regenerated files.
