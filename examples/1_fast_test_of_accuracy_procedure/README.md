# PEWO Example n°1 : Fast test of Accuracy procedure

## Overview

This demo measures placement accuracy in terms of Node Distance (ND)
and expected Node Distance (eND) for a reference dataset
of 150 16S-rRNA barcodes.

EPA-ng, PPlacer and RAPPAS are run using their default parameters.
Only two prunings are launched to yield results in less than 20 minutes.

A better analysis would require for >50 prunings to generate a wide
range of topologies (1 leaf pruned, large clades pruned, ...).


## How to launch

Download pipeline.
``` bash
git clone --recursive https://github.com/phylo42/PEWO.git
cd PEWO
```

Execute installation script
``` bash
chmod u+x INSTALL.sh
./INSTALL.sh
```

After installation, load environement.
``` bash
conda activate PEWO
```

Test workflow before launch.
``` bash
snakemake -np \
--snakefile eval_accuracy.smk \
--configfile examples/1_fast_test_of_accuracy_procedure/config.yaml
```

Execute workflow, using 1 CPU core.
``` bash
snakemake -p --cores 1 \
--snakefile eval_accuracy.smk \
--configfile examples/1_fast_test_of_accuracy_procedure/config.yaml
```

## Comments

Raw results will be written in
`examples/1_fast_test_of_accuracy_procedure/run`.

Results summaries and plots will be written in
`examples/1_fast_test_of_accuracy_procedure/run`.

See PEWO wiki for a more detailed explanation of the results:
https://github.com/phylo42/PEWO/wiki/IV.-Tutorials-and-results-interpretation
