# PEWO demo n°2

## Overview

This demo measures placement accuracy in terms of Node Distance (ND)
and expected Node Distance (eND) for a reference dataset
of 150 16S-rRNA barcodes.

EPA-ng, PPlacer, RAPPAS and Apples are tested.

Only 10 prunings are executeed and for a set of parameters in each program.
This analysis will require around 2 hours of computation.

A better analysis would require for >50 prunings to generate a wide
range of topologies (1 leaf pruned, large clades pruned, ...).


## How to run the pipeline
Download the pipeline.
``` bash
git clone --recursive https://github.com/phylo42/PEWO.git
cd PEWO
```

Execute the installation script
``` bash
chmod u+x INSTALL.sh
./INSTALL.sh
```

After installation, load the environment.
``` bash
conda activate PEWO
```

Test workflow before execution.
``` bash
snakemake -n \
--snakefile eval_accuracy.smk \
--configfile examples/2_placement_accuracy_for_a_bacterial_taxonomy/config.yaml
```

Execute workflow, using 2 CPU cores and 8Gb of RAM.
```
snakemake --cores 2 --resources mem_mb=8000 \
--snakefile eval_accuracy.smk \
--configfile examples/2_placement_accuracy_for_a_bacterial_taxonomy/config.yaml
```

## Comments

Raw results will be written in
`examples/2_placement_accuracy_for_a_bacterial_taxonomy/run`.

Results summaries and plots will be written in
`examples/2_placement_accuracy_for_a_bacterial_taxonomy/run`.

See PEWO wiki for a more detailed explanation of the results:
https://github.com/phylo42/PEWO/wiki/IV.-Tutorials-and-results-interpretation
