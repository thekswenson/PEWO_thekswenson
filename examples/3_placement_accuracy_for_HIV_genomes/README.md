# PEWO demo n°3

## Overview

This demo measures placement accuracy in terms of Node Distance (ND)
and expected Node Distance (eND)for a reference dataset
of 104 HIV complete genomes.

EPA-ng, PPlacer, RAPPAS are tested.

Only 10 pruning are launched and only default parameters are tested.
This analysis will require around 1 hours of computation.

A better analysis would ask for >50 prunings; to generate a wide
range of topologies (1 leaf pruned, large clades pruned, ...).


## How to launch

Download pipeline.
```
git clone --recursive https://github.com/phylo42/PEWO.git
cd PEWO
```

Execute installation script.
```
chmod u+x INSTALL.sh
./INSTALL.sh
```

After installation, load environement.
```
conda activate PEWO
```

Test workflow before launch.
```
snakemake -n \
--snakefile eval_accuracy.smk \
--configfile examples/3_placement_accuracy_for_HIV_genomes/config.yaml
```

Execute workflow, using 2 CPU cores and 8Gb of RAM.
```
snakemake --cores 2 --resources mem_mb=8000 \
--snakefile eval_accuracy.smk \
--configfile examples/3_placement_accuracy_for_HIV_genomes/config.yaml
```

## Comments

Raw results will be written in
'examples/3_placement_accuracy_for_HIV_genomes/run'.

Results summaries and plots will be written in
'examples/3_placement_accuracy_for_HIV_genomes/run'.

See PEWO wiki for a more detailed explanation of the results:
https://github.com/phylo42/PEWO/wiki/IV.-Tutorials-and-results-interpretation
