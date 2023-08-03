# PEWO demo n°4

## Overview

This demo measures placement accuracy in terms of Node Distance (ND)
and expected Node Distance (eND) for two different reference trees
based on the same 1000 coleopteran mitochondria but three different
loci: cox1, 12S, cytb and 16S rRNA genes.

Four accuracy procedures are launched independantly for each locus.
The goal is to compare placement accuray measured for each locus and
deduce which reference tree appears to be a better reference for future
placements. 

EPA-ng, PPlacer, RAPPAS are tested with default parameters.

Only 10 pruning are launched and for a set of parameters in each program.
This analysis will require around 2 hours of computation.

A better analysis would require for >50 prunings to generate a wide
range of topologies (1 leaf pruned, large clades pruned, ...).


## How to launch

Download pipeline.
```
git clone --recursive https://github.com/phylo42/PEWO.git
cd PEWO
```

Execute installation script
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
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_12S/config_12S.yaml

snakemake -n \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_16S/config_16S.yaml

snakemake -n \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_cox1/config_cox1.yaml

snakemake -n \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_cob/config_cob.yaml
```

Execute workflow for 12S gene, using 2 CPU cores and 8Gb of RAM.
```
snakemake --cores 2 --resources mem_mb=8000 \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_12S/config_12S.yaml
```

Execute workflow for 16S gene, using 2 CPU cores and 8Gb of RAM.
```
snakemake --cores 2 --resources mem_mb=8000 \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_16S/config_16S.yaml
```

Execute workflow for cox1 gene, using 2 CPU cores and 8Gb of RAM.
```
snakemake --cores 2 --resources mem_mb=8000 \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_cox1/config_cox1.yaml
```

Execute workflow for cytb gene, using 2 CPU cores and 8Gb of RAM.
```
snakemake --cores 2 --resources mem_mb=8000 \
--snakefile eval_accuracy.smk \
--configfile examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_cob/config_cob.yaml
```


## Comments

Raw results will be written in

'examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_16S/run'
'examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_12S/run'
'examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_cox1/run'

Results summaries and plots will be written in

'examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_16S/run'
'examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_12S/run'
'examples/4_search_for_most_accurate_taxonomic_marker/coleoptera_cox1/run'

See PEWO wiki for a more detailed explanation of the results:

https://github.com/phylo42/PEWO/wiki/IV.-Tutorials-and-results-interpretation
