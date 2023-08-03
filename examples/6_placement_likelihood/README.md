# PEWO demo №6

## Overview

This demo measures average likelihood of modified trees produced by extending the original tree with placed sequences. 
For each query sequence, an extended tree is being constructed by inserting a new node, thus splitting the most likely
branch according to the placement of the sequence. Then, the likelihood of every extended tree is calculated. 
A reference dataset of bacterial 16S-rRNA barcodes is used.

This example will produce this evaluation from only 100 query reads.
A better analysis would involve thousands of query reads expected to be related to every regions of the reference tree.

## How to launch

Download pipeline:
```
git clone --recursive https://github.com/phylo42/PEWO.git
cd PEWO
```

Execute installation script:
```
chmod u+x INSTALL.sh
./INSTALL.sh
```

Load the environement:
```
conda activate PEWO
```

Test the workflow:
```
snakemake -n \
--snakefile eval_likelihood.smk \
--configfile examples/6_placement_likelihood/config.yaml
```

Execute workflow, using 2 CPU cores and 16Gb of RAM:
```
snakemake --cores 2 --resources mem_mb=16000 \
--snakefile eval_likelihood.smk \
--configfile examples/6_placement_likelihood/config.yaml
```

## Comments

Raw results will be written in 'examples/6_placement_likelihood/run'.

Results summaries and plots will be written in
'examples/6_placement_likelihood/run'.

See PEWO wiki for a more detailed explanation of the results:
https://github.com/phylo42/PEWO/wiki/IV.-Tutorials-and-results-interpretation

