"""
WORKFLOW TO EVALUATE PLACEMENT ACCURACY, GIVEN PARAMETERS SET IN "config.yaml"
This snakefile loads all necessary modules and builds the evaluation workflow itself
based on the setup defined in the config file.
"""

__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"

if not os.environ.get('CONDA_PREFIX'):
    print("WARNING: activate conda environment with 'conda activate PEWO'!",
          file=sys.stderr)

configfile: "config.yaml"


config["mode"] = "accuracy"

#utils
include:
    "rules/utils/workflow.smk"
#prunings
include:
    "rules/op/operate_prunings.smk"


rule all:
    '''
    top snakemake rule, necessary to launch the workflow
    '''
    input:
        build_reads_workflow()
