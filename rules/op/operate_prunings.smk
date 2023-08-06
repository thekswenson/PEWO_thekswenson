"""
This module computes different prunings using the alignment/tree set from the config file.
"""

__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"

# TODO: keep or remove read length sd ?

import os
import glob
from pathlib import Path
from typing import Dict
from Bio import SeqIO
import pewo.config as cfg
from pewo.templates import get_common_queryname_template

_work_dir = cfg.get_work_dir(config)

_READS_TMPPREF = Path("generated_reads")

def get_input_reads():
    """
    Creates an input read file parameter.
    """
    return config["dataset_reads"] if "dataset_reads" in config else []


def get_pruning_output_read_files():
    """
    Creates a list of output read files.
    """
    return [os.path.join(_work_dir, "R",
                         "{0}_r{1}.fasta".format(pruning, length))
            for pruning in range(config["pruning_count"])
            for length in config["read_length"]]


def _generate_reads(config: Dict) -> bool:
    return cfg.get_mode(config) == cfg.Mode.ACCURACY


def get_params_length():
    """
    Creates {params.length} parameter.
    """
    return str(config["read_length"]).replace("[","").replace("]","").replace(" ","")


rule operate_pruning:
    """
    Runs PrunedTreeGenerator to create prunings.
    """
    input:
        a = config["dataset_align"],
        t = config["dataset_tree"],
        #r = get_input_reads()
    output:
        #os.path.join(_work_dir, "Dtx.csv"),   #Won't work for likelihood
        #os.path.join(_work_dir, "D2tx.csv"),
        a = expand(_work_dir + "/A/{pruning}.align", pruning=range(config["pruning_count"])),
        t = expand(_work_dir + "/T/{pruning}.tree", pruning=range(config["pruning_count"])),
        g = expand(_work_dir + "/G/{pruning}.fasta", pruning=range(config["pruning_count"])),
        r = get_pruning_output_read_files() if _generate_reads(config) else []
    log:
        _work_dir + "/logs/operate_pruning.log"
    version:
        "1.00"
    params:
        wd = _work_dir,
        pcount = config["pruning_count"],
        states = config["states"],
        jar = config["pewo_jar"],
        length = get_params_length()
        #length_sd=config["read_length_sd"],
        #bpe=config["bpe"],
    run:
        if _generate_reads(config):
            shell(
                "java -cp {params.jar} PrunedTreeGenerator "
                "{params.wd} {input.a} {input.t} "
                "{params.pcount} {params.length} 0 1 {params.states} "
                "&> {log}"
            )
        else:
            shell(
                "mkdir -p {params.wd}/A {params.wd}/T {params.wd}/G {params.wd}/R;"
                "cp {input.a} {params.wd}/A/0.align;"
                "cp {input.t} {params.wd}/T/0.tree;"
                "touch {params.wd}/G/0.fasta;"
            )


ruleorder: rename_PrunedTreeGenerator_reads > simulate_ART_reads


rule rename_PrunedTreeGenerator_reads:
    """
    PrunedTreeGenerator creates files with names hardcoded based on the
    read length (e.g. 0_r150.fasta).  This rule renames them to those expected
    by the snakemake workflow (e.g. 0_rgenPARTITION_r150.fasta).
    """
    input:
        "{d}/{counter}_r{length}.fasta"

    output:
        "{d}/{counter}_rgen" + cfg.ReadGen.PARTITION.value + "_r{length}.fasta"

    shell:
        "mv {input} {output}"


paramre = re.compile(r"([a-z]+)([A-Z,\d]+)$")
def getSubParameters(paramstr):
    """
    Creates a dictionary of subparameters from a string delimited by '-'.
    Each paramter has a lowercase name and and uppercase value.
    """
    retval = {}
    params = paramstr.split("-")
    for param in params:
        m = paramre.match(param)
        if m:
            retval[m.group(1)] = m.group(2)
        else:
            raise(ValueError(f"Expected all lowercase before all uppercase in "
                             f"(sub)parameter format: "
                             f'"{param}" in "{paramstr}".'))

    return retval


rule simulate_ART_reads:
    """
    This uses simulation tools to generate ART reads, rather than
    paritioning the leaf sequences as done in by the PrunedTreeGenerator.
    """
    input:
        _work_dir + "/G/{pruning}.fasta"

    output:
        _work_dir + "/R/{pruning}_rgen" + cfg.ReadGen.ART.value +"-{parameters}_r{length}.fasta"

    log:
        _work_dir + "/logs/read_simulation/{pruning}_rgen" + cfg.ReadGen.ART.value + "ART-{parameters}_r{length}.log"

    params:
        wd = _work_dir,
        reps = config["ART_gen"]["reps"],
        length = get_params_length()

    run:
        subparams = getSubParameters(wildcards.parameters)

        out_dir = Path(params.wd, "R")
        tmp_pref = str(out_dir / _READS_TMPPREF)

        if subparams["co"] == "ILLUMINA":
            platform = cfg.ILLUM_PL[subparams["pl"]]

            tmp_pref += (f"{wildcards.pruning}_rgen{cfg.ReadGen.ART.value}"
                         f"-co{subparams['co']}-pl{subparams['pl']}"
                         f"_r{wildcards.length}")
            shell(
                "art_illumina -na -sam -ss {platform} -i {input} "
                "-l {wildcards.length} -c {params.reps} -o {tmp_pref} "
                " &> {log}"
            )
        else:
            raise(ValueError(f"Unexpected company specifier: "
                             f'"{wilcards.company}".'))

        seqs = SeqIO.parse(tmp_pref + '.fq', format='fastq')
        SeqIO.write(seqs, output[0], format='fasta')

        for filename in glob.glob(f"{tmp_pref}*"):
            os.remove(filename) 
