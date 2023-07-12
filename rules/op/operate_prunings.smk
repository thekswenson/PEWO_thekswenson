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
from pewo.io.fasta import split_fasta

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
                "java -cp {params.jar} PrunedTreeGenerator_LITE "
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
    read length (e.g. 0_r150.fasta).  We convert them to names expected
    by the snakemake workflow (e.g. 0_partition_r150.fasta).
    """
    input:
        "{d}/{counter}_r{length}.fasta"

    output:
        "{d}/{counter}_partition_r{length}.fasta"

    shell:
        "mv {input} {output}"


rule simulate_ART_reads:
    """
    This uses simulation tools to generate ART reads, rather than
    paritioning the leaf sequences as done in by the PrunedTreeGenerator.
    """
    input:
        _work_dir + "/G/{pruning}.fasta"

    output:
        _work_dir + "/R/{pruning}_ART-{platform}-{sequencer}_r{length}.fasta"

    log:
        _work_dir + "/logs/read_simulation/{pruning}_ART-{platform}-{sequencer}_r{length}.log"

    params:
        wd = _work_dir,
        reps = config["ART_gen"]["reps"],
        length = get_params_length()

    run:
        out_dir = Path(params.wd, "R")
        tmp_pref = str(out_dir / _READS_TMPPREF)
        tmp_pref += (f"{wildcards.pruning}_ART-{wildcards.platform}-"
                     f"{wildcards.sequencer}_r{wildcards.length}")
        if wildcards.platform == "illumina":
            shell(
                "art_illumina -na -sam -ss {wildcards.sequencer} -i {input} "
                "-l {params.length} -c {params.reps} -o {tmp_pref} "
                " &> {log}"
            )
        else:
            raise(ValueError(f"Unexpected platform specifier: "
                             f'"{wilcards.platform}".'))

        seqs = SeqIO.parse(tmp_pref + '.fq', format='fastq')
        SeqIO.write(seqs, output[0], format='fasta')

        for filename in glob.glob(f"{tmp_pref}*"):
            os.remove(filename) 