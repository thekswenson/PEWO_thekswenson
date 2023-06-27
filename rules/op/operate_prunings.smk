"""
This module computes different prunings using the alignment/tree set from the config file.
"""

__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"

# TODO: keep or remove read length sd ?

import os
import pewo.config as cfg
from pathlib import Path
from typing import Dict
from pewo.templates import get_common_queryname_template, get_common_queryname_re

_work_dir = cfg.get_work_dir(config)


def get_input_reads():
    """
    Creates an input read file parameter.
    """
    return config["dataset_reads"] if "dataset_reads" in config else []


def get_pruning_output_read_files():
    """
    Creates a list of output read files.
    """
    output_directory = os.path.join(_work_dir, "R")

    filename = f"{get_common_queryname_template(config)}.fasta"
    return [os.path.join(output_directory,
                         filename.format(pruning=pruning, length=length))
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
    Runs PrunedTreeGenerator to creates prunings.
    """
    input:
        a = config["dataset_align"],
        t = config["dataset_tree"],
        r = get_input_reads()
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
        count = config["pruning_count"],
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
                "{params.count} {params.length} 0 1 {params.states} "
                "&> {log}"
            )
                    #Rename the read files appropriately:
            namere = get_common_queryname_re(config)
            for readfile in get_pruning_output_read_files():
                m = namere.match(readfile)
                if m:
                    prunnum, length = m.groups()
                    assert length == params.length, f"bad length in {readfile}."
                else:
                    raise ValueError(f"Unexpected read file name {readfile}.")

                Path(f"{params.wd}/R/{prunnum}_r{length}.fasta").rename(readfile)

        else:
            shell(
                "mkdir -p {params.wd}/A {params.wd}/T {params.wd}/G {params.wd}/R;"
                "cp {input.a} {params.wd}/A/0.align;"
                "cp {input.t} {params.wd}/T/0.tree;"
                "touch {params.wd}/G/0.fasta;"
            )
