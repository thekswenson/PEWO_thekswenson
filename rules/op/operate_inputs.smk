"""
Prepare directories to run a resource test.
Basically, it uses same directory structure than a pruning test, but with a unique run '0', which represents
a full, non-pruned tree and for which repeated placements runs will be launched for a set of reads set
(set as R/O/O_r0.fasta).

This rule is called once when starting the "resources" workflow.
"""


__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"

import os
import pewo.config as cfg
from pewo.templates import get_common_queryname_template


_work_dir = cfg.get_work_dir(config)

# TODO: script of function to compute reads from alignment

def queries():
    if config["query_type"]=="user":
        return cfg.get_query_user(config)
    else:
        #compute queries from alignment
        return []



rule define_resource_inputs:
    input:
        a=config["dataset_align"],
        t=config["dataset_tree"],
        r=queries()
    output:
        aout=_work_dir+"/A/0.align",
        tout=_work_dir+"/T/0.tree",
        gout=_work_dir+"/G/0.fasta",
        rout=_work_dir+"/R/"+get_common_queryname_template(config).format(query=0,
                                                                          generator="RESOURCES",
                                                                          length=0)+".fasta"
    run:
        if not os.path.isdir(_work_dir):
            os.mkdir(_work_dir)
        if not os.path.isdir(_work_dir+"/A"):
            os.mkdir(_work_dir+"/A")
        if not os.path.isdir(_work_dir+"/T"):
            os.mkdir(_work_dir+"/T")
        if not os.path.isdir(_work_dir+"/G"):
            os.mkdir(_work_dir+"/G")
        if not os.path.isdir(_work_dir+"/R"):
            os.mkdir(_work_dir+"/R")
        shell(
            """
            cp {input.a} {output.aout}
            cp {input.t} {output.tout}
            """
        )
        if config["query_type"]=='user':
            shell(
                """
                cp {input.r} {output.rout}
                cp {input.r} {output.gout}
                """
            )
