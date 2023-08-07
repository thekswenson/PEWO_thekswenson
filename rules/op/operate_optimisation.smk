"""
This module optimises pruned trees.
"""

__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"

# TODO: add model/parameters selection in the config file, this config needs to be propagated to all placement software

import os
from pewo.config import get_work_dir

_work_dir = get_work_dir(config)


rule optimise:
    input:
        a=_work_dir+"/A/{pruning}.align",
        t=_work_dir+"/T/{pruning}.tree"
    output:
        temp(_work_dir+"/T/RAxML_binaryModelParameters.{pruning}"),
        temp(_work_dir+"/T/RAxML_log.{pruning}"),
        _work_dir+"/T/{pruning}_optimised.tree",
        _work_dir+"/T/{pruning}_optimised.info"
    log:
        _work_dir+"/logs/optimisation/{pruning}.log"
    version: "1.00"
    params:
        m=select_model_raxmlstyle(),
        c=config["phylo_params"]["categories"],
        name="{pruning}",
        raxmlname=_work_dir+"/T/RAxML_result.{pruning}",
        outname=_work_dir+"/T/{pruning}_optimised.tree",
        raxmlinfoname=_work_dir+"/T/RAxML_info.{pruning}",
        outinfoname=_work_dir+"/T/{pruning}_optimised.info",
        reduction=_work_dir+"/A/{pruning}.align.reduced",
        outdir= os.path.join(_work_dir,"T")
    shell:
        """
        raxmlHPC-SSE3 -f e -w {params.outdir} -m {params.m} -c {params.c} -s {input.a} -t {input.t} -n {params.name} &> {log}
        mv {params.raxmlname} {params.outname}
        mv {params.raxmlinfoname} {params.outinfoname}
        rm -f {params.reduction}
        """
