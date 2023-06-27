"""
This module computes extended trees and AR files, which are inputs to RAPPAS db_build.
this first is generally managed by RAPPAS, but separated here to facilitate analysis of very large datasets
which may generate ressource-consuming ARs (e.g., phyml requiring lots of RAM)
"""

__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"

from pewo.software import get_ar_binary
from pewo.templates import get_ar_output_templates
import pewo.config as cfg

_work_dir = cfg.get_work_dir(config)

rule compute_ar_inputs:
    """
    Prepares AR outputs using RAPPAS.
    """
    input:
        a=_work_dir+"/A/{pruning}.align",
        t=_work_dir+"/T/{pruning}.tree"
    output:
        _work_dir+"/RAPPAS/{pruning}/red{reduction}_ar{arsoft}/extended_trees/extended_align.fasta",
        _work_dir+"/RAPPAS/{pruning}/red{reduction}_ar{arsoft}/extended_trees/extended_align.phylip",
        _work_dir+"/RAPPAS/{pruning}/red{reduction}_ar{arsoft}/extended_trees/extended_tree_withBL.tree",
        _work_dir+"/RAPPAS/{pruning}/red{reduction}_ar{arsoft}/extended_trees/extended_tree_withBL_withoutInterLabels.tree"
    log:
        _work_dir+"/logs/ar_inputs/{pruning}_red{reduction}_ar{arsoft}.log"
    params:
        states = ["nucl"] if config["states"]==0 else ["amino"],
        reduc = config["config_rappas"]["reduction"],
        workdir = _work_dir+"/RAPPAS/{pruning}/red{reduction}_ar{arsoft}",
        arbin = lambda wildcards: get_ar_binary(config, wildcards.arsoft)
    version: "1.00"
    shell:
        "java -Xms2G -jar $(which RAPPAS.jar) -p b -b $(which {params.arbin}) "
        "-t {input.t} -r {input.a} "
        "-w {params.workdir} -s {params.states} --ratio-reduction {wildcards.reduction} "
        "--use_unrooted --arinputonly &> {log}"

rule ar_phyml:
    input:
        a = _work_dir+"/RAPPAS/{pruning}/red{red}_arPHYML/extended_trees/extended_align.phylip",
        t = _work_dir+"/RAPPAS/{pruning}/red{red}_arPHYML/extended_trees/extended_tree_withBL_withoutInterLabels.tree",
        s = _work_dir+"/T/{pruning}_optimised.info"
    output: get_ar_output_templates(config, "PHYML")
    log:
        _work_dir+"/logs/ar/{pruning}_red{red}_arPHYML.log"
    benchmark:
        repeat(_work_dir+"/benchmarks/{pruning}_red{red}_arPHYML_ansrec_benchmark.tsv", config["repeats"])
    params:
        outname = os.path.join(_work_dir, "RAPPAS", "{pruning}", "red{red}_arPHYML"),
        c = config["phylo_params"]["categories"],
    run:
        phylo_params = extract_params(input.s)  #launch 1 extract per pruning
        states = "nt" if config["states"] == 0 else "aa"
        arbin = get_ar_binary(config,"PHYML")

        ar_command = arbin + \
            " --ancestral " \
            " --no_memory_check " \
            " --leave_duplicates " \
            " -d " + states + \
            " -f e " + \
            " -o r " + \
            " -b 0 " + \
            " -v 0.0 " + \
            " -i " + input.a + \
            " -u " + input.t + \
            " -c " + str(params.c) + \
            " -m " + select_model_phymlstyle() + \
            " -a " + str(phylo_params['alpha']) + \
            " &> " + str(log)
        commands = [
            ar_command,
            "mkdir -p {params.outname}/AR",
            "mv {params.outname}/extended_trees/extended_align.phylip_phyml_ancestral_seq.txt {params.outname}/AR/extended_align.phylip_phyml_ancestral_seq.txt",
            "mv {params.outname}/extended_trees/extended_align.phylip_phyml_ancestral_tree.txt {params.outname}/AR/extended_align.phylip_phyml_ancestral_tree.txt",
            "mv {params.outname}/extended_trees/extended_align.phylip_phyml_stats.txt {params.outname}/AR/extended_align.phylip_phyml_stats.txt",
            "mv {params.outname}/extended_trees/extended_align.phylip_phyml_tree.txt {params.outname}/AR/extended_align.phylip_phyml_tree.txt",
        ]
        shell(";\n".join(command for command in commands))

rule ar_raxmlng:
    input:
        a=_work_dir+"/RAPPAS/{pruning}/red{red}_arRAXMLNG/extended_trees/extended_align.phylip",
        t=_work_dir+"/RAPPAS/{pruning}/red{red}_arRAXMLNG/extended_trees/extended_tree_withBL_withoutInterLabels.tree",
        s=_work_dir+"/T/{pruning}_optimised.info"
    output:
        get_ar_output_templates(config, "RAXMLNG")
    log:
        _work_dir+"/logs/ar/{pruning}_red{red}_arRAXMLNG.log"
    benchmark:
        repeat(_work_dir+"/benchmarks/{pruning}_red{red}_arRAXMLNG_ansrec_benchmark.tsv", config["repeats"])
    params:
        outname=_work_dir+"/RAPPAS/{pruning}/red{red}_arRAXMLNG",
        c=config["phylo_params"]["categories"],
    run:
        phylo_params=extract_params(input.s)  #launch 1 extract per pruning
        states="DNA" if config["states"]==0 else "AA"
        arbin=get_ar_binary(config, "RAXMLNG")
        model=select_model_phymlstyle()+"+G"+str(config["phylo_params"]["categories"])+"{{"+str(phylo_params['alpha'])+"}}+IU{{0}}+FC"
        shell(
            arbin+" --ancestral --redo --precision 9 --seed 1 --force msa --data-type "+states+" "
            "--threads "+str(config["config_rappas"]["arthreads"])+" "                                                                                                   
            "--msa {input.a} --tree {input.t} --model "+model+" "
            "--blopt nr_safe --opt-model on --opt-branches off &> {log} ;"
            """
            mkdir -p {params.outname}/AR
            mv {params.outname}/extended_trees/extended_align.phylip.raxml.log {params.outname}/AR/extended_align.phylip.raxml.log
            mv {params.outname}/extended_trees/extended_align.phylip.raxml.ancestralTree {params.outname}/AR/extended_align.phylip.raxml.ancestralTree
            mv {params.outname}/extended_trees/extended_align.phylip.raxml.ancestralProbs {params.outname}/AR/extended_align.phylip.raxml.ancestralProbs
            mv {params.outname}/extended_trees/extended_align.phylip.raxml.startTree {params.outname}/AR/extended_align.phylip.raxml.startTree
            mv {params.outname}/extended_trees/extended_align.phylip.raxml.ancestralStates {params.outname}/AR/extended_align.phylip.raxml.ancestralStates
            mv {params.outname}/extended_trees/extended_align.phylip.raxml.rba {params.outname}/AR/extended_align.phylip.raxml.rba
            """
        )

rule ar_paml:
    input:
        a=_work_dir+"/RAPPAS/{pruning}/red{red}_arPAML/extended_trees/extended_align.phylip",
        t=_work_dir+"/RAPPAS/{pruning}/red{red}_arPAML/extended_trees/extended_tree_withBL_withoutInterLabels.tree",
        s=_work_dir+"/T/{pruning}_optimised.info"
    output:
        get_ar_output_templates(config, "PAML")
    log:
        _work_dir+"/logs/ar/{pruning}_red{red}_arPAML.log"
    benchmark:
        repeat(_work_dir+"/benchmarks/{pruning}_red{red}_arPAML_ansrec_benchmark.tsv", config["repeats"])
    params:
        outname=_work_dir+"/RAPPAS/{pruning}/red{red}_arPAML",
        c=config["phylo_params"]["categories"],
    run:
        arbin = get_ar_binary(config, "PAML")
        shell(
            "mkdir -p {params.outname}/AR ; "
            "cd {params.outname}/AR ;"
            " " + arbin + " "+arbin+".ctl --stdout-no-buf &> {log} ;"
        )
