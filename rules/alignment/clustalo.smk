"""
module to operate muscle alignments between pruned leaves and pruned alignments
1) build hmm profile from pruned alignment
2) align reads to profile
3) convert psiblast output alignment to fasta alignment
"""

__author__ = "Benjamin Linard, Nikolai Romashchenko"
__license__ = "MIT"


import pewo.config as cfg

_work_dir = cfg.get_workdir(config)


rule clustalo_profile_alignment:
    '''
    align to profile
    '''
    input:
        align=_work_dir+"/A/{pruning}.align",
        reads=_work_dir+"/R/{pruning}_r{length}.fasta"
    output:
        temp(_work_dir+"/HMM/{pruning}_r{length}.clustalo.fasta")
    version: "1.0"
    log:
        _work_dir+"/logs/clustalo/{pruning}_r{length}.log"
    benchmark:
        repeat(_work_dir+"/benchmarks/{pruning}_r{length}_clustalo_benchmark.tsv", config["repeats"])
    params:
        states=["DNA"] if config["states"]==0 else ["Protein"]
    shell:
        #for some reason, clustalo require profile-profile when query is just 1 sequence
        """
        n=$( cat input.fasta | grep '>' ) ;
        if [ $n -gt 1 ];
        then
            clustalo --infmt=fa --outfmt=fa --seqtype={params.states}
            -i {input.reads} -p1 {input.align} -o {output}  
            ;
        else
            clustalo --infmt=fa --outfmt=fa --seqtype={params.states}
            -p2 {input.reads} -p1 {input.align} -o {output}  
        fi
        """

rule clustalo_split_alignment:
    '''
    split hmm alignment results in "query only" and "reference alignment only" sub-alignments
    contrary to other placement software, such input is required by epa-ng
    '''
    input:
        align=_work_dir+"/HMM/{pruning}_r{length}.clustalo.fasta",
        reads=_work_dir+"/R/{pruning}_r{length}.fasta"
    output:
        _work_dir+"/HMM/{pruning}_r{length}.fasta_queries",
        _work_dir+"/HMM/{pruning}_r{length}.fasta_refs"
    version: "1.0"
    shell:
        "pewo/split_hmm_alignment.py {input.reads} {input.align}"
