#!/usr/bin/env python

rule genome_coverage:
    input:
        {True:
            {   True : "alignment/{sample}_net-seq-noPCRduplicates-{species}.bam",
                False: "alignment/{sample}_net-seq-noPCRduplicates.bam"
            },
        False:
            {   True: "alignment/{sample}_net-seq-uniquemappers-{species}.bam",
                False: "alignment/{sample}_net-seq-uniquemappers.bam"
            }
        }.get(config["random-hexamer"]).get(len(SISAMPLES)>0)
    output:
        "coverage/{counttype}/{sample}_netseq-{readtype}-{counttype}-{strand}.bedgraph",
    params:
        strand = lambda wc: {"plus": "+", "minus": "-"}.get(wc.strand),
        split = lambda wc: {"5end": "", "wholeread": "-split"}.get(wc.readtype),
        end = lambda wc: {"5end": "-5", "wholeread": ""}.get(wc.readtype)
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/genome_coverage/genome_coverage_{sample}-{counttype}-{readtype}-{strand}.log"
    shell: """
        (genomeCoverageBed -bga {params.end} -strand {params.strand} {params.split} -ibam {input} > {output}) &> {log}
        """

# rule normalize:
#     input:
#         counts = "coverage/counts/{sample}-netseq-counts-{readtype}-{strand}.bedgraph",
#         plmin = lambda wc: "coverage/counts/" + wc.sample + "-netseq-counts-5end-plmin.bedgraph" if wc.norm=="libsizenorm" else "coverage/sicounts/" + wc.sample + "-netseq-sicounts-5end-plmin.bedgraph"
#     params:
#         scalefactor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
#     output:
#         normalized = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bedgraph",
#     wildcard_constraints:
#         norm="libsizenorm|spikenorm",
#         strand="plus|minus"
#     log: "logs/normalize/normalize-{sample}-{norm}-{readtype}.log"
#     shell: """
#         (bash scripts/libsizenorm.sh {input.plmin} {input.counts} {params.scalefactor} > {output.normalized}) &> {log}
#         """

# rule make_stranded_bedgraph:
#     input:
#         plus = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-plus.bedgraph",
#         minus = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-minus.bedgraph",
#     output:
#         sense = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-SENSE.bedgraph",
#         antisense = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-ANTISENSE.bedgraph",
#     log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}-{readtype}.log"
#     shell: """
#         (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
#         (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
#         """

# rule bg_to_bw:
#     input:
#         bedgraph = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bedgraph",
#         chrsizes = selectchrom
#     output:
#         "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bw",
#     log : "logs/bg_to_bw/bg_to_bw-{sample}-{norm}-{readtype}-{strand}.log"
#     shell: """
#         (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
#         """

