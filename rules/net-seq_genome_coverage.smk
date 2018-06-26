#!/usr/bin/env python

rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam" if config["random-hexamer"] else "alignment/{sample}-unique.bam"
    params:
        prefix = lambda wc: config["combinedgenome"]["experimental_prefix"] if wc.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
    output:
        plmin = "coverage/{counttype}/{sample}-netseq-{counttype}-5end-plmin.bedgraph",
        plus = "coverage/{counttype}/{sample}-netseq-{counttype}-5end-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}-netseq-{counttype}-5end-minus.bedgraph",
        pluswr = "coverage/{counttype}/{sample}-netseq-{counttype}-wholeread-plus.bedgraph",
        minuswr = "coverage/{counttype}/{sample}-netseq-{counttype}-wholeread-minus.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/get_coverage/get_coverage-{sample}-{counttype}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plmin}) &> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        (genomeCoverageBed -bga -strand - -split -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.pluswr}) &>> {log}
        (genomeCoverageBed -bga -strand + -split -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.minuswr}) &>> {log}
        """

rule normalize:
    input:
        counts = "coverage/counts/{sample}-netseq-counts-{readtype}-{strand}.bedgraph",
        plmin = lambda wc: "coverage/counts/" + wc.sample + "-netseq-counts-5end-plmin.bedgraph" if wc.norm=="libsizenorm" else "coverage/sicounts/" + wc.sample + "-netseq-sicounts-5end-plmin.bedgraph"
    params:
        scalefactor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
    output:
        normalized = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log: "logs/normalize/normalize-{sample}-{norm}-{readtype}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.plmin} {input.counts} {params.scalefactor} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-minus.bedgraph",
    output:
        sense = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-ANTISENSE.bedgraph",
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}-{readtype}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
        """

rule bg_to_bw:
    input:
        bedgraph = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bedgraph",
        chrsizes = selectchrom
    output:
        "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bw",
    log : "logs/bg_to_bw/bg_to_bw-{sample}-{norm}-{readtype}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

