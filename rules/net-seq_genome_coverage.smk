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
        counttype="counts|sicounts",
        strand="plus|minus"
    log: "logs/genome_coverage/genome_coverage_{sample}-{counttype}-{readtype}-{strand}.log"
    shell: """
        (bedtools genomecov -bga {params.end} -strand {params.strand} {params.split} -ibam {input} > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_netseq-{readtype}-counts-{strand}.bedgraph",
        bam = lambda wc: {True:
                            {True :
                                {"libsizenorm": "alignment/{sample}_net-seq-noPCRduplicates-experimental.bam",
                                 "spikenorm"  : "alignment/{sample}_net-seq-noPCRduplicates-spikein.bam"},
                             False:
                                {"libsizenorm": "alignment/{sample}_net-seq-noPCRduplicates.bam"}
                            },
                          False:
                            {True :
                                {"libsizenorm": "alignment/{sample}_net-seq-uniquemappers-experimental.bam",
                                 "spikenorm"  : "alignment/{sample}_net-seq-uniquemappers-spikein.bam"},
                             False:
                                {"libsizenorm": "alignment/{sample}_net-seq-uniquemappers.bam"}
                            }
                         }.get(config["random-hexamer"]).get(len(SISAMPLES)>0).get(wc.norm).format(**wc)
    output:
        normalized = "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-{strand}.bedgraph",
    params:
        scale_factor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log: "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{readtype}-{strand}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam} | paste -d "" - <(echo "/({params.scale_factor}*1000000)") | bc -l) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-minus.bedgraph",
    output:
        sense = "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-ANTISENSE.bedgraph",
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{readtype}-{norm}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus} > {output.antisense}) &>> {log}
        """

def select_chromsizes(wc):
    if wc.strand in ["plus", "minus"]:
        if wc.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wc.norm=="sicounts":
        return os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv"
    return os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-{strand}.bedgraph",
        chrsizes = select_chromsizes
    output:
        "coverage/{norm}/{sample}_netseq-{readtype}-{norm}-{strand}.bw",
    log : "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{readtype}-{norm}-{strand}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} {input.chrsizes} {output}) &> {log}
        """

