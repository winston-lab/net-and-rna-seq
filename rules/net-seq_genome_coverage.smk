#!/usr/bin/env python

rule genome_coverage:
    input:
        lambda wc:  {True:
                        {   True : f"alignment/{wc.sample}_{ASSAY}-noPCRduplicates-" + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
                            False: f"alignment/{wc.sample}_{ASSAY}-noPCRduplicates.bam"
                        },
                    False:
                        {   True: f"alignment/{wc.sample}_{ASSAY}-uniquemappers-" + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
                            False: f"alignment/{wc.sample}_{ASSAY}-uniquemappers.bam"
                        }
                    }.get(config["random-hexamer"]).get(len(SISAMPLES)>0)
    output:
        "coverage/{counttype}/{sample}_{ASSAY}-{readtype}-{counttype}-{strand}.bedgraph",
    params:
        strand = lambda wc: {"plus": "-", "minus": "+"}.get(wc.strand),
        split = lambda wc: {"5end": "", "wholeread": "-split"}.get(wc.readtype),
        end = lambda wc: {"5end": "-5", "wholeread": ""}.get(wc.readtype)
    wildcard_constraints:
        counttype="counts|sicounts",
        strand="plus|minus"
    log: "logs/genome_coverage/genome_coverage_{sample}-{counttype}-{readtype}-{strand}-{ASSAY}.log"
    shell: """
        (bedtools genomecov -bga {params.end} -strand {params.strand} {params.split} -ibam {input} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_{ASSAY}-{readtype}-counts-{strand}.bedgraph",
        bam = lambda wc: {True:
                            {True :
                                {"libsizenorm": "alignment/{sample}_{ASSAY}-noPCRduplicates-experimental.bam",
                                 "spikenorm"  : "alignment/{sample}_{ASSAY}-noPCRduplicates-spikein.bam"},
                             False:
                                {"libsizenorm": "alignment/{sample}_{ASSAY}-noPCRduplicates.bam"}
                            },
                          False:
                            {True :
                                {"libsizenorm": "alignment/{sample}_{ASSAY}-uniquemappers-experimental.bam",
                                 "spikenorm"  : "alignment/{sample}_{ASSAY}-uniquemappers-spikein.bam"},
                             False:
                                {"libsizenorm": "alignment/{sample}_{ASSAY}-uniquemappers.bam"}
                            }
                         }.get(config["random-hexamer"]).get(len(SISAMPLES)>0).get(wc.norm).format(**wc)
    output:
        normalized = "coverage/{norm}/{sample}_{ASSAY}-{readtype}-{norm}-{strand}.bedgraph",
    params:
        scale_factor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log: "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{readtype}-{strand}-{ASSAY}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam} | paste -d "" - <(echo "/({params.scale_factor}*1000000)") | bc -l) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = lambda wc: "coverage/{norm}/{sample}_{ASSAY}-{readtype}-{norm}-".format(**wc) + ("plus" if wc.strand=="SENSE" else "minus") +".bedgraph",
        minus = lambda wc: "coverage/{norm}/{sample}_{ASSAY}-{readtype}-{norm}-".format(**wc) + ("minus" if wc.strand=="SENSE" else "plus") +".bedgraph",
    output:
        "coverage/{norm}/{sample}_{ASSAY}-{readtype}-{norm}-{strand}.bedgraph",
    wildcard_constraints:
        strand="SENSE|ANTISENSE"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{readtype}-{norm}-{strand}-{ASSAY}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus} > {output}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = "coverage/{norm}/{sample}_{ASSAY}-{readtype}-{norm}-{strand}.bedgraph",
        fasta = lambda wc: config["genome"]["spikein_fasta"] if wc.norm=="sicounts" else config["genome"]["fasta"]
    output:
        "coverage/{norm}/{sample}_{ASSAY}-{readtype}-{norm}-{strand}.bw",
    params:
        stranded = lambda wc: [] if wc.strand in ["plus", "minus"] else """| awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2; print $1"-minus", $2}}' | LC_COLLATE=C sort -k1,1"""
    log : "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{readtype}-{norm}-{strand}-{ASSAY}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes {params.stranded}) {output}) &> {log}
        """

