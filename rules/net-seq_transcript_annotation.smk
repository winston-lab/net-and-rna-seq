#!/usr/bin/env python

localrules:
    merge_called_transcripts,
    annotate_merged_transcripts,

rule call_transcripts:
    input:
        bam = lambda wc: {True:
                            {   True: f"alignment/{wc.sample}_{ASSAY}-noPCRduplicates-{wc.species}.bam",
                                False: f"alignment/{wc.sample}_{ASSAY}-noPCRduplicates.bam"
                            },
                          False:
                            {   True: f"alignment/{wc.sample}_{ASSAY}-uniquemappers-{wc.species}.bam",
                                False: f"alignment/{wc.sample}_{ASSAY}-uniquemappers.bam"
                            }
                        }.get(config["molecular-barcode"]).get(len(SISAMPLES)>0)
    output:
        gtf = "transcript_annotation/{sample}_{ASSAY}-{species}-transcripts.gtf",
    params:
        library_type = "--fr" if ASSAY=="rnaseq" and config["sequence-from-5prime"] else "--rf",
        min_isoform_abundance = config["stringtie"]["min-isoform-abundance"],
        min_transcript_length = config["stringtie"]["min-transcript-length"],
        min_splice_distance = config["stringtie"]["min-splice-distance"],
        min_splicejunction_coverage = config["stringtie"]["min-splicejunction-coverage"],
        min_transcript_coverage = config["stringtie"]["min-transcript-coverage"],
        min_gap_length = config["stringtie"]["min-gap-length"]
    log:
        "logs/call_transcripts/call_transcripts-{sample}-{species}-{ASSAY}.log"
    conda:
        "../envs/stringtie.yaml"
    shell: """
        (stringtie {input.bam} -v -o {output.gtf} -p 1 {params.library_type} -l {wildcards.sample} -f {params.min_isoform_abundance} -m {params.min_transcript_length} -a {params.min_splice_distance} -j {params.min_splicejunction_coverage} -c {params.min_transcript_coverage} -g {params.min_gap_length} -M 0) &> {log}
        """

rule merge_called_transcripts:
    input:
        called = lambda wc: ["transcript_annotation/" + x + f"_{ASSAY}-{wc.species}-transcripts.gtf" for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])],
    output:
        gff = "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts.gff",
        bed = "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts.bed",
    params:
        min_transcript_length = config["stringtie"]["min-transcript-length"],
        min_gap_length = config["stringtie"]["min-gap-length"]
    conda:
        "../envs/stringtie.yaml"
    log:
        "logs/merge_called_transcripts/merge_called_transcripts_{condition}-v-{control}_{species}.log"
    shell: """
        (stringtie --merge {input.called} -v -m {params.min_transcript_length} -g {params.min_gap_length} -i -l stringtie | \
         tee {output.gff} | \
         awk 'BEGIN{{FS="\t|gene_id \\"|\\"; transcript_id";OFS="\t"}} NR>2 && $3=="transcript" {{print $1, $4-1, $5, $10, $6, $7 }}' | \
         sort -k1,1 -k2,2n | \
         bedtools merge -s -c 4,5,6 -o first -i stdin > {output.bed}) &> {log}
        """

rule annotate_merged_transcripts:
    input:
        called = "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts.bed",
        reference = lambda wc: os.path.abspath(build_annotations(config["genome"]["transcript_annotation"])) if wc.species=="experimental" else config["spike_in"]["transcript_annotation"]
    output:
        "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts-annotated.bed"
    log:
        "logs/annotate_merged_transcripts/annotate_merged_transcripts_{condition}-v-{control}_{species}.log"
    shell: """
        (bedtools intersect -s -v -a {input.called} -b {input.reference} | \
         cat - {input.reference} | \
         sort -k1,1 -k2,2n > {output}) &> {log}
        """

