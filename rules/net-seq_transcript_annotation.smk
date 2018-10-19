#!/usr/bin/env python

rule call_transcripts:
    input:
        bam = f"alignment/{{sample}}_{ASSAY}-noPCRduplicates-{{species}}.bam" if config["random-hexamer"] else f"alignment/{{sample}}_{ASSAY}-uniquemappers-{{species}}.bam"
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
    conda: "../envs/stringtie.yaml"
    shell: """
        stringtie {input.bam} -v -o {output.gtf} -p 1 {params.library_type} -l {wildcards.sample} -f {params.min_isoform_abundance} -m {params.min_transcript_length} -a {params.min_splice_distance} -j {params.min_splicejunction_coverage} -c {params.min_transcript_coverage} -g {params.min_gap_length} -M 0
        """

# rule build_reference_gff:
#     input:
#         bed = lambda wc: config["genome"]["transcripts"] if wc.species=="experimental" else config["genome"]["spikein_transcripts"]
#     output:
#         gff = "transcript_annotation/{species}_reference-transcripts.gff"
#     shell: """
#         awk 'BEGIN{{FS=OFS="\t"; print "# reference transcripts from {input.bed}"}} {{print $1, "{input.bed}", "transcript", $2+1, $3, $5, $6, ".", "gene_id \\""$4"\\"; transcript_id \\""$4"\\";"}}' {input.bed} > {output.gff}
#         """

rule merge_called_transcripts:
    input:
        called = lambda wc: ["transcript_annotation/" + x + f"_{ASSAY}-{wc.species}-transcripts.gtf" for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])],
        reference = lambda wc: config["genome"]["gff"] if wc.species=="experimental" else config["genome"]["spikein_gff"]
    output:
        "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts.gff"
    params:
        min_transcript_length = config["stringtie"]["min-transcript-length"],
        min_gap_length = config["stringtie"]["min-gap-length"]
    conda: "../envs/stringtie.yaml"
    shell: """
        stringtie --merge {input.called} -v -o {output} -G {input.reference} -m {params.min_transcript_length} -g {params.min_gap_length} -i -l stringtie
        """

#TODO: give transcripts from annotation higher priority than those called by stringtie when merging overlapping
rule collapse_transcripts_to_genes:
    input:
        "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts.gff"
    output:
        "transcript_annotation/{condition}-v-{control}/{condition}-v-{control}_{species}-merged-transcripts-collapsed.bed"
    params:
        cluster_dist = -30
    shell: """
        # transform transcript entries in GFF to BED format without gene name
        awk 'BEGIN{{FS=OFS="\t"}} NR>2 && $3=="transcript" {{print $1, $4-1, $5, $6, $7 }}' {input} | \
        # extract the same transcript rows again for gene name information
        paste - <(awk 'BEGIN{{FS=OFS="\t"}}NR>2 && $3=="transcript"' {input} | \
            # if no gene name, use stringtie gene ID as gene name
            sed -e '/gene_name/!s/gene_id/gene_name/' | \
            # extract gene name from GFF tags
            grep -oP '(?<=gene_name \\").*?(?=\\";)') | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $6, $4, $5}}' | \
        # cluster annotations that overlap by cluster dist
        sort -k1,1 -k2,2n | bedtools cluster -i stdin -s -d {params.cluster_dist} | \
        # keep the highest scoring annotation from each cluster
        sort -k7,7n -k5,5n | bedtools groupby -g 7 -c 1,2,3,4,5,6 -o first | cut --complement -f1 > {output}
        """












