#!/usr/bin/env python

localrules:
    classify_genic_diffexp_transcript_annotation,
    classify_antisense_diffexp_transcript_annotation,
    classify_convergent_diffexp_transcript_annotation,
    classify_divergent_diffexp_transcript_annotation,
    classify_intergenic_diffexp_transcript_annotation,

# genic transcript_annotation are those which were in
# the original annotation and thus aren't named "stringtie.xx"
rule classify_genic_diffexp_transcript_annotation:
    input:
        "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{direction}.tsv",
    output:
        table = "diff_exp/transcripts/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-genic-{direction}.tsv",
        bed = "diff_exp/transcripts/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-genic-{direction}.bed",
    log:
        "logs/classify_genic_diffexp_transcript_annotation/classify_genic_diffexp_transcript_annotation_{condition}-v-{control}_{norm}-{direction}-{assay}.log"
    shell: """
        (awk 'NR==1 || ($4 !~ /stringtie./)' {input} | \
         tee {output.table} | \
         tail -n +2 | \
         cut -f1-6 > {output.bed}) &> {log}
        """

# 0. cut header off results
# 1. get "TSS" information from start of transcript
# 2. get "TSSs" that are antisense to transcript_annotation
# 3. add distance of sense to antisense "TSS"
# 4. remove "TSS" information
# 5. add new header
rule classify_antisense_diffexp_transcript_annotation:
    input:
        results = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{direction}.tsv",
        transcript_annotation = os.path.abspath(build_annotations(config["genome"]["transcript_annotation"]))
    output:
        table = "diff_exp/transcripts/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-antisense-{direction}.tsv",
        bed = "diff_exp/transcripts/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-antisense-{direction}.bed",
    log:
        "logs/classify_antisense_diffexp_transcript_annotation/classify_antisense_diffexp_transcript_annotation_{condition}-v-{control}_{norm}-{direction}-{assay}.log"
    shell: """
        (tail -n +2 {input.results} | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6, $0}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6, $0}}' | \
         bedtools intersect -wo -S -a stdin -b <(cut -f1-6 {input.transcript_annotation}) | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{$27=$23-$3}} $6=="-"{{$27=$2-$22}} {{print $0}}' | \
         cut --complement -f1-6 | \
         cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_anti_tss_dist")) - | \
         tee {output.table} | \
         tail -n +2 | \
         cut -f1-6 > {output.bed}) &> {log}
        """

# 0. cut header off results
# 1. get "TSS" information from start of transcript
# 2. exclude genic "TSSs"
# 3. get "TSSs" in convergent regions
# 4. join to information of transcript that the convergent transcript is convergent to (rather than the convergent region info)
# 5. re-sort by significance
# 6. add distance of sense to convergent "TSS"
# 4. remove "TSS" information
# 5. add new header
rule classify_convergent_diffexp_transcript_annotation:
    input:
        transcript_anno = os.path.abspath(build_annotations(config["genome"]["transcript_annotation"])),
        conv_anno = build_annotations("annotations/" + config["genome"]["name"] + "_convergent-regions.bed"),
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        results = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{direction}.tsv",
    output:
        table = "diff_exp/transcripts/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-convergent-{direction}.tsv",
        bed = "diff_exp/transcripts/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-convergent-{direction}.bed",
    log:
        "logs/classify_convergent_diffexp_transcript_annotation/classify_convergent_diffexp_transcript_annotation_{condition}-v-{control}_{norm}-{direction}-{assay}.log"
    shell: """
        (tail -n +2 {input.results} | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6, $0}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6, $0}}' | \
         bedtools intersect -a stdin -b {input.genic_anno} -v -s | \
         bedtools intersect -wo -s -a stdin -b {input.conv_anno} | \
         sort -k24,24 | \
         join -1 24 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4,4 {input.transcript_anno}) | \
         sort -k17,17nr -k16,16nr | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{$27=$23-$3}} $6=="-"{{$27=$2-$22}} {{print $0}}' | \
         cut --complement -f1-6 | \
         cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_conv_tss_dist")) - | \
         tee {output.table} | \
         tail -n +2 | \
         cut -f1-6 > {output.bed}) &> {log}
        """

# 0. cut header off results
# 1. get "TSS" information from start of transcript
# 2. exclude genic "TSSs"
# 3. get "TSSs" in divergent regions
# 4. join to information of transcript that the divergent transcript is divergent to (rather than the divergent region info)
# 5. re-sort by significance
# 6. add distance of sense to divergent "TSS"
# 4. remove "TSS" information
# 5. add new header
rule classify_divergent_diffexp_transcript_annotation:
    input:
        transcript_anno = os.path.abspath(build_annotations(config["genome"]["transcript_annotation"])),
        div_anno = build_annotations("annotations/" + config["genome"]["name"] + "_divergent-regions.bed"),
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        results = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{direction}.tsv",
    output:
        table = "diff_exp/transcripts/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-divergent-{direction}.tsv",
        bed = "diff_exp/transcripts/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-divergent-{direction}.bed",
    log:
        "logs/classify_divergent_diffexp_transcript_annotation/classify_divergent_diffexp_transcript_annotation_{condition}-v-{control}_{norm}-{direction}-{assay}.log"
    shell: """
        (tail -n +2 {input.results} | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6, $0}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6, $0}}' | \
         bedtools intersect -a stdin -b {input.genic_anno} -v -s | \
         bedtools intersect -wo -s -a stdin -b {input.div_anno} | \
         sort -k24,24 | \
         join -1 24 -2 4 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,2.1,2.2,2.3,2.4,2.5,2.6 - <(sort -k4,4 {input.transcript_anno}) | \
         sort -k17,17nr -k16,16nr | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{$27=$3-$23}} $6=="-"{{$27=$22-$2}} {{print $0}}' | \
         cut --complement -f1-6 | \
         cat <(paste <(head -n 1 {input.results}) <(echo -e "transcript_chrom\ttranscript_start\ttranscript_end\ttranscript_name\ttranscript_score\ttranscript_strand\tsense_tss_to_div_tss_dist")) - | \
         tee {output.table} | \
         tail -n +2 | \
         cut -f1-6 > {output.bed}) &> {log}
        """

# 0. exclude transcript_annotation in reference annotation
# 1. get "TSS" information from start of transcript
# 2. exclude transcript_annotation starting in transcript or genic regions
# 3. get "TSSs" in intergenic regions
# 4. remove "TSS" information
# 5. add back header
rule classify_intergenic_diffexp_transcript_annotation:
    input:
        intergenic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_intergenic-regions.bed"),
        transcript_anno = os.path.abspath(build_annotations(config["genome"]["transcript_annotation"])),
        genic_anno = build_annotations("annotations/" + config["genome"]["name"] + "_genic-regions.bed"),
        results = "diff_exp/transcripts/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-{direction}.tsv",
    output:
        table = "diff_exp/transcripts/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-intergenic-{direction}.tsv",
        bed = "diff_exp/transcripts/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{assay}-{norm}-transcripts-diffexp-results-intergenic-{direction}.bed",
    log:
        "logs/classify_intergenic_diffexp_transcript_annotation/classify_intergenic_diffexp_transcript_annotation_{condition}-v-{control}_{norm}-{direction}-{assay}.log"
    shell: """
        (awk '$4 ~ /stringtie./' {input.results} | \
         awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $2+1, $4, $5, $6, $0}} $6=="-"{{print $1, $3-1, $3, $4, $5, $6, $0}}' | \
         bedtools intersect -a stdin -b {input.transcript_anno} {input.genic_anno} -v | \
         bedtools intersect -u -a stdin -b {input.intergenic_anno} | \
         cut --complement -f1-6 | \
         cat <(head -n 1 {input.results}) - | \
         tee {output.table} | \
         tail -n +2 | \
         cut -f1-6 > {output.bed}) &> {log}
        """

