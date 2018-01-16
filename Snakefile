#!/usr/bin/env python
import os
from math import log2
from math import log10

configfile: "config.yaml"

SAMPLES = config["samples"]
sisamples = {k:v for (k,v) in SAMPLES.items() if v["spikein"]=="y"}
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
sipassing = {k:v for (k,v) in PASSING.items() if v["spikein"] == "y"}

controlgroups = config["comparisons"]["libsizenorm"]["controls"]
conditiongroups = config["comparisons"]["libsizenorm"]["conditions"]

if sipassing:
    controlgroups_si = config["comparisons"]["spikenorm"]["controls"]
    conditiongroups_si = config["comparisons"]["spikenorm"]["conditions"]

# CATEGORIES = ["genic", "intragenic", "intergenic", "antisense", "convergent", "divergent"]
COUNTTYPES = ["counts", "sicounts"] if sisamples else ["counts"]
NORMS = ["libsizenorm", "spikenorm"] if sisamples else ["libsizenorm"]

localrules:
    all,
    fastqc_aggregate,
    get_si_pct, plot_si_pct,
    make_stranded_genome, make_stranded_annotations,
    cat_matrices

rule all:
    input:
        #FastQC
        'qual_ctrl/fastqc/per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.{fmt}", norm=NORMS+COUNTTYPES, sample=SAMPLES, readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE", "plus", "minus"], fmt=["bedgraph", "bw"]),
        #quality control
        "qual_ctrl/read_processing-loss.svg",
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]) if sisamples else [],
        # expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all", "passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-netseq-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"], windowsize=config["corr-windowsizes"]) + expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-netseq-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status = ["all", "passing"], windowsize=config["corr-windowsizes"]) if sisamples else
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-netseq-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status = ["all", "passing"], windowsize=config["corr-windowsizes"]),
        #datavis
        # expand("datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{readtype}-{strand}.tsv.gz", annotation=config["annotations"], norm=NORMS, readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE"]),
        expand(expand("datavis/{{annotation}}/spikenorm/netseq-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-{{strand}}-{{plottype}}-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], status=["all","passing"], readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE"], plottype=["heatmap", "metagene"]) +
        expand(expand("datavis/{{annotation}}/libsizenorm/netseq-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-{{strand}}-{{plottype}}-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], status=["all","passing"], readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE"], plottype=["heatmap", "metagene"]) if sisamples else
        expand(expand("datavis/{{annotation}}/libsizenorm/netseq-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-{{strand}}-{{plottype}}-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], status=["all","passing"], readtype=["5end", "wholeread"], strand=["SENSE", "ANTISENSE"], plottype=["heatmap", "metagene"]),

def plotcorrsamples(wildcards):
    dd = SAMPLES if wildcards.status=="all" else PASSING
    if wildcards.condition=="all":
        if wildcards.norm=="libsizenorm": #condition==all,norm==lib
            return list(dd.keys())
        else: #condition==all,norm==spike
            return [k for k,v in dd.items() if v["spikein"]=="y"]
    elif wildcards.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in dd.items() if v["group"]==wildcards.control or v["group"]==wildcards.condition]
    else: #condition!=all;norm==spike
        return [k for k,v in dd.items() if (v["group"]==wildcards.control or v["group"]==wildcards.condition) and v["spikein"]=="y"]

rule fastqc_raw:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc-raw-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/raw/{wildcards.sample}) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw/{wildcards.sample} {input}) &>> {log}
        """

rule clean_reads:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        fq = temp("fastq/cleaned/{sample}-trim.fastq"),
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"]
    shell: """
        (cutadapt --cut=-1 --trim-n -a {params.adapter} --nextseq-trim={params.trim_qual} -m 12 --length-tag 'length=' -o {output.fq} {input}) &> {output.log}
        """

rule remove_molec_barcode:
    input:
        "fastq/cleaned/{sample}-trim.fastq"
    output:
        fq = "fastq/cleaned/{sample}-clean.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molec_barcode/removeMBC-{sample}.log"
    shell: """
        (python scripts/extractMolecularBarcode.py {input} fastq/cleaned/{wildcards.sample}-clean.fastq {output.barcodes} {output.ligation}) &> {log}
        (pigz -f fastq/cleaned/{wildcards.sample}-clean.fastq) &>> {log}
        """

rule fastqc_cleaned:
    input:
        "fastq/cleaned/{sample}-clean.fastq.gz"
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/cleaned/{sample}-clean_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/cleaned/fastqc-cleaned-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/cleaned) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/cleaned {input}) &>> {log}
        """

rule fastqc_aligned:
    input:
        lambda wildcards: "alignment/" + wildcards.sample + "-noPCRdup.bam" if wildcards.fqtype=="aligned_noPCRdup" else "alignment/" + wildcards.sample + "/unmapped.bam",
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}-{fqtype}_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/{fqtype}/fastqc-{fqtype}-{sample}.log"
    wildcard_constraints:
        fqtype="aligned_noPCRdup|unaligned"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}) &> {log}
        (bedtools bamtofastq -fq qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}-{wildcards.fqtype}.fastq -i {input}) &>> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}-{wildcards.fqtype}.fastq) &>> {log}
        (rm qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}-{wildcards.fqtype}.fastq) &>> {log}
        """

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt", zip, sample=SAMPLES, fname=[os.path.split(v["fastq"])[1].split(".fastq")[0] + "_fastqc" for k,v in SAMPLES.items()]),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}-clean_fastqc/fastqc_data.txt", sample=SAMPLES),
        aligned_noPCRdup = expand("qual_ctrl/fastqc/aligned_noPCRdup/{sample}-aligned_noPCRdup_fastqc/fastqc_data.txt", sample=SAMPLES),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}-unaligned_fastqc/fastqc_data.txt", sample=SAMPLES),
    output:
        'qual_ctrl/fastqc/per_base_quality.tsv',
        'qual_ctrl/fastqc/per_tile_quality.tsv',
        'qual_ctrl/fastqc/per_sequence_quality.tsv',
        'qual_ctrl/fastqc/per_base_sequence_content.tsv',
        'qual_ctrl/fastqc/per_sequence_gc.tsv',
        'qual_ctrl/fastqc/per_base_n.tsv',
        'qual_ctrl/fastqc/sequence_length_distribution.tsv',
        'qual_ctrl/fastqc/sequence_duplication_levels.tsv',
        'qual_ctrl/fastqc/adapter_content.tsv',
        'qual_ctrl/fastqc/kmer_content.tsv'
    run:
        shell("rm -f {output}")
        #for each statistic
        for outpath, stat, header in zip(output, ["Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Total Deduplicated Percentage", "Adapter Content", "Kmer Content"], ["base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus", "tile\tbase\tmean\tsample\tstatus",
        "quality\tcount\tsample\tstatus", "base\tg\ta\tt\tc\tsample\tstatus", "gc_content\tcount\tsample\tstatus", "base\tn_count\tsample\tstatus", "length\tcount\tsample\tstatus", "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus", "position\tpct\tsample\tstatus",
        "sequence\tcount\tpval\tobs_over_exp_max\tmax_position\tsample\tstatus" ]):
            for input_type in ["raw", "cleaned", "aligned_noPCRdup", "unaligned"]:
                for sample_id, fqc in zip(SAMPLES.keys(), input[input_type]):
                    shell("""awk 'BEGIN{{FS=OFS="\t"}} /{stat}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, "{sample_id}", "{input_type}"}}' {fqc} | tail -n +2 >> {outpath}""")
            shell("""sed -i "1i {header}" {outpath}""")

rule plot_fastqc_summary:
    input:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.tsv',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.tsv',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.tsv',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.tsv',
        per_base_n = 'qual_ctrl/fastqc/per_base_n.tsv',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.tsv',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.tsv',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.tsv',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.tsv',
        kmer = 'qual_ctrl/fastqc/kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.svg',
        kmer = 'qual_ctrl/fastqc/kmer_content.svg',
    script: "scripts/fastqc_summary.R"

#align to combined genome with Tophat2 (single genome only if no samples have spike-ins), WITHOUT reference transcriptome (i.e., the -G gff)
#(because we don't always have a reference gff and it doesn't make much difference)
rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"] if sisamples else config["genome"]["fasta"]
    output:
        expand(config["tophat2"]["bowtie2-index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["tophat2"]["bowtie2-index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["tophat2"]["bowtie2-index-path"],
        prefix = config["combinedgenome"]["experimental_prefix"]
    log: "logs/bowtie2_build.log"
    run:
        if sisamples:
            shell("(bowtie2-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}")
        else:
            shell("(sed -e 's/>/>{params.prefix}/g' {input.fasta} > .{params.prefix}.fa; bowtie2-build .{params.prefix}.fa {params.idx_path}/{wildcards.basename}; rm .{params.prefix}.fa) &> {log}")

rule align:
    input:
        expand(config["tophat2"]["bowtie2-index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.bt2", num = [1,2,3,4]) if sisamples else expand(config["tophat2"]["bowtie2-index-path"] + "/" + config["genome"]["name"] + ".{num}.bt2", num = [1,2,3,4]),
        expand(config["tophat2"]["bowtie2-index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.bt2", num=[1,2]) if sisamples else expand(config["tophat2"]["bowtie2-index-path"] + "/" + config["genome"]["name"] + ".rev.{num}.bt2", num=[1,2]),
        fastq = "fastq/cleaned/{sample}-clean.fastq.gz"
    output:
        aligned = "alignment/{sample}/accepted_hits.bam",
        unaligned = "alignment/{sample}/unmapped.bam",
        summary = "alignment/{sample}/align_summary.txt",
    params:
        idx_path = config["tophat2"]["bowtie2-index-path"],
        basename = config["combinedgenome"]["name"] if sisamples else config["genome"]["name"],
        read_mismatches = config["tophat2"]["read-mismatches"],
        read_gap_length = config["tophat2"]["read-gap-length"],
        read_edit_dist = config["tophat2"]["read-edit-dist"],
        min_anchor_length = config["tophat2"]["min-anchor-length"],
        splice_mismatches = config["tophat2"]["splice-mismatches"],
        min_intron_length = config["tophat2"]["min-intron-length"],
        max_intron_length = config["tophat2"]["max-intron-length"],
        max_insertion_length = config["tophat2"]["max-insertion-length"],
        max_deletion_length = config["tophat2"]["max-deletion-length"],
        max_multihits = config["tophat2"]["max-multihits"],
        segment_mismatches = config["tophat2"]["segment-mismatches"],
        segment_length = config["tophat2"]["segment-length"],
        min_coverage_intron = config["tophat2"]["min-coverage-intron"],
        max_coverage_intron = config["tophat2"]["max-coverage-intron"],
        min_segment_intron = config["tophat2"]["min-segment-intron"],
        max_segment_intron = config["tophat2"]["max-segment-intron"],
    conda:
        "envs/tophat2.yaml"
    threads : config["threads"]
    log: "logs/align/align-{sample}.log"
    shell:
        """
        (tophat2 --read-mismatches {params.read_mismatches} --read-gap-length {params.read_gap_length} --read-edit-dist {params.read_edit_dist} -o alignment/{wildcards.sample} --min-anchor-length {params.min_anchor_length} --splice-mismatches {params.splice_mismatches} --min-intron-length {params.min_intron_length} --max-intron-length {params.max_intron_length} --max-insertion-length {params.max_insertion_length} --max-deletion-length {params.max_deletion_length} --num-threads {threads} --max-multihits {params.max_multihits} --library-type fr-firststrand --segment-mismatches {params.segment_mismatches} --no-coverage-search --segment-length {params.segment_length} --min-coverage-intron {params.min_coverage_intron} --max-coverage-intron {params.max_coverage_intron} --min-segment-intron {params.min_segment_intron} --max-segment-intron {params.max_segment_intron} --b2-sensitive {params.idx_path}/{params.basename} {input.fastq}) &> {log}
        """

rule select_unique_mappers:
    input:
        "alignment/{sample}/accepted_hits.bam"
    output:
        temp("alignment/{sample}-unique.bam")
    threads: config["threads"]
    log: "logs/select_unique_mappers/select_unique_mappers-{sample}.log"
    shell: """
        (samtools view -buh -q 50 -@ {threads} {input} | samtools sort -T .{wildcards.sample} -@ {threads} - > {output}) &> {log}
        """

rule remove_PCR_duplicates:
    input:
        "alignment/{sample}-unique.bam"
    output:
        "alignment/{sample}-noPCRdup.bam"
    log: "logs/remove_PCR_duplicates/removePCRduplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        """

rule read_processing_numbers:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("alignment/{sample}/align_summary.txt", sample=SAMPLES),
        nodups = expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES)
    output:
        "qual_ctrl/read_processing_summary.tsv"
    log: "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map\tnoPCRdup" > {output}) &> {log}""")
        for sample, adapter, align, nodups in zip(SAMPLES.keys(), input.adapter, input.align, input.nodups):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(awk 'BEGIN{{ORS="\t"}} NR==3 || NR==4{{print $3}}' {align} >> {output}) &> {log}""")
            shell("""(samtools view -c {nodups} | awk '{{print $1}}' >> {output}) &> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$6=$5-$6; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing_summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing-loss.svg",
    script: "scripts/processing_summary.R"

rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam"
    params:
        prefix = lambda wildcards: config["combinedgenome"]["experimental_prefix"] if wildcards.counttype=="counts" else config["combinedgenome"]["spikein_prefix"]
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
        plmin = lambda wildcards: "coverage/counts/" + wildcards.sample + "-netseq-counts-5end-plmin.bedgraph" if wildcards.norm=="libsizenorm" else "coverage/sicounts/" + wildcards.sample + "-netseq-counts-5end-plmin.bedgraph"
    params:
        scalefactor = lambda wildcards: config["spikein-pct"] if wildcards.norm=="spikenorm" else 1
    output:
        normalized = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus"
    log: "logs/normalize/normalize-{sample}-{norm}-{readtype}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.plmin} {input.counts} {params.scalefactor} > {output.normalized}) &> {log}
        """

rule get_si_pct:
    input:
        plmin = expand("coverage/counts/{sample}-tss-counts-plmin.bedgraph", sample=sisamples),
        SIplmin = expand("coverage/sicounts/{sample}-tss-sicounts-plmin.bedgraph", sample=sisamples)
    params:
        group = [v["group"] for k,v in sisamples.items()]
    output:
        "qual_ctrl/all/spikein-counts.tsv"
    log: "logs/get_si_pct.log"
    run:
        shell("rm -f {output}")
        for name, exp, si, g in zip(sisamples.keys(), input.plmin, input.SIplmin, params.group):
            shell("""(echo -e "{name}\t{g}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {si} {exp}) >> {output}) &> {log}""")

rule plot_si_pct:
    input:
        "qual_ctrl/all/spikein-counts.tsv"
    output:
        plot = "qual_ctrl/{status}/{status}-spikein-plots.svg",
        stats = "qual_ctrl/{status}/{status}-spikein-stats.tsv"
    params:
        samplelist = lambda wildcards : [k for k,v in sisamples.items() if v["spikein"]=="y"] if wildcards.status=="all" else [k for k,v in sipassing.items() if v["spikein"]=="y"],
        conditions = config["comparisons"]["spikenorm"]["conditions"],
        controls = config["comparisons"]["spikenorm"]["controls"],
    script: "scripts/plotsipct.R"

rule make_stranded_genome:
    input:
        exp = config["genome"]["chrsizes"],
        si = config["genome"]["sichrsizes"]
    output:
        exp = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
        si = os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv",
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.exp} > {output.exp}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.si} > {output.si}) &>> {log}
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

rule map_to_windows:
    input:
        bg = "coverage/{norm}/{sample}-netseq-{norm}-5end-SENSE.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    output:
        exp = temp("coverage/{norm}/{sample}-netseq-window-{windowsize}-coverage-{norm}.bedgraph"),
    shell: """
        bedtools makewindows -g {input.chrsizes} -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output.exp}
        """

rule join_window_counts:
    input:
        exp = expand("coverage/{{norm}}/{sample}-netseq-window-{{windowsize}}-coverage-{{norm}}.bedgraph", sample=SAMPLES),
    output:
        exp = "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    params:
        names = list(SAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output.exp}) &> {log}
        """

rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    output:
        "qual_ctrl/{status}/{condition}-v-{control}-netseq-{status}-window-{windowsize}-{norm}-correlations.svg"
    params:
        pcount = lambda wildcards: 0.1*int(wildcards.windowsize),
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "{annopath}/stranded/{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

def selectchrom(wildcards):
    if wildcards.strand in ["plus", "minus"]:
        if wildcards.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wildcards.norm=="sicounts":
        return os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv"
    return os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"

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

rule deeptools_matrix:
    input:
        annotation = lambda wildcards: os.path.dirname(config["annotations"][wildcards.annotation]["path"]) + "/stranded/" + wildcards.annotation + "-STRANDED" + os.path.splitext(config["annotations"][wildcards.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bw",
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{readtype}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{readtype}-{strand}.tsv"),
        matrix_gz = "datavis/{annotation}/{norm}/{annotation}-{sample}-{norm}-{readtype}-{strand}.tsv.gz"

rule compute_matrix:
    input:
        annotation = lambda wildcards: os.path.dirname(config["annotations"][wildcards.annotation]["path"]) + "/stranded/" + wildcards.annotation + "-STRANDED" + os.path.splitext(config["annotations"][wildcards.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}-netseq-{norm}-{readtype}-{strand}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}.tsv"),
        melted = "datavis/{annotation}/{norm}/{annotation}_{sample}_{norm}-{readtype}-{strand}-melted.tsv.gz"
    params:
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"] + config["annotations"][wildcards.annotation]["binsize"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"] + config["annotations"][wildcards.annotation]["binsize"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"],
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{annotation}_{sample}_{norm}-{strand}.log"
    run:
        if config["annotations"][wildcards.annotation]["type"]=="absolute":
            refpoint = config["annotations"][wildcards.annotation]["refpoint"]
            if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
            else:
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            refpoint = "TSS"
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {refpoint} --group {params.group} -s {wildcards.sample} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/{{annotation}}_{sample}_{{norm}}-{{readtype}}-{{strand}}-melted.tsv.gz", sample=SAMPLES)
    output:
        "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{readtype}-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}-{norm}-{readtype}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_heatmaps:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{readtype}-{strand}.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-heatmap-bygroup.svg"
    params:
        samplelist = plotcorrsamples,
        mtype = lambda wildcards : config["annotations"][wildcards.annotation]["type"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        cluster = lambda wildcards : config["annotations"][wildcards.annotation]["cluster"],
        nclust = lambda wildcards: config["annotations"][wildcards.annotation]["nclusters"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    run:
        if config["annotations"][wildcards.annotation]["type"]=="scaled":
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            endlabel = config["annotations"][wildcards.annotation]["three_prime_label"]
        else:
            scaled_length=0
            endlabel = "HAIL SATAN!"
        shell("""Rscript scripts/plot_netseq_heatmaps.R -i {input.matrix} -s {params.samplelist} -t {params.mtype} -u {params.upstream} -d {params.dnstream} -c {params.pct_cutoff} -z {params.cluster} -k {params.nclust} -r {params.refpointlabel} -f {wildcards.strand} -l {scaled_length} -e {endlabel} -y {params.ylabel} -m {params.heatmap_cmap} -o {output.heatmap_sample} -p {output.heatmap_group}""")

rule plot_metagenes:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{norm}-{readtype}-{strand}.tsv.gz"
    output:
        meta_sample = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-metagene-bysample.svg",
        meta_sample_overlay = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-metagene-overlay-bysample.svg",
        meta_heatmap_sample = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-metaheatmap-bysample.svg",
        meta_group = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-metagene-bygroup.svg",
        meta_group_overlay = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-metagene-overlay-bygroup.svg",
        meta_heatmap_group = "datavis/{annotation}/{norm}/netseq-{annotation}-{norm}-{status}_{condition}-v-{control}_{readtype}-{strand}-metaheatmap-bygroup.svg",
    params:
        samplelist = plotcorrsamples,
        mtype = lambda wildcards : config["annotations"][wildcards.annotation]["type"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        trim_pct = lambda wildcards : config["annotations"][wildcards.annotation]["trim_pct"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    run:
        if config["annotations"][wildcards.annotation]["type"]=="scaled":
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            endlabel = config["annotations"][wildcards.annotation]["three_prime_label"]
        else:
            scaled_length=0
            endlabel = "HAIL SATAN!"
        shell("""Rscript scripts/plot_netseq_metagenes.R -i {input.matrix} -s {params.samplelist} -t {params.mtype} -f {wildcards.strand} -u {params.upstream} -d {params.dnstream} -c {params.trim_pct} -r {params.refpointlabel} -l {scaled_length} -e {endlabel} -y {params.ylabel} --out1 {output.meta_sample} --out2 {output.meta_sample_overlay} --out3 {output.meta_heatmap_sample} --out4 {output.meta_group} --out5 {output.meta_group_overlay} --out6 {output.meta_heatmap_group}""")

# rule map_counts_to_transcripts:
#     input:
#         bed = config["genome"]["transcripts"] if wildcards.type=="exp" else config["genome"]["spikein-transcripts"],
#         bg = lambda wildcards: "coverage/counts/" + wildcards.sample + "-netseq-counts-5end-SENSE.bedgraph" if wildcards.type=="exp" else "coverage/sicounts/" + wildcards.sample + "-netseq-sicounts-5end-SENSE.bedgraph"
#     output:
#         temp("diff_exp/{condition}-v-{control}/{sample}-{type}-allpeakcounts.tsv")
#     log: "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}-{sample}-{type}.log"
#     shell: """
#         (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum | awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-"$2"-"$3, $4}}' &> {output}) &> {log}
#         """

def getsamples(ctrl, cond):
    return [k for k,v in PASSING.items() if (v["group"]==ctrl or v["group"]==cond)]

# rule get_peak_counts:
#     input:
#         lambda wildcards : ["diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + x + "-" + wildcards.type + "-allpeakcounts.tsv" for x in getsamples(wildcards.control, wildcards.condition)]
#     output:
#         "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{type}-peak-counts.tsv"
#     params:
#         n = lambda wildcards: 2*len(getsamples(wildcards.control, wildcards.condition)),
#         names = lambda wildcards: "\t".join(getsamples(wildcards.control, wildcards.condition))
#     log: "logs/get_peak_counts/get_peak_counts-{condition}-v-{control}-{type}.log"
#     shell: """
#         (paste {input} | cut -f$(paste -d, <(echo "1") <(seq -s, 2 2 {params.n})) | cat <(echo -e "name\t" "{params.names}" ) - > {output}) &> {log}
#         """

# rule call_de_peaks:
#     input:
#         expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-exp-peak-counts.tsv",
#         sicounts = lambda wildcards: "diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-si-peak-counts.tsv" if wildcards.norm=="spikenorm" else "diff_exp/" + wildcards.condition + "-v-" + wildcards.control + "/" + wildcards.condition + "-v-" + wildcards.control + "-exp-peak-counts.tsv"
#     params:
#         samples = lambda wildcards : getsamples(wildcards.control, wildcards.condition),
#         groups = lambda wildcards : [PASSING[x]["group"] for x in getsamples(wildcards.control, wildcards.condition)],
#         alpha = config["deseq"]["fdr"],
#         lfc = log2(config["deseq"]["fold-change-threshold"])
#     output:
#         results = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
#         #need to write out norm counts here or just in the total qc?
#         normcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-counts-sfnorm-{norm}.tsv",
#         rldcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-counts-rlog-{norm}.tsv",
#         qcplots = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-qcplots-{norm}.svg"
#     script:
#         "scripts/call_de_peaks.R"

# rule separate_de_peaks:
#     input:
#         "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
#     output:
#         up = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-up.tsv",
#         down = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-down.tsv",
#     params:
#         fdr = -log10(config["deseq"]["fdr"])
#     log: "logs/separate_de_peaks/separate_de_peaks-{condition}-v-{control}-{norm}.log"
#     shell: """
#         (awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}} NR==1{{print > "{output.up}"; print > "{output.down}" }} NR>1 && $10>afdr && $7>0 {{print > "{output.up}"}} NR>1 && $10>afdr && $7<0 {{print > "{output.down}"}}' {input}) &> {log}
#         """

# rule de_peaks_to_bed:
#     input:
#         "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.tsv",
#     output:
#         "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-{direction}.bed",
#     log: "logs/de_peaks_to_bed/de_peaks_to_bed-{condition}-v-{control}-{norm}-{direction}.log"
#     shell: """
#         (tail -n +2 {input} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, $7":"$11, $3}}' > {output}) &> {log}
#         """

# rule summarise_de_results:
#     input:
#         total = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-results-{norm}-all.tsv",
#         genic = "diff_exp/{condition}-v-{control}/genic/{condition}-v-{control}-results-{norm}-all-genic.tsv",
#         intragenic = "diff_exp/{condition}-v-{control}/intragenic/{condition}-v-{control}-results-{norm}-all-intragenic.tsv",
#         antisense = "diff_exp/{condition}-v-{control}/antisense/{condition}-v-{control}-results-{norm}-all-antisense.tsv",
#         convergent = "diff_exp/{condition}-v-{control}/convergent/{condition}-v-{control}-results-{norm}-all-convergent.tsv",
#         divergent = "diff_exp/{condition}-v-{control}/divergent/{condition}-v-{control}-results-{norm}-all-divergent.tsv",
#         intergenic = "diff_exp/{condition}-v-{control}/intergenic/{condition}-v-{control}-results-{norm}-all-intergenic.tsv",
#     params:
#         lfc = log2(config["deseq"]["fold-change-threshold"]),
#         alpha = config["deseq"]["fdr"]
#     output:
#         summary = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{norm}-diffexp-summary.svg",
#         maplot = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{norm}-diffexp-maplot.svg",
#         volcano = "diff_exp/{condition}-v-{control}/{condition}-v-{control}-{norm}-diffexp-volcano.svg",
#     script: "scripts/de_summary.R"
