#!/usr/bin/env python

rule fastqc_raw:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
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

rule fastqc_cleaned:
    input:
        lambda wc: "fastq/cleaned/" + wc.sample + "-clean.fastq.gz" if wc.fqtype=="clean" else "fastq/cleaned/" + wc.sample + "-trim.fastq"
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/cleaned/{sample}-{fqtype}_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/cleaned/fastqc-cleaned-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/cleaned) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/cleaned {input}) &>> {log}
        """

rule fastqc_aligned:
    input:
        lambda wc: "alignment/" + wc.sample + "-noPCRdup.bam" if wc.fqtype=="aligned_noPCRdup" else "alignment/" + wc.sample + "-unique.bam" if wc.fqtype=="aligned" else "alignment/" + wc.sample + "/unmapped.bam",
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}-{fqtype}_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/{fqtype}/fastqc-{fqtype}-{sample}.log"
    wildcard_constraints:
        fqtype="aligned|aligned_noPCRdup|unaligned"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}) &> {log}
        (bedtools bamtofastq -fq qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}-{wildcards.fqtype}.fastq -i {input}) &>> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}-{wildcards.fqtype}.fastq) &>> {log}
        (rm qual_ctrl/fastqc/{wildcards.fqtype}/{wildcards.sample}-{wildcards.fqtype}.fastq) &>> {log}
        """

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt", zip, sample=SAMPLES, fname=[os.path.split(v["fastq"])[1].split(".fastq")[0] + "_fastqc" for k,v in SAMPLES.items()]),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}-clean_fastqc/fastqc_data.txt", sample=SAMPLES) if config["random-hexamer"] else expand("qual_ctrl/fastqc/cleaned/{sample}-trim_fastqc/fastqc_data.txt", sample=SAMPLES),
        aligned = expand("qual_ctrl/fastqc/aligned_noPCRdup/{sample}-aligned_noPCRdup_fastqc/fastqc_data.txt", sample=SAMPLES) if config["random-hexamer"] else expand("qual_ctrl/fastqc/aligned/{sample}-aligned_fastqc/fastqc_data.txt", sample=SAMPLES),
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
        # 'qual_ctrl/fastqc/kmer_content.tsv'
    run:
        shell("rm -f {output}")
        aln = "aligned_noPCRdup" if config["random-hexamer"] else "aligned"
        #for each statistic
        for outpath, stat, header in zip(output, ["Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Total Deduplicated Percentage", "Adapter Content"], ["base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus", "tile\tbase\tmean\tsample\tstatus",
        "quality\tcount\tsample\tstatus", "base\tg\ta\tt\tc\tsample\tstatus", "gc_content\tcount\tsample\tstatus", "base\tn_count\tsample\tstatus", "length\tcount\tsample\tstatus", "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus", "position\tpct\tsample\tstatus"]):
            for input_type in ["raw", "cleaned", aln, "unaligned"]:
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
        # kmer = 'qual_ctrl/fastqc/kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.svg',
        # kmer = 'qual_ctrl/fastqc/kmer_content.svg',
    script: "scripts/fastqc_summary.R"

