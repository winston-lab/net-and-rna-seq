#!/usr/bin/env python

rule read_processing_numbers:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("alignment/{sample}/align_summary.txt", sample=SAMPLES),
        nodups = expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES) if config["random-hexamer"] else expand("alignment/{sample}-unique.bam", sample=SAMPLES)
    output:
        "qual_ctrl/read_processing_summary.tsv"
    log: "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map\tnoPCRdup" > {output}) &> {log}""")
        for sample, adapter, align, nodups in zip(SAMPLES.keys(), input.adapter, input.align, input.nodups):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(awk 'BEGIN{{ORS="\t"}} NR==3 || NR==4{{print $3}}' {align} >> {output}) &> {log}""")
            shell("""(samtools view -c {nodups} | awk '{{print $1}}' >> {output}) &> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$5=$4-$5; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing_summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing-loss.svg",
    script: "scripts/processing_summary.R"

rule get_si_pct:
    input:
        plmin = expand("coverage/counts/{sample}-netseq-counts-5end-plmin.bedgraph", sample=SISAMPLES),
        SIplmin = expand("coverage/sicounts/{sample}-netseq-sicounts-5end-plmin.bedgraph", sample=SISAMPLES)
    params:
        group = [v["group"] for k,v in SISAMPLES.items()]
    output:
        "qual_ctrl/all/spikein-counts.tsv"
    log: "logs/get_si_pct.log"
    run:
        shell("rm -f {output}")
        for name, exp, si, g in zip(SISAMPLES.keys(), input.plmin, input.SIplmin, params.group):
            shell("""(echo -e "{name}\t{g}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {si} {exp}) >> {output}) &> {log}""")

rule plot_si_pct:
    input:
        "qual_ctrl/all/spikein-counts.tsv"
    output:
        plot = "qual_ctrl/{status}/{status}-spikein-plots.svg",
        stats = "qual_ctrl/{status}/{status}-spikein-stats.tsv"
    params:
        samplelist = lambda wc : list(SISAMPLES.keys()) if wc.status=="all" else list(SIPASSING.keys()),
        conditions = conditiongroups_si if SISAMPLES else [],
        controls = controlgroups_si if SISAMPLES else [],
    script: "scripts/plotsipct.R"

