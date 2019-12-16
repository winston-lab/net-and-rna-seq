#!/usr/bin/env python

"""
Date : April 24th, 2014
Last updated: Dec 16, 2019

extract molecular barcode from reads and write out barcode and ligation counts

use: python extract_molecular_barcode.py input.fastq.gz output.fastq.gz 6 barcodes.tsv ligation.tsv
"""

def main():
    import sys, itertools, subprocess
    from collections import Counter

    input_fastq_path, output_fastq_path, barcode_length, output_barcodes_path, output_ligation_path = sys.argv[1:6]
    barcode_length = int(barcode_length)

    nucleotides = 'ACTGN'
    barcodes_dict = Counter({"".join(nts) : 0 for nts in itertools.product(nucleotides, repeat=barcode_length)})
    ligations_dict = Counter({"".join(nts) : 0 for nts in itertools.product(nucleotides, repeat=barcode_length)})

    with open(output_fastq_path, 'w') as out_fastq_file:
        in_fastq= subprocess.Popen(['pigz', '-cdfq', input_fastq_path, stdout=subprocess.PIPE, encoding='utf-8')
        out_fastq= subprocess.Popen(['pigz', '-fq', '-'], stdin=subprocess.PIPE, stdout=out_fastq_file, encoding='utf-8')
        header = in_fastq.stdout.readline().rstrip()
        while header != '':
            full_seq = in_fastq.stdout.readline()
            plus_line = in_fastq.stdout.readline()
            qual_scores = in_fastq.stdout.readline()
            barcode = full_seq[0:barcode_length]
            ligation = full_seq[barcode_length-3 : barcode_length+3]
            trimmed_seq = full_seq[barcode_length:]
            out_fastq.stdin.write(header.split(" ")[0]+'_MolecularBarcode:'+barcode+' '+header.split(" ")[1]+'\n')
            out_fastq.stdin.write(trimmed_seq)
            out_fastq.stdin.write(plus_line)
            out_fastq.stdin.write(qual_scores[barcode_length:])
            header = in_fastq.stdout.readline().rstrip()

            if len(barcode) == barcode_length:
                barcodes_dict.update([barcode])
            if len(trimmed_seq) >= 4 :
                ligations_dict.update([ligation])

    with open(output_barcodes_path, 'w') as out_barcode, open(output_ligation_path, 'w') as out_ligation:
        out_barcode.write("barcode\tcount\n")
        for barcode, count in barcodes_dict.items():
            out_barcode.write(f"{barcode}\t{count}\n")

        out_ligation.write("ligation\tcount\n")
        for ligation, count in ligations_dict.items():
            out_ligation.write(f"{ligation}\t{count}\n")

if __name__ == '__main__':
    main()

