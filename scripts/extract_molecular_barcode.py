#!/usr/bin/env python

"""
Date : April 24th, 2014
Last updated: June 27th, 2018

extract molecular barcode from reads and write out barcode and ligation counts

use: python extract_molecular_barcode.py in_fastq out_fastq out_barcode out_ligation
"""

def main():
    import sys, itertools, subprocess
    from collections import Counter

    nucleotides = 'ACTGN'
    barcodes_dict = Counter({"".join(nts) : 0 for nts in itertools.product(nucleotides, repeat=6)})
    ligations_dict = Counter({"".join(nts) : 0 for nts in itertools.product(nucleotides, repeat=6)})

    with open(sys.argv[2], 'w') as out_fastq_file:
        in_fastq= subprocess.Popen(['pigz', '-cdfq', sys.argv[1]], stdout=subprocess.PIPE, encoding='utf-8')
        out_fastq= subprocess.Popen(['pigz', '-fq', '-'], stdin=subprocess.PIPE, stdout=out_fastq_file, encoding='utf-8')
        header = in_fastq.stdout.readline().rstrip()
        while header != '':
            full_seq = in_fastq.stdout.readline()
            plus_line = in_fastq.stdout.readline()
            qual_scores = in_fastq.stdout.readline()
            barcode = full_seq[0:6]
            ligation = full_seq[3:9]
            trimmed_seq = full_seq[6:]
            out_fastq.stdin.write(header.split(" ")[0]+'_MolecularBarcode:'+barcode+' '+header.split(" ")[1]+'\n')
            out_fastq.stdin.write(trimmed_seq)
            out_fastq.stdin.write(plus_line)
            out_fastq.stdin.write(qual_scores[6:])
            header = in_fastq.stdout.readline().rstrip()

            barcodes_dict.update([barcode])
            if len(trimmed_seq) >= 4 :
                ligations_dict.update([ligation])

    with open(sys.argv[3], 'w') as out_barcode, open(sys.argv[4], 'w') as out_ligation:
        out_barcode.write("barcode\tcount\n")
        for barcode, count in barcodes_dict.items():
            out_barcode.write(f"{barcode}\t{count}\n")

        out_ligation.write("ligation\tcount\n")
        for ligation, count in ligations_dict.items():
            out_ligation.write(f"{ligation}\t{count}\n")

if __name__ == '__main__':
    main()

