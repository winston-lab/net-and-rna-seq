#!/usr/bin/env python

"""
Date : April 24th, 2014
Last updated: June 26th, 2018

extract molecular barcode from reads and write out barcode and ligation counts

use: python extract_molecular_barcode.py in_fastq out_fastq out_barcode out_ligation
"""

def main:
    import sys, itertools, gzip
    from collections import Counter

    nucleotides = 'ACTGN'
    barcodes_dict = Counter({"".join(nts) : 0 for nts in itertools.product(nucleotides, repeat=6)})
    ligations_dict = Counter({"".join(nts) : 0 for nts in itertools.product(nucleotides, repeat=6)})

    with gzip.open(sys.argv[1], 'rt') as in_fastq, gzip.open(sys.argv[2], 'wt') as out_fastq, open(sys.argv[3], 'w') as out_barcode, open(sys.argv[4], 'w') as out_ligation:

        header = in_fastq.readline().rstrip()
        while header != '':
            full_seq = in_fastq.readline()
            plus_line = in_fastq.readline()
            qual_scores = in_fastq.readline()
            barcode = full_seq[0:6]
            ligation = full_seq[3:9]
            trimmed_seq = full_seq[6:]
            out_fastq.write(header.split(" ")[0]+'_MolecularBarcode:'+barcode+' '+header.split(" ")[1]+'\n')
            out_fastq.write(trimmed_seq)
            out_fastq.write(plus_line)
            out_fastq.write(qual_scores[6:])
            header = in_fastq.readline().rstrip()

            barcodes_dict.update([barcode])
            if len(trimmed_seq) >= 4 :
                ligations_dict.update([ligation])

        out_barcode.write("barcode\tcount\n")
        for barcode, count in barcodes_dict.items():
            out_barcode.write(f"{barcode}\t{count}\n")

        out_ligation.write("ligation\tcount\n")
        for ligation, count in ligations_dict.items():
            out_ligation.write(f"{ligation}\t{count}\n")

if __name__ == '__main__':
    main()

