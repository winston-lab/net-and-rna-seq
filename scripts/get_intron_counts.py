#!/usr/bin/env python

import argparse
import collections
import HTSeq as htseq

def main(intron_annotation_path = "S_cerevisiae.R64-2-1_introns_verifiedcoding_nomito_no5utr.gff",
         swap_alignment_strand=True,
         closed_intervals = True,
         alignment_path = "test.bam",
         min_overhang = 4,
         sample_id = "test_sample",
         group_id = "test_group",
         output_path = "results.tsv"):

    intron_annotation = htseq.GFF_Reader(intron_annotation_path,
            end_included=closed_intervals)
    alignments = htseq.BAM_Reader(alignment_path)


    # Make a counter for the number of times an intron name is
    # seen, so we can look up when to append a count to the
    # intron name in order to keep intron names unique
    intron_id_counts = collections.Counter()
    for feature in intron_annotation:
        intron_id_counts.update([feature.attr["Name"]])

    seen_intron_ids = collections.Counter()
    intron_intervals = dict()
    intron_map = htseq.GenomicArrayOfSets("auto", stranded=True)

    for feature in intron_annotation:
        intron_id = feature.attr["Name"]

        # Append a count to non-unique intron names
        if intron_id_counts[intron_id] > 1:
            seen_intron_ids.update([intron_id])
            intron_id = f"{intron_id}_{str(seen_intron_ids[intron_id])}"
            feature.attr["Name"] = intron_id

        intron_intervals[intron_id] = feature.iv

        # Update map of genomic position to overlapping annotations with intron name
        intron_map[feature.iv] += intron_id

    read_counts = {alignment_type: collections.Counter() for \
            alignment_type in ("ambiguous",
                               "spliced",
                               "intronic",
                               "junction_five",
                               "junction_three")}

    # bam_writers = {"ambiguous": htseq.BAM_Writer.from_BAM_Reader("ambiguous.bam", alignments),
    #                "spliced": htseq.BAM_Writer.from_BAM_Reader("spliced.bam", alignments),
    #                "intronic": htseq.BAM_Writer.from_BAM_Reader("intronic.bam", alignments),
    #                "junction_five": htseq.BAM_Writer.from_BAM_Reader("junction_five.bam", alignments),
    #                "junction_three": htseq.BAM_Writer.from_BAM_Reader("junction_three.bam", alignments)}

    for alignment in alignments:
        # when sequencing from 3' end, the strand is swapped
        if swap_alignment_strand:
            alignment.iv.strand = {"+": "-",
                                   "-": "+"}.get(alignment.iv.strand)

        # for each alignment, find overlapping introns
        overlapped_introns = set()
        for interval, value in intron_map[alignment.iv].steps():
            overlapped_introns |= value

        # ignore alignments not spanning any introns
        if len(overlapped_introns) == 0:
            continue

        # mark alignments spanning multiple introns as ambiguous
        if len(overlapped_introns) > 1:
            read_counts["ambiguous"].update(overlapped_introns)
            # bam_writers["ambiguous"].write(alignment)
            continue

        cigar = alignment.cigar
        cigar_length = len(cigar)

        # mark alignments with complex CIGAR strings as ambiguous
        if cigar_length not in [1, 3]:
            read_counts["ambiguous"].update(overlapped_introns)
            # bam_writers["ambiguous"].write(alignment)
            continue

        overlapped_intron = list(overlapped_introns)[0]

        # handle potentially spliced alignments
        if cigar_length == 3:
            if [x.type for x in cigar] != ["M", "N", "M"]:
                read_counts["ambiguous"].update(overlapped_introns)
                # bam_writers["ambiguous"].write(alignment)
                continue
            if cigar[1].ref_iv.start != intron_intervals[overlapped_intron].start or \
                    cigar[1].ref_iv.end != intron_intervals[overlapped_intron].end or \
                    cigar[0].ref_iv.end != intron_intervals[overlapped_intron].start or \
                    cigar[2].ref_iv.start != intron_intervals[overlapped_intron].end or \
                    cigar[0].size < min_overhang or \
                    cigar[2].size < min_overhang:
                read_counts["ambiguous"].update(overlapped_introns)
                # bam_writers["ambiguous"].write(alignment)
                continue
            read_counts["spliced"].update(overlapped_introns)
            # bam_writers["spliced"].write(alignment)
            continue

        # handle potential junction or intronic reads
        if cigar[0].type != "M":
            read_counts["ambiguous"].update(overlapped_introns)
            # bam_writers["ambiguous"].write(alignment)
            continue

        if cigar[0].ref_iv.start >= intron_intervals[overlapped_intron].start and \
                cigar[0].ref_iv.end <= intron_intervals[overlapped_intron].end:
            read_counts["intronic"].update(overlapped_introns)
            # bam_writers["intronic"].write(alignment)
            continue

        if cigar[0].ref_iv.start <= (intron_intervals[overlapped_intron].start - min_overhang) and \
                cigar[0].ref_iv.end >= (intron_intervals[overlapped_intron].start + min_overhang):
            ({"+": read_counts["junction_five"],
              "-": read_counts["junction_three"]}.get(alignment.iv.strand)).update(overlapped_introns)
            # {"+": bam_writers["junction_five"],
            #  "-": bam_writers["junction_three"]}.get(alignment.iv.strand).write(alignment)
            continue
        if cigar[0].ref_iv.start <= (intron_intervals[overlapped_intron].end - min_overhang) and \
                cigar[0].ref_iv.end >= (intron_intervals[overlapped_intron].end + min_overhang):
            ({"+": read_counts["junction_three"],
              "-": read_counts["junction_five"]}.get(alignment.iv.strand)).update(overlapped_introns)
            # {"+": bam_writers["junction_three"],
            #  "-": bam_writers["junction_five"]}.get(alignment.iv.strand).write(alignment)
            continue
        read_counts["ambiguous"].update(overlapped_introns)

    # for writer in bam_writers.values():
    #     writer.close()

    with open(output_path, "w") as output_file:
        output_file.write("\t".join(["chrom",
                                     "start",
                                     "end",
                                     "name",
                                     "score",
                                     "strand",
                                     "sample_id",
                                     "group_id",
                                     "spliced",
                                     "junction_5",
                                     "junction_3",
                                     "intronic",
                                     "ambiguous"]) + "\n")
        for key, interval in intron_intervals.items():
            output_string = "\t".join([interval.chrom,
                                       str(interval.start),
                                       str(interval.end),
                                       key,
                                       "0",
                                       interval.strand,
                                       sample_id,
                                       group_id,
                                       str(read_counts["spliced"][key]),
                                       str(read_counts["junction_five"][key]),
                                       str(read_counts["junction_three"][key]),
                                       str(read_counts["intronic"][key]),
                                       str(read_counts["ambiguous"][key])]) + "\n"
            output_file.write(output_string)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Given a GFF of intron annotatations and BAM of (filtered) RNA-seq alignments, count spliced, junction, and intronic alignments for each intron.")
    parser.add_argument(
            '-i', '--introns',
            dest='intron_annotation_path',
            type=str,
            help='GFF of intron annotations.')
    parser.add_argument(
            '-a', '--alignment',
            dest='alignment_path',
            type=str,
            help='BAM of (filtered) RNA-seq alignments.')
    parser.add_argument('-sw', '--swap_alignment_strand',
            dest='swap_alignment_strand',
            action="store_true",
            help='Flag to indicate that alignment strand is opposite of RNA strand, which occurs when sequencing from the 3-prime end of the RNA.')
    parser.add_argument('-op', '--open_intervals',
            dest='closed_intervals',
            action="store_false",
            help='Flag to indicate that regions in intron GFF are 0-indexed, half-open (as in BED format) versus 1-indexed, fully-closed (default, as in Ensembl GTFs.)')
    parser.add_argument(
            '-s', '--sample_id',
            dest='sample_id',
            type=str,
            help='Experimental sample name for metadata.')
    parser.add_argument(
            '-g', '--group_id',
            dest='group_id',
            type=str,
            help='Experimental group name for metadata.')
    parser.add_argument('-m', '--min_overhang',
            dest='min_overhang',
            type=int,
            default=4,
            help="Minimum number of bases an alignment must span on either side of a junction to be assigned unambiguously.")
    parser.add_argument(
            '-o', '--output',
            dest='output_path',
            type=str,
            help='Path to write tsv of results.')
    args = parser.parse_args()
    main(intron_annotation_path=args.intron_annotation_path,
         alignment_path=args.alignment_path,
         swap_alignment_strand=args.swap_alignment_strand,
         closed_intervals=args.closed_intervals,
         sample_id=args.sample_id,
         group_id=args.group_id,
         min_overhang=args.min_overhang,
         output_path=args.output_path)

