#!/usr/bin/env python3
import sys

def parse_fasta_lengths(fasta_path):
    """
    Parse a FASTA file and return a dictionary of {chrom: length}.
    Assumes the FASTA may contain multiple sequences (chromosomes).
    """
    lengths = {}
    current_chrom = None
    current_len = 0

    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # New chromosome header
                if current_chrom is not None:
                    lengths[current_chrom] = current_len
                current_chrom = line[1:].strip().split()[0]
                current_len = 0
            else:
                # Sequence line
                seq = line.strip()
                current_len += len(seq)
        # Last chromosome
        if current_chrom is not None:
            lengths[current_chrom] = current_len

    return lengths

def parse_gtf(gtf_path):
    """
    Parse the GTF file and return two dictionaries:
    coding_intervals[chrom] = list of (start, end) for exons
    noncoding_intervals[chrom] = list of (start, end) for all other features
    """
    coding_intervals = {}
    noncoding_intervals = {}

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attr = fields
            start = int(start)
            end = int(end)
            if feature.lower() == "exon":
                coding_intervals.setdefault(chrom, []).append((start, end))
            else:
                noncoding_intervals.setdefault(chrom, []).append((start, end))
    return coding_intervals, noncoding_intervals

def merge_intervals(intervals):
    """Merge overlapping intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for curr in intervals[1:]:
        prev = merged[-1]
        if curr[0] <= prev[1] + 1:
            merged[-1] = (prev[0], max(prev[1], curr[1]))
        else:
            merged.append(curr)
    return merged

def intervals_to_nucleotides(intervals_dict):
    """
    Convert intervals into a dict of sets, each containing all covered positions.
    """
    pos_dict = {}
    for chrom, intervals in intervals_dict.items():
        merged = merge_intervals(intervals)
        pos_set = set()
        for (start, end) in merged:
            pos_set.update(range(start, end+1))
        pos_dict[chrom] = pos_set
    return pos_dict

def genome_as_nucleotides(genome_lengths):
    """
    Create a dictionary of sets representing every nucleotide in the genome.
    {chrom: set_of_positions_for_that_chrom}
    """
    universe = {}
    for chrom, length in genome_lengths.items():
        # This can be huge for large genomes.
        universe[chrom] = set(range(1, length+1))
    return universe

def dict_size(pos_dict):
    """Return the total number of nucleotides across all chromosomes in pos_dict."""
    total = 0
    for s in pos_dict.values():
        total += len(s)
    return total

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: {} genome.fa reference.gtf prediction.gtf".format(sys.argv[0]))
        sys.exit(1)

    genome_fa = sys.argv[1]
    ref_gtf = sys.argv[2]
    pred_gtf = sys.argv[3]

    # Parse the genome to get chromosome lengths
    genome_lengths = parse_fasta_lengths(genome_fa)

    # Turn the entire genome into a set of positions
    genome_universe = genome_as_nucleotides(genome_lengths)

    # Parse GTFs
    ref_coding, ref_noncoding = parse_gtf(ref_gtf)
    pred_coding, pred_noncoding = parse_gtf(pred_gtf)

    # Convert intervals to sets of nucleotides for coding
    ref_coding_set = intervals_to_nucleotides(ref_coding)
    pred_coding_set = intervals_to_nucleotides(pred_coding)

    # Now assign all nucleotides not in coding sets as non-coding, using the full genome as baseline
    ref_coding_full = {}
    ref_noncoding_full = {}
    pred_coding_full = {}
    pred_noncoding_full = {}

    for chrom, uni_positions in genome_universe.items():
        r_c = ref_coding_set.get(chrom, set())
        p_c = pred_coding_set.get(chrom, set())

        # Non-coding is universe minus coding
        r_nc = uni_positions - r_c
        p_nc = uni_positions - p_c

        ref_coding_full[chrom] = r_c
        ref_noncoding_full[chrom] = r_nc
        pred_coding_full[chrom] = p_c
        pred_noncoding_full[chrom] = p_nc

    # Calculate totals
    ref_coding_total = dict_size(ref_coding_full)
    ref_noncoding_total = dict_size(ref_noncoding_full)
    ref_total = ref_coding_total + ref_noncoding_total

    pred_coding_total = dict_size(pred_coding_full)
    pred_noncoding_total = dict_size(pred_noncoding_full)
    pred_total = pred_coding_total + pred_noncoding_total

    # Print totals
    print("Reference GTF (with gaps as non-coding, entire genome considered):")
    print(f"  Total coding nucleotides:    {ref_coding_total}")
    print(f"  Total non-coding nucleotides:{ref_noncoding_total}")
    print(f"  Sum (coding + non-coding):   {ref_total}")

    print("\nPrediction GTF (with gaps as non-coding, entire genome considered):")
    print(f"  Total coding nucleotides:    {pred_coding_total}")
    print(f"  Total non-coding nucleotides:{pred_noncoding_total}")
    print(f"  Sum (coding + non-coding):   {pred_total}")

    # Since now we consider the entire genome, ref_total and pred_total should match
    # (assuming the same chromosomes referenced in both).
    # If there are chromosomes in the genome not referenced at all, they are still included,
    # ensuring both cover the full genome.

    # Compute TP, FP, FN, TN
    # Positive = coding, Negative = non-coding
    TP = 0
    FP = 0
    FN = 0
    TN = 0

    all_chroms = genome_universe.keys()
    for chrom in all_chroms:
        ref_coding_pos = ref_coding_full[chrom]
        ref_non_coding_pos = ref_noncoding_full[chrom]
        pred_coding_pos = pred_coding_full[chrom]
        pred_non_coding_pos = pred_noncoding_full[chrom]

        # TP: coding in both ref and pred
        TP += len(ref_coding_pos & pred_coding_pos)
        # FP: predicted coding but reference non-coding
        FP += len(pred_coding_pos & ref_non_coding_pos)
        # FN: reference coding but predicted non-coding
        FN += len(ref_coding_pos & pred_non_coding_pos)
        # TN: non-coding in both
        TN += len(ref_non_coding_pos & pred_non_coding_pos)

    # Compute metrics
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    total_eval = TP + FP + FN + TN
    TP_pct = (TP / total_eval * 100) if total_eval > 0 else 0.0
    FP_pct = (FP / total_eval * 100) if total_eval > 0 else 0.0
    FN_pct = (FN / total_eval * 100) if total_eval > 0 else 0.0
    TN_pct = (TN / total_eval * 100) if total_eval > 0 else 0.0

    precision_pct = precision * 100
    recall_pct = recall * 100
    f1_pct = f1 * 100

    print("\nEvaluation (Coding = Positive, Non-coding = Negative):")
    print(f"  TP: {TP} ({TP_pct:.2f}%)")
    print(f"  FP: {FP} ({FP_pct:.2f}%)")
    print(f"  FN: {FN} ({FN_pct:.2f}%)")
    print(f"  TN: {TN} ({TN_pct:.2f}%)")

    print("\nMetrics:")
    print(f"  Precision: {precision:.4f} ({precision_pct:.2f}%)")
    print(f"  Recall:    {recall:.4f} ({recall_pct:.2f}%)")
    print(f"  F1-score:  {f1:.4f} ({f1_pct:.2f}%)")
