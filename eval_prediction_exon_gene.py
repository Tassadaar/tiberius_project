#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os
import re

def parse_gtf(gtf_file):
    """
    Parses a GTF file and returns a DataFrame with relevant columns.

    Parameters:
        gtf_file (str): Path to the GTF file.

    Returns:
        pd.DataFrame: DataFrame containing parsed GTF data.
    """
    try:
        # Define column names as per GTF specification
        column_names = [
            'seqname', 'source', 'feature', 'start', 'end',
            'score', 'strand', 'frame', 'attribute'
        ]

        # Read the GTF file using pandas
        gtf_df = pd.read_csv(
            gtf_file,
            sep='\t',
            comment='#',
            header=None,
            names=column_names,
            dtype={
                'seqname': str,
                'source': str,
                'feature': str,
                'start': int,
                'end': int,
                'score': str,
                'strand': str,
                'frame': str,
                'attribute': str
            }
        )

        return gtf_df

    except Exception as e:
        print(f"Error parsing GTF file {gtf_file}: {e}", file=sys.stderr)
        sys.exit(1)

def filter_by_strand(gtf_df, allowed_strands=None):
    """
    Filters the GTF DataFrame to include only specified strands.

    Parameters:
        gtf_df (pd.DataFrame): DataFrame containing GTF data.
        allowed_strands (list, optional): List of strands to include (e.g., ['+', '.']).
                                          Defaults to ['+', '.'].

    Returns:
        pd.DataFrame: Filtered DataFrame containing only the specified strands.
    """
    if allowed_strands is None:
        allowed_strands = ['+']
    # Ensure strands are strings and uppercase
    gtf_df['strand'] = gtf_df['strand'].astype(str).str.upper()
    filtered_df = gtf_df[gtf_df['strand'].isin(allowed_strands)].copy()
    return filtered_df

def parse_genome_fasta(genome_fasta):
    """
    Parses a genome FASTA file and returns a dictionary mapping chromosome names to their lengths.

    Parameters:
        genome_fasta (str): Path to the genome FASTA file.

    Returns:
        dict: Dictionary mapping chromosome names to their sizes.
    """
    try:
        chrom_sizes = {}
        chrom = None
        size = 0
        with open(genome_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if chrom:
                        chrom_sizes[chrom] = size
                    chrom = line[1:].split()[0]  # Take the first word after '>'
                    size = 0
                else:
                    size += len(line)
            if chrom:
                chrom_sizes[chrom] = size
        return chrom_sizes
    except Exception as e:
        print(f"Error parsing genome FASTA file {genome_fasta}: {e}", file=sys.stderr)
        sys.exit(1)

def extract_features(gtf_df, feature_types):
    """
    Extracts features of specified types from the GTF DataFrame.

    Parameters:
        gtf_df (pd.DataFrame): DataFrame containing GTF data.
        feature_types (list): List of feature types to extract (e.g., ['gene', 'exon']).

    Returns:
        dict: Dictionary with feature types as keys and sets of feature tuples as values.
    """
    features = {}
    for feature in feature_types:
        subset = gtf_df[gtf_df['feature'] == feature]
        # Create tuples of (seqname, start, end, strand)
        feature_set = set(
            zip(
                subset['seqname'],
                subset['start'],
                subset['end'],
                subset['strand']
            )
        )
        features[feature] = feature_set
    return features

def calculate_metrics(ref_features, pred_features):
    """
    Calculates precision, recall, and F1-score based on reference and predicted features.

    Parameters:
        ref_features (set): Set of reference feature tuples.
        pred_features (set): Set of predicted feature tuples.

    Returns:
        tuple: (TP, FP, FN, precision, recall, f1_score)
    """
    TP = len(ref_features & pred_features)
    FP = len(pred_features - ref_features)
    FN = len(ref_features - pred_features)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0.0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0.0
    f1_score = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    return TP, FP, FN, precision, recall, f1_score

def infer_introns(ref_gtf_df):
    """
    Infers intron regions from exons within each gene.

    Parameters:
        ref_gtf_df (pd.DataFrame): Filtered reference GTF DataFrame.

    Returns:
        set: Set of intron feature tuples (seqname, start, end, strand).
    """
    intron_set = set()
    # Group by gene
    genes = ref_gtf_df[ref_gtf_df['feature'] == 'gene']['attribute']
    # Extract gene IDs
    gene_ids = ref_gtf_df[ref_gtf_df['feature'] == 'gene']['attribute'].apply(
        lambda x: re.search('gene_id "([^"]+)"', x).group(1) if re.search('gene_id "([^"]+)"', x) else None
    )
    ref_gtf_df = ref_gtf_df.copy()
    ref_gtf_df['gene_id'] = gene_ids
    # Iterate over each gene
    for gene_id, group in ref_gtf_df.groupby('gene_id'):
        exons = group[group['feature'] == 'exon'].sort_values('start')
        if len(exons) < 2:
            continue  # No intron if less than 2 exons
        for i in range(len(exons) - 1):
            intron_start = exons.iloc[i]['end'] + 1
            intron_end = exons.iloc[i + 1]['start'] - 1
            if intron_start <= intron_end:
                intron_tuple = (
                    exons.iloc[i]['seqname'],
                    intron_start,
                    intron_end,
                    exons.iloc[i]['strand']
                )
                intron_set.add(intron_tuple)
    return intron_set

def compute_intergenic_regions(chrom_sizes, gene_features):
    """
    Computes intergenic regions based on chromosome sizes and gene features.

    Parameters:
        chrom_sizes (dict): Dictionary mapping chromosome names to their sizes.
        gene_features (set): Set of gene feature tuples (seqname, start, end, strand).

    Returns:
        set: Set of intergenic feature tuples (seqname, start, end, strand).
    """
    intergenic_set = set()
    # Organize genes by chromosome
    genes_by_chrom = {}
    for gene in gene_features:
        chrom = gene[0]
        start = gene[1]
        end = gene[2]
        strand = gene[3]
        if chrom not in genes_by_chrom:
            genes_by_chrom[chrom] = []
        genes_by_chrom[chrom].append((start, end))
    # For each chromosome, compute intergenic regions
    for chrom, size in chrom_sizes.items():
        if chrom not in genes_by_chrom:
            # Entire chromosome is intergenic
            intergenic_set.add((chrom, 1, size, '.'))
            continue
        sorted_genes = sorted(genes_by_chrom[chrom], key=lambda x: x[0])
        prev_end = 0
        for gene_start, gene_end in sorted_genes:
            if prev_end + 1 < gene_start:
                intergenic_set.add((chrom, prev_end + 1, gene_start - 1, '.'))
            prev_end = max(prev_end, gene_end)
        if prev_end < size:
            intergenic_set.add((chrom, prev_end + 1, size, '.'))
    return intergenic_set

def main():
    parser = argparse.ArgumentParser(
        description="Calculate feature counts and evaluation metrics between reference and predicted GTF annotations, including intronic and intergenic regions."
    )
    parser.add_argument(
        '-r', '--reference',
        required=True,
        help="Path to the reference GTF file."
    )
    parser.add_argument(
        '-p', '--predicted',
        required=True,
        help="Path to the predicted GTF file."
    )
    parser.add_argument(
        '-g', '--genome',
        required=True,
        help="Path to the genome FASTA file."
    )
    parser.add_argument(
        '-o', '--output',
        required=False,
        default="annotation_evaluation_report.txt",
        help="Path to the output report file."
    )
    args = parser.parse_args()

    # Check if files exist
    for file_path, desc in [(args.reference, "Reference GTF"),
                            (args.predicted, "Predicted GTF"),
                            (args.genome, "Genome FASTA")]:
        if not os.path.isfile(file_path):
            print(f"{desc} file not found: {file_path}", file=sys.stderr)
            sys.exit(1)

    # Parse genome FASTA to get chromosome sizes
    print("Parsing genome FASTA file to extract chromosome sizes...")
    chrom_sizes = parse_genome_fasta(args.genome)
    print(f"Extracted chromosome sizes for {len(chrom_sizes)} chromosomes.")

    # Parse GTF files
    print("\nParsing reference GTF file...")
    ref_gtf = parse_gtf(args.reference)
    print("Parsing predicted GTF file...")
    pred_gtf = parse_gtf(args.predicted)

    # Filter by strand: only '+' and '.' strands
    allowed_strands = ['+', '.']
    print(f"\nFiltering reference GTF for strands: {allowed_strands}")
    ref_gtf_filtered = filter_by_strand(ref_gtf, allowed_strands=allowed_strands)
    print(f"Reference GTF after filtering: {len(ref_gtf_filtered)} features")

    print(f"\nFiltering predicted GTF for strands: {allowed_strands}")
    pred_gtf_filtered = filter_by_strand(pred_gtf, allowed_strands=allowed_strands)
    print(f"Predicted GTF after filtering: {len(pred_gtf_filtered)} features")

    # Define feature types to evaluate
    feature_types = ['gene', 'exon']

    # Extract features
    print("\nExtracting reference features...")
    ref_features = extract_features(ref_gtf_filtered, feature_types)
    print("Extracting predicted features...")
    pred_features = extract_features(pred_gtf_filtered, feature_types)

    # Calculate total features
    total_ref_features = len(ref_gtf_filtered)
    total_pred_features = len(pred_gtf_filtered)

    # Count gene and non-gene features
    ref_gene_count = len(ref_features.get('gene', set()))
    pred_gene_count = len(pred_features.get('gene', set()))
    ref_non_gene_count = total_ref_features - ref_gene_count
    pred_non_gene_count = total_pred_features - pred_gene_count

    # Count exon and non-exon features
    ref_exon_count = len(ref_features.get('exon', set()))
    pred_exon_count = len(pred_features.get('exon', set()))
    ref_non_exon_count = total_ref_features - ref_exon_count
    pred_non_exon_count = total_pred_features - pred_exon_count

    # Infer introns from reference GTF
    print("\nInferring intron regions from reference GTF...")
    ref_introns = infer_introns(ref_gtf_filtered)
    ref_intron_count = len(ref_introns)
    print(f"Number of inferred introns in reference: {ref_intron_count}")

    # Infer introns from predicted GTF
    print("Inferring intron regions from predicted GTF...")
    pred_introns = infer_introns(pred_gtf_filtered)
    pred_intron_count = len(pred_introns)
    print(f"Number of inferred introns in predicted: {pred_intron_count}")

    # Compute intergenic regions from reference GTF
    print("\nComputing intergenic regions from reference GTF...")
    ref_intergenic = compute_intergenic_regions(chrom_sizes, ref_features.get('gene', set()))
    ref_intergenic_count = len(ref_intergenic)
    print(f"Number of intergenic regions in reference: {ref_intergenic_count}")

    # Compute intergenic regions from predicted GTF
    print("Computing intergenic regions from predicted GTF...")
    pred_intergenic = compute_intergenic_regions(chrom_sizes, pred_features.get('gene', set()))
    pred_intergenic_count = len(pred_intergenic)
    print(f"Number of intergenic regions in predicted: {pred_intergenic_count}")

    # Adjust non-gene and non-exon counts
    # Non-gene features now include introns and intergenic regions
    ref_non_gene_adjusted = ref_non_gene_count - ref_exon_count  # Remove exons from non-gene
    pred_non_gene_adjusted = pred_non_gene_count - pred_exon_count

    # Calculate metrics for genes
    print("\nEvaluating Gene Features...")
    ref_genes = ref_features.get('gene', set())
    pred_genes = pred_features.get('gene', set())
    TP_gene, FP_gene, FN_gene, precision_gene, recall_gene, f1_gene = calculate_metrics(ref_genes, pred_genes)

    # Calculate metrics for exons
    print("Evaluating Exon Features...")
    ref_exons = ref_features.get('exon', set())
    pred_exons = pred_features.get('exon', set())
    TP_exon, FP_exon, FN_exon, precision_exon, recall_exon, f1_exon = calculate_metrics(ref_exons, pred_exons)

    # Calculate metrics for introns
    print("Evaluating Intron Features...")
    TP_intron, FP_intron, FN_intron, precision_intron, recall_intron, f1_intron = calculate_metrics(ref_introns, pred_introns)

    # Calculate metrics for intergenic regions
    print("Evaluating Intergenic Features...")
    TP_intergenic, FP_intergenic, FN_intergenic, precision_intergenic, recall_intergenic, f1_intergenic = calculate_metrics(ref_intergenic, pred_intergenic)

    # Calculate overall metrics (all features)
    print("\nEvaluating All Features...")
    # For overall, include feature type in matching to ensure exact matches
    # Create sets with feature type included
    ref_all_features = set(
        zip(
            ref_gtf_filtered['seqname'],
            ref_gtf_filtered['feature'],
            ref_gtf_filtered['start'],
            ref_gtf_filtered['end'],
            ref_gtf_filtered['strand']
        )
    )
    pred_all_features = set(
        zip(
            pred_gtf_filtered['seqname'],
            pred_gtf_filtered['feature'],
            pred_gtf_filtered['start'],
            pred_gtf_filtered['end'],
            pred_gtf_filtered['strand']
        )
    )
    TP_all, FP_all, FN_all, precision_all, recall_all, f1_all = calculate_metrics(ref_all_features, pred_all_features)

    # Initialize report
    report_lines = []
    report_lines.append("### Annotation Evaluation Report ###\n")
    report_lines.append(f"Reference GTF File: {args.reference}")
    report_lines.append(f"Predicted GTF File: {args.predicted}")
    report_lines.append(f"Genome FASTA File: {args.genome}\n")

    # Total Features
    report_lines.append("#### Total Features ####")
    report_lines.append(f"Total Reference Features (Strand '+' or '.'): {total_ref_features}")
    report_lines.append(f"Total Predicted Features (Strand '+' or '.'): {total_pred_features}\n")

    # Gene and Non-Gene Counts
    report_lines.append("#### Gene and Non-Gene Features ####")
    report_lines.append(f"Reference Genes: {ref_gene_count}")
    report_lines.append(f"Predicted Genes: {pred_gene_count}")
    report_lines.append(f"Reference Non-Gene Features (Exons excluded): {ref_non_gene_adjusted}")
    report_lines.append(f"Predicted Non-Gene Features (Exons excluded): {pred_non_gene_adjusted}\n")

    # Exon and Non-Exon Counts
    report_lines.append("#### Exon and Non-Exon Features ####")
    report_lines.append(f"Reference Exons: {ref_exon_count}")
    report_lines.append(f"Predicted Exons: {pred_exon_count}")
    report_lines.append(f"Reference Non-Exon Features: {ref_non_exon_count}")
    report_lines.append(f"Predicted Non-Exon Features: {pred_non_exon_count}\n")

    # Intron and Intergenic Counts
    report_lines.append("#### Intron and Intergenic Features ####")
    report_lines.append(f"Reference Introns: {ref_intron_count}")
    report_lines.append(f"Predicted Introns: {pred_intron_count}")
    report_lines.append(f"Reference Intergenic Regions: {ref_intergenic_count}")
    report_lines.append(f"Predicted Intergenic Regions: {pred_intergenic_count}\n")

    # Metrics for Genes
    report_lines.append("#### Gene Feature Metrics ####")
    report_lines.append(f"True Positives (TP): {TP_gene}")
    report_lines.append(f"False Positives (FP): {FP_gene}")
    report_lines.append(f"False Negatives (FN): {FN_gene}")
    report_lines.append(f"Precision: {precision_gene:.4f}")
    report_lines.append(f"Recall: {recall_gene:.4f}")
    report_lines.append(f"F1-Score: {f1_gene:.4f}\n")

    # Metrics for Exons
    report_lines.append("#### Exon Feature Metrics ####")
    report_lines.append(f"True Positives (TP): {TP_exon}")
    report_lines.append(f"False Positives (FP): {FP_exon}")
    report_lines.append(f"False Negatives (FN): {FN_exon}")
    report_lines.append(f"Precision: {precision_exon:.4f}")
    report_lines.append(f"Recall: {recall_exon:.4f}")
    report_lines.append(f"F1-Score: {f1_exon:.4f}\n")

    # Metrics for Introns
    report_lines.append("#### Intron Feature Metrics ####")
    report_lines.append(f"True Positives (TP): {TP_intron}")
    report_lines.append(f"False Positives (FP): {FP_intron}")
    report_lines.append(f"False Negatives (FN): {FN_intron}")
    report_lines.append(f"Precision: {precision_intron:.4f}")
    report_lines.append(f"Recall: {recall_intron:.4f}")
    report_lines.append(f"F1-Score: {f1_intron:.4f}\n")

    # Metrics for Intergenic Regions
    report_lines.append("#### Intergenic Feature Metrics ####")
    report_lines.append(f"True Positives (TP): {TP_intergenic}")
    report_lines.append(f"False Positives (FP): {FP_intergenic}")
    report_lines.append(f"False Negatives (FN): {FN_intergenic}")
    report_lines.append(f"Precision: {precision_intergenic:.4f}")
    report_lines.append(f"Recall: {recall_intergenic:.4f}")
    report_lines.append(f"F1-Score: {f1_intergenic:.4f}\n")

    # Metrics for All Features
    report_lines.append("#### All Features Metrics ####")
    report_lines.append(f"True Positives (TP): {TP_all}")
    report_lines.append(f"False Positives (FP): {FP_all}")
    report_lines.append(f"False Negatives (FN): {FN_all}")
    report_lines.append(f"Precision: {precision_all:.4f}")
    report_lines.append(f"Recall: {recall_all:.4f}")
    report_lines.append(f"F1-Score: {f1_all:.4f}\n")

    # Write report to file
    try:
        with open(args.output, 'w') as f:
            for line in report_lines:
                f.write(line + '\n')
        print(f"\nEvaluation report saved to {args.output}")
    except Exception as e:
        print(f"Error writing report to {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    # Optionally, print report to console
    print("\n=== Annotation Evaluation Report ===\n")
    with open(args.output, 'r') as f:
        print(f.read())

if __name__ == "__main__":
    main()