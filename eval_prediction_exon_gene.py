#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import os

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

def main():
    parser = argparse.ArgumentParser(
        description="Calculate feature counts and evaluation metrics between reference and predicted GTF annotations."
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
        '-o', '--output',
        required=False,
        default="annotation_evaluation_report.txt",
        help="Path to the output report file."
    )
    args = parser.parse_args()

    # Check if files exist
    if not os.path.isfile(args.reference):
        print(f"Reference GTF file not found: {args.reference}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.predicted):
        print(f"Predicted GTF file not found: {args.predicted}", file=sys.stderr)
        sys.exit(1)

    # Parse GTF files
    print("Parsing reference GTF file...")
    ref_gtf = parse_gtf(args.reference)
    print("Parsing predicted GTF file...")
    pred_gtf = parse_gtf(args.predicted)

    # Define feature types to evaluate
    feature_types = ['gene', 'exon']

    # Extract features
    print("Extracting reference features...")
    ref_features = extract_features(ref_gtf, feature_types)
    print("Extracting predicted features...")
    pred_features = extract_features(pred_gtf, feature_types)

    # Calculate total features
    total_ref_features = len(ref_gtf)
    total_pred_features = len(pred_gtf)

    # Count gene and non-gene features
    ref_gene_count = len(ref_features.get('gene', set()))
    pred_gene_count = len(pred_features.get('gene', set()))
    ref_non_gene_count = total_ref_features - ref_gene_count
    pred_non_gene_count = total_pred_features - pred_gene_count

    # Count exon and non-exon features
    ref_exon_count = len(ref_features.get('exon', set()))
    pred_exon_count = len(pred_features.get('exon', set()))
    ref_non_exon_count = total_ref_features - ref_exon_count
    pred_non_exon_count = total_pred_features - ref_exon_count  # Note: Using ref_exon_count for non-exon
    # Correction: It should be total_ref_features - exon_count
    ref_non_exon_count = total_ref_features - ref_exon_count
    pred_non_exon_count = total_pred_features - pred_exon_count

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

    # Calculate overall metrics (all features)
    print("Evaluating All Features...")
    # For overall, include feature type in matching to ensure exact matches
    # Create sets with feature type included
    ref_all_features = set(
        zip(
            ref_gtf['seqname'],
            ref_gtf['feature'],
            ref_gtf['start'],
            ref_gtf['end'],
            ref_gtf['strand']
        )
    )
    pred_all_features = set(
        zip(
            pred_gtf['seqname'],
            pred_gtf['feature'],
            pred_gtf['start'],
            pred_gtf['end'],
            pred_gtf['strand']
        )
    )
    TP_all, FP_all, FN_all, precision_all, recall_all, f1_all = calculate_metrics(ref_all_features, pred_all_features)

    # Initialize report
    report_lines = []
    report_lines.append("### Annotation Evaluation Report ###\n")
    report_lines.append(f"Reference GTF File: {args.reference}")
    report_lines.append(f"Predicted GTF File: {args.predicted}\n")

    # Total Features
    report_lines.append("#### Total Features ####")
    report_lines.append(f"Total Reference Features: {total_ref_features}")
    report_lines.append(f"Total Predicted Features: {total_pred_features}\n")

    # Gene and Non-Gene Counts
    report_lines.append("#### Gene and Non-Gene Features ####")
    report_lines.append(f"Reference Genes: {ref_gene_count}")
    report_lines.append(f"Predicted Genes: {pred_gene_count}")
    report_lines.append(f"Reference Non-Gene Features: {ref_non_gene_count}")
    report_lines.append(f"Predicted Non-Gene Features: {pred_non_gene_count}\n")

    # Exon and Non-Exon Counts
    report_lines.append("#### Exon and Non-Exon Features ####")
    report_lines.append(f"Reference Exons: {ref_exon_count}")
    report_lines.append(f"Predicted Exons: {pred_exon_count}")
    report_lines.append(f"Reference Non-Exon Features: {ref_non_exon_count}")
    report_lines.append(f"Predicted Non-Exon Features: {pred_non_exon_count}\n")

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