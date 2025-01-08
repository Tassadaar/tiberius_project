#!/usr/bin/env python3

"""
GFF3 to GTF Converter

This script converts a GFF3 file to GTF format, ensuring compatibility with tools that require GTF annotations.

Usage:
    python gff3_to_gtf.py -i input.gff3 -o output.gtf

Requirements:
    - Python 3.x
    - BioPython

Author: [Your Name]
Date: [Date]
"""

import argparse
import sys
import re

def parse_arguments():
    """
    Parses command-line arguments.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Convert GFF3 file to GTF format."
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help="Path to the input GFF3 file."
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help="Path to the output GTF file."
    )
    return parser.parse_args()

def parse_attributes(attribute_string):
    """
    Parses the attribute column from GFF3 and converts it to GTF format.

    GFF3 attributes are semicolon-separated key-value pairs.
    GTF requires attributes to be in the format key "value";

    Args:
        attribute_string (str): The attribute string from GFF3.

    Returns:
        str: Formatted attribute string for GTF.
    """
    attributes = []
    # Split attributes by semicolon, considering possible escaped semicolons
    attr_pairs = re.findall(r'(\S+?="[^"]*")', attribute_string)
    for pair in attr_pairs:
        key, value = pair.split('=', 1)
        # Remove quotes if present and handle spaces
        value = value.strip('"')
        # GTF requires a trailing semicolon and quotes around values
        attributes.append(f'{key} "{value}";')
    # Join all attributes with space
    return ' '.join(attributes)

def convert_gff3_to_gtf(input_file, output_file):
    """
    Converts a GFF3 file to GTF format.

    Args:
        input_file (str): Path to the input GFF3 file.
        output_file (str): Path to the output GTF file.
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) != 9:
                    print(f"Skipping malformed line (does not have 9 columns): {line}", file=sys.stderr)
                    continue
                seqname, source, feature, start, end, score, strand, frame, attribute = parts

                # Parse attributes to extract gene_id and transcript_id
                attr_dict = {}
                # GFF3 attributes are separated by semicolons
                attrs = attribute.split(';')
                for attr in attrs:
                    if '=' not in attr:
                        continue
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

                # Ensure required attributes are present
                gene_id = attr_dict.get('gene_id') or attr_dict.get('geneID') or attr_dict.get('geneName') or 'gene_unknown'
                transcript_id = attr_dict.get('transcript_id') or attr_dict.get('transcriptID') or attr_dict.get('transcriptName') or 'transcript_unknown'

                # Create GTF attribute string
                gtf_attributes = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'

                # GTF requires specific fields; mapping GFF3 fields to GTF
                # GTF fields: seqname, source, feature, start, end, score, strand, frame, attribute
                # Assuming frame is the same; GTF usually expects frame to be '0', '1', '2', or '.' for non-CDS features
                # Ensure frame is set correctly
                if feature.lower() not in ['cds', 'start_codon', 'stop_codon']:
                    gtf_frame = '.'
                else:
                    gtf_frame = frame if frame in ['0', '1', '2'] else '.'

                gtf_line = '\t'.join([
                    seqname,
                    source,
                    feature,
                    start,
                    end,
                    score,
                    strand,
                    gtf_frame,
                    gtf_attributes
                ]) + '\n'

                outfile.write(gtf_line)
        print(f"Conversion completed successfully. GTF file saved to {output_file}")
    except FileNotFoundError:
        print(f"Error: File {input_file} not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during conversion: {e}", file=sys.stderr)

def main():
    args = parse_arguments()
    convert_gff3_to_gtf(args.input, args.output)

if __name__ == "__main__":
    main()