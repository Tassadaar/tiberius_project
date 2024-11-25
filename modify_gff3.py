#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_attributes(attributes_str):
    """Parse the GFF3 attributes column and return a dictionary"""
    attributes = {}
    for attribute in attributes_str.strip().split(';'):
        if attribute:
            key_value = attribute.strip().split('=', 1)
            if len(key_value) == 2:
                key, value = key_value
                attributes[key.strip()] = value.strip()
    return attributes

def process_gff3(input_file, output_file):
    features_by_seqid = defaultdict(list)  # {seqid: [feature_dicts]}
    seqid_lengths = {}  # {seqid: length}
    fasta_section = False
    fasta_lines = []
    gff_version_line = ''

    # Read the GFF3 file and collect all features
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()

            if line.startswith('##gff-version'):
                gff_version_line = line
                continue
            if line.startswith('##sequence-region'):
                # Parse sequence-region directive to get sequence lengths
                parts = line.split()
                if len(parts) >= 4:
                    seqid = parts[1]
                    seq_start = int(parts[2])
                    seq_end = int(parts[3])
                    seqid_lengths[seqid] = seq_end
                continue
            if line == '##FASTA':
                fasta_section = True
                continue
            if fasta_section:
                fasta_lines.append(line)
                continue
            if not line or line.startswith('#'):
                continue

            fields = line.split('\t')
            if len(fields) != 9:
                continue  # Skip malformed lines

            seqid, source, feature_type, start, end, score, strand, phase, attributes_str = fields

            feature = {
                'seqid': seqid,
                'source': source,
                'type': feature_type,
                'start': int(start),
                'end': int(end),
                'score': score,
                'strand': strand,
                'phase': phase,
                'attributes': attributes_str,
                'line': line
            }

            features_by_seqid[seqid].append(feature)

    # Now process the features to add intergenic regions, start/stop codons, and introns
    output_features = []

    for seqid, features in features_by_seqid.items():
        # Sort features by start coordinate
        features.sort(key=lambda x: x['start'])

        # Collect gene/transcript features to identify exons
        exons_by_parent = defaultdict(list)

        for feature in features:
            if feature['type'] == 'exon':
                attributes = parse_attributes(feature['attributes'])
                parent_id = attributes.get('Parent')
                if parent_id:
                    exons_by_parent[parent_id].append(feature)

        # Generate start codons, stop codons, and introns for each transcript
        gene_regions = []

        for parent_id, exons in exons_by_parent.items():
            # Sort exons by genomic coordinates
            exons.sort(key=lambda x: x['start'])
            strand = exons[0]['strand']
            source = exons[0]['source']

            # Determine start codon and stop codon positions
            if strand == '+':
                # Start codon at start of first exon
                start_codon_start = exons[0]['start']
                start_codon_end = start_codon_start + 2  # Start codon is 3 bases
                # Stop codon at end of last exon
                stop_codon_end = exons[-1]['end']
                stop_codon_start = stop_codon_end - 2
            elif strand == '-':
                # Start codon at end of first exon
                start_codon_end = exons[0]['end']
                start_codon_start = start_codon_end - 2
                # Stop codon at start of last exon
                stop_codon_start = exons[-1]['start']
                stop_codon_end = stop_codon_start + 2
            else:
                continue  # Skip if strand information is missing

            # Add start codon feature
            start_codon_attributes = f"ID=start_codon:{parent_id};Parent={parent_id}"
            start_codon_feature = {
                'seqid': seqid,
                'source': source,
                'type': 'start_codon',
                'start': start_codon_start,
                'end': start_codon_end,
                'score': '.',
                'strand': strand,
                'phase': '0',
                'attributes': start_codon_attributes
            }
            output_features.append(start_codon_feature)
            gene_regions.append((start_codon_start, start_codon_end))

            # Add stop codon feature
            stop_codon_attributes = f"ID=stop_codon:{parent_id};Parent={parent_id}"
            stop_codon_feature = {
                'seqid': seqid,
                'source': source,
                'type': 'stop_codon',
                'start': stop_codon_start,
                'end': stop_codon_end,
                'score': '.',
                'strand': strand,
                'phase': '0',
                'attributes': stop_codon_attributes
            }
            output_features.append(stop_codon_feature)
            gene_regions.append((stop_codon_start, stop_codon_end))

            # Add intron features
            for i in range(len(exons) - 1):
                intron_start = exons[i]['end'] + 1
                intron_end = exons[i + 1]['start'] - 1
                if intron_start <= intron_end:
                    intron_attributes = f"ID=intron:{parent_id}:{i+1};Parent={parent_id}"
                    intron_feature = {
                        'seqid': seqid,
                        'source': source,
                        'type': 'intron',
                        'start': intron_start,
                        'end': intron_end,
                        'score': '.',
                        'strand': strand,
                        'phase': '.',
                        'attributes': intron_attributes
                    }
                    output_features.append(intron_feature)
                    gene_regions.append((intron_start, intron_end))

            # Collect exon regions
            for exon in exons:
                gene_regions.append((exon['start'], exon['end']))

        # Merge gene regions to find intergenic regions
        gene_regions.sort()
        merged_gene_regions = []
        for start, end in gene_regions:
            if not merged_gene_regions:
                merged_gene_regions.append([start, end])
            else:
                last_start, last_end = merged_gene_regions[-1]
                if start <= last_end + 1:
                    merged_gene_regions[-1][1] = max(last_end, end)
                else:
                    merged_gene_regions.append([start, end])

        # Determine sequence length
        seq_length = seqid_lengths.get(seqid)
        if not seq_length:
            # If sequence length is not provided, estimate from features
            seq_length = max(feature['end'] for feature in features)

        # Find intergenic regions
        prev_end = 1
        for region_start, region_end in merged_gene_regions:
            if prev_end < region_start:
                intergenic_start = prev_end
                intergenic_end = region_start - 1
                intergenic_attributes = f"ID=intergenic_region:{seqid}:{intergenic_start}-{intergenic_end}"
                intergenic_feature = {
                    'seqid': seqid,
                    'source': '.',
                    'type': 'intergenic_region',
                    'start': intergenic_start,
                    'end': intergenic_end,
                    'score': '.',
                    'strand': '.',
                    'phase': '.',
                    'attributes': intergenic_attributes
                }
                output_features.append(intergenic_feature)
            prev_end = region_end + 1
        if prev_end <= seq_length:
            # Add intergenic region after the last gene region
            intergenic_start = prev_end
            intergenic_end = seq_length
            intergenic_attributes = f"ID=intergenic_region:{seqid}:{intergenic_start}-{intergenic_end}"
            intergenic_feature = {
                'seqid': seqid,
                'source': '.',
                'type': 'intergenic_region',
                'start': intergenic_start,
                'end': intergenic_end,
                'score': '.',
                'strand': '.',
                'phase': '.',
                'attributes': intergenic_attributes
            }
            output_features.append(intergenic_feature)

    # Now that we've added intergenic regions, start/stop codons, and introns,
    # we proceed to remove all features that are not exons, start/stop codons, introns, or intergenic regions.

    # Collect the features to keep
    features_to_keep = []

    for seqid, features in features_by_seqid.items():
        for feature in features:
            if feature['type'] in {'exon'}:
                features_to_keep.append(feature)

    # Combine features from output_features (which includes the new features added)
    features_to_keep.extend(output_features)

    # Filter out any features that are not the desired types
    desired_feature_types = {'exon', 'start_codon', 'stop_codon', 'intron', 'intergenic_region'}
    filtered_features = [f for f in features_to_keep if f['type'] in desired_feature_types]

    # Sort features by seqid and start position
    filtered_features.sort(key=lambda x: (x['seqid'], x['start']))

    # Write the output GFF3 file
    with open(output_file, 'w') as outfile:
        if gff_version_line:
            outfile.write(f"{gff_version_line}\n")

        # Optionally, write sequence-region directives
        for seqid, length in seqid_lengths.items():
            outfile.write(f"##sequence-region {seqid} 1 {length}\n")

        for feature in filtered_features:
            fields = [
                feature['seqid'],
                feature['source'],
                feature['type'],
                str(feature['start']),
                str(feature['end']),
                feature['score'],
                feature['strand'],
                feature['phase'],
                feature['attributes']
            ]
            outfile.write('\t'.join(fields) + '\n')

        # Write any FASTA sequences if present
        if fasta_lines:
            outfile.write('##FASTA\n')
            for line in fasta_lines:
                outfile.write(line + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python modify_gff3.py input.gff3 output.gff3")
        sys.exit(1)

    input_gff3 = sys.argv[1]
    output_gff3 = sys.argv[2]
    process_gff3(input_gff3, output_gff3)