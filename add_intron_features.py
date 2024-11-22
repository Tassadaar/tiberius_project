#!/usr/bin/env python

"""
Copyright 2023 Joran Martijn.

The development of this and other bioinformatic tools in the Roger lab was funded by
Discovery grant RGPIN-2022-05430 from the Natural Sciences and Engineering Research Council of Canada
awarded to Andrew J. Roger.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.
If not, see <https://www.gnu.org/licenses/>.
"""

# default libraries
import argparse, sys

# stuff you need to install
import gffutils
import sqlite3

# check if you have version 0.12 or later
from packaging.version import Version
if Version(gffutils.__version__) < Version("0.12"):
    print('For this script, gffutils needs to be at least version 0.12')
    sys.exit(1)

# NOTE:
#   The script assumes that exon features have an ID 
#   that includes the word 'exon'

# make command line interface
parser = argparse.ArgumentParser(description="Add intron features to a GFF3 file")

parser.add_argument(
    "-g", 
    dest="gff3_file",
    metavar="GFF3_FILE",
    required=True, 
    help="Input GFF3 file"
)

args = parser.parse_args()


# custom function to ensure unique ID's in some cases
def id_func(f) -> str:
    '''
    If merge_strategy='create_unique' fails,
    use this function to create unique id's
    as a fallback strategy.

    Particularly useful for Liftoff-generated
    GFF3's, which like gffutils also appends
    _1 to identical feature IDs
    '''
    if f.featuretype in ['gene', 'mRNA']:
        return f.attributes['ID'][0]
    elif f.featuretype in ['exon', 'CDS']:
        return '-'.join( [f.attributes['ID'][0], f.seqid, str(f.start), str(f.end)] )


def main(args) -> None:

    # load GFF3
    try:
        db = gffutils.create_db(args.gff3_file, dbfn=':memory:', merge_strategy='create_unique')
    except sqlite3.IntegrityError:
        print( "Default loading strategy 'create_unique' failed,",
               "changing to custom id_spec function strategy",
              file=sys.stderr )
        db = gffutils.create_db(args.gff3_file, dbfn=':memory:', id_spec=id_func)

    # if the GFF3 has some introns explicitly mentioned,
    # but not all, it's easiest to just purge them
    # and let gffutils recreate them
    for i in db.features_of_type('intron'):
        db.delete(i)

    # infer and add intron features
    introns = list(db.create_introns())

    prev_mrna = ''
    for i in introns:

        # change source field from 'gffutils_derived' to 'gffutils'
        i.source = 'gffutils'

        # if exons don't have IDs, resultant introns won't either,
        # so we'll have to create them here
        if not 'ID' in i.attributes:
            mrna = i.attributes['Parent'][0]
            if prev_mrna != mrna:
                prev_mrna = mrna
                n = 1
            else:
                n += 1
            i.attributes['ID'] = [f'{mrna}.exon{n}',f'{mrna}.exon{n+1}']


    # update db with new introns
    db.update(introns)

    # print updated gff record in logical order
    for g in db.features_of_type('gene', order_by=('seqid', 'start')):
        print()
        print(g)
        for m in db.children(g, featuretype=('mRNA','tRNA')):
            print(m)
            for f in db.children(m, order_by='start'):
                print(f)



# # we need a natural sort function
# # to sort exon10,exon9 -> exon9,exon10
# # I don't understand it, got it from StackOverflow
# def natural_sort_key(s):
#     _nsre = re.compile('([0-9]+)')
#     return [int(text) if text.isdigit() else text.lower()
#             for text in re.split(_nsre, s)]



if __name__ == '__main__':
    main(args)
