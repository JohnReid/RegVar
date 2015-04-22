#!/usr/bin/env python

import pandas as pd
import numpy as np
import requests
import cStringIO
from bx.align import maf
from collections import defaultdict


def parseS1(csvfile):
    #
    # Load the CSV file
    S1 = pd.read_csv(csvfile, header=6).ix[1:]

    #
    # Rename the columns
    S1.columns = [
        u'Empty',
        u'TagSNP',
        u'Gene',
        u'Ref',
        u'ProxySNP',
        u'r2',
        u'D',
        u'Chr',
        u'Position',
        u'Risk.allele',
        u'Non.risk.allele',
        u'Risk.allele.freq',
        u'Type',
        u'Nearest.Genes',
        u'Omega.TFBS',
        u'p.est.TFBS',
        u'Omega.modules',
        u'p.est.modules',
        u'Omega.TFBS.in.modules',
        u'p.est.TFBS.in.modules',
        u'S.all',
        u'Omega.restr.TFBS',
        u'p.est.restr.TFBS',
        u'Omega.restr.modules',
        u'p.est.restr.modules',
        u'Omega.restr.TFBS.in.modules',
        u'p.est.restr.TFBS.in.modules',
        u'S.restr.all',
        u'PMCA.result'
    ]

    #
    # Remove empty first column
    S1 = S1.drop(u'Empty', axis=1)

    #
    # Drop empty rows
    S1 = S1.loc[~ np.isnan(S1.Position)]

    #
    # Retype data
    S1.Position = S1.Position.apply(int)

    return S1


def get_sequences_from_maf(reader, centre_species):
    """Extract the sequences from the MAF file."""
    # get sequences from MAF
    offset = None
    sequences = defaultdict(str)
    last_start = None
    for align in reader:
        # check alignments are increasing in position
        this_start = align.components[0].start
        if None != last_start and this_start < last_start:
            raise RuntimeError('Expecting alignments to increase in position.')
        last_start = this_start
        for comp in align.components:
            species, chromosome = comp.src.split('.')
            if None == offset and centre_species == species:
                offset = comp.start
            sequences[species] += comp.text

    # we must have got the offset
    if None == offset:
        raise RuntimeError('Never found offset for centre species.')

    return sequences, offset


#
# Parse Table S1 from the Claussnitzer paper
S1 = parseS1('S1.csv')

#
# Show the data
print S1.head()
print S1.columns

#
# Iterate through each SNP
alignmentsURL = 'http://localhost:9083/alignment'
# genome = 'hg19'
genome = 'hg38'
alignment = 'multiz20way'
spread = 60
for row_it in S1.iterrows():
    row = row_it[1]
    print row.ProxySNP, row.Chr, row.Position
    url = '{0}/{1}/{2}/{3}/{4}/{5}'.format(alignmentsURL, genome, alignment,
                                           row.Chr, row.Position - spread,
                                           row.Position + spread)
    r = requests.get(url)
    reader = maf.Reader(cStringIO.StringIO(r.text))
    sequences, offset = get_sequences_from_maf(reader, genome)
    len(sequences)
    offset
    print sequences[genome]
    break
