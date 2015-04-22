#!/usr/bin/env python

import pandas as pd
import numpy as np

#
# Load the CSV file
S1 = pd.read_csv('S1.csv', header=6).ix[1:]

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

#
# Show the data
print S1.head()
print S1.columns

# S1.Position.apply(int)
for row in S1.iterrows():
    print row[1].Position
