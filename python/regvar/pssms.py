"""
Choose PSSMs to use in analysis
"""

import biopsy.transfac as tr


def get_pssms():
    pssm_filter = tr.PssmFilter(
        use_consensus_sequences=False,
        species_filter='V')
    return tr.Matrix.pssms(pssm_filter)
