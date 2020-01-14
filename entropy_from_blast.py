# Jonathan Oribello, 2020-01-13
"""
Calculates Shannon entropy from sequence homology using blast files.
"""
import argparse

from Bio import SearchIO
from Bio import SeqIO

def configure():
    """ configure variables from command line and defaults """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-blast_input", "-b", help="input file, must be xml format of blast alignment")
    parser.add_argument("-fasta_input", "-f", help="input file, must be xml format of blast alignment")
    parser.add_argument("-output", "-o", help="output file, csv format")
    args = parser.parse_args()
    return args

def main():
    args = configure()

    blasts = SearchIO.parse(args.blast_input, "blast-xml")
    fastas = SeqIO.parse(args.fasta_input, "fasta")

    for sequence in fastas:
        blast = next(blasts)
        if blast.id != sequence.id:
            raise ValueError(f"Unexpected mismatch between fasta and blast file: {blast.id} vs {sequence.id}")
        compute_entropy(sequence, blast)

def compute_entropy(sequence: SeqIO.SeqRecord, blast: SearchIO.QueryResult):
    """

    Args:
        sequence: primary sequence information of query
        blast: alignment results for query
    """
    # TODO: implement next
    pass


if __name__ == '__main__':
    main()