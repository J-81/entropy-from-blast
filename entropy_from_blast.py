# Jonathan Oribello, 2020-01-13
"""
Calculates Shannon entropy from sequence homology using blast files.
"""
import argparse
from collections import namedtuple

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
        alignment_arrays, bit_scores = extract_aligned_sequences(sequence, blast)
        write_alignments(alignment_arrays, bit_scores, filename=f"test/{''.join(sequence.id.split('|')[1:3])}.alignments")
        raise Exception
        compute_entropy()

def extract_aligned_sequences(sequence, blast):
    """ Extracts aligned letters based for each position in the query sequence

    Args:
        sequence: Bio.SeqRecord.SeqRecord
        blast: Bio.SearchIO._model.query.QueryResult

    Returns:
        alignment_arrays: list of named tuples for each query sequence position
        bit_scores: list of float, the bit score for each alignment recorded
    """
    AlignmentArray = namedtuple("AlignmentArray",
                                 field_names=["index",
                                              "residue",
                                              "alignments"]
                                 )
    # initialize containers for alignments
    alignment_arrays = list()
    for i, letter in enumerate(sequence):
        alignment_arrays.append(AlignmentArray(index=i,
                                               residue=letter,
                                               alignments=list()))
    # populate containers by iterating through each alignment
    bit_scores = list()
    for hit in blast:
        for hsp in hit:
            bit_scores.append(hsp.bitscore)
            aligned = dict()
            for i, letter in enumerate(hsp.hit):
                aligned[i+hsp.query_start] = letter
            for array in alignment_arrays:
                try:
                    array.alignments.append(aligned[array.index])
                except KeyError:
                    array.alignments.append(".") # * will denote not present in the hsp

    return alignment_arrays, bit_scores

def write_alignments(alignment_arrays, bit_scores, filename):
    """ Saves alignment data to file.

    Note: This is not needed for computing entropy, but is nice to have.

    Args:
        alignment_arrays: AlignmentArray, namedtuple
        bit_scores: [float] for bit_scores
        filename: filepath to save to
    """
    with open(filename, "w") as f:
        f.write(",".join([str(score) for score in bit_scores]))
        for array in alignment_arrays:
            f.write(f"{array.index},{','.join(array.alignments)}")


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