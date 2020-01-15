# Jonathan Oribello, 2020-01-13
"""
Calculates Shannon entropy from sequence homology using blast files.
"""
import argparse
from collections import namedtuple, Counter
from itertools import chain
import math

from Bio import SearchIO
from Bio import SeqIO

def configure():
    """ configure variables from command line and defaults """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-blast-input", "-b", help="input file, must be xml format of blast alignment")
    parser.add_argument("-fasta-input", "-f", help="input file, must be xml format of blast alignment")
    parser.add_argument("-output", "-o", help="output file, csv format")
    parser.add_argument("--save-alignments", "-s", dest="save_alignments", action="store_true", help="Save alignments to csv format, not necessary "
                                                                        "for calculation of entropy",)
    parser.add_argument("--no-save-alignments", dest="save_alignments", action="store_false")
    parser.set_defaults(save_alignments=False)
    args = parser.parse_args()
    return args

def main():
    args = configure()

    blasts = SearchIO.parse(args.blast_input, "blast-xml")
    fastas = SeqIO.parse(args.fasta_input, "fasta")

    for sequence in fastas:
        pdb5id = ''.join(sequence.id.split('|')[1:3])
        print(f"Computing entropy for {pdb5id}")
        blast = next(blasts)
        if blast.id != sequence.id:
            raise ValueError(f"Unexpected mismatch between fasta and blast file: {blast.id} vs {sequence.id}")
        alignment_arrays, bit_scores = extract_aligned_sequences(sequence, blast)
        if args.save_alignments:
            write_alignments(alignment_arrays, bit_scores, filename=f"test/alignments/{pdb5id}_alignments.csv")
        entropies = compute_entropy(alignment_arrays, bit_scores, bit_score_cutoff_percent=40)
        write_entropies(entropies, filename=f"test/40_percent_cutoff/{pdb5id}_entropy.csv")

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
        f.write(f"{','.join([str(score) for score in bit_scores])}\n")
        for array in alignment_arrays:
            f.write(f"{array.index+1},{','.join(array.alignments)}\n")

def __recategorize(alignment_list, categories):
    category_list = alignment_list.copy()
    for i, category in enumerate(categories):
        category_list = [
            str(i) if (residue in category)
            else residue
            for residue in category_list
        ]
    return category_list

def __filter_aligned(alignment_list, allowed):
    return [
        residue
        for residue in alignment_list
        if residue in allowed
    ]

def __calculate_entropy(sequence):
    counts = Counter(sequence)

    entropy = 0
    for residue, count in counts.items():
        prob = float(count / len(sequence))
        entropy = entropy - (prob * (math.log(prob) / math.log(2)))
    return entropy

E6_CATEGORIES = [list("AVLIMC"), # aliphatic
                  list("FWYH"), # aromatic
                  list("STNQ"), # polor
                  list("KR"), # positive
                  list("DE"), # negative
                  list("GP")] # special

E20_CATEGORIES = list("AVLIMCFWYHSTNQKRDEGP") # all 20 canonical

def compute_entropy(alignment_arrays, bit_scores, bit_score_cutoff_percent, schemes={"E20":E20_CATEGORIES,"E6":E6_CATEGORIES}):
    Entropy = namedtuple("Entropy", field_names=["index","residue","bit_score_cutoff_percent","entropy"])
    entropies = list()
    for array in alignment_arrays:
        minimum_score = bit_scores[0] * bit_score_cutoff_percent / 100
        aligned_with_minimum_score = [
            residue for i, residue in enumerate(array.alignments)
            if bit_scores[i] >= minimum_score
        ]

        entropy = Entropy(index=array.index,
                          residue=array.residue,
                          entropy=dict(),
                          bit_score_cutoff_percent=bit_score_cutoff_percent)
        for scheme, categories in schemes.items():
            if scheme == "E20":

                calculation_sequence = __filter_aligned(alignment_list=aligned_with_minimum_score, allowed=list(chain(*categories)))
                calculation_sequence = __recategorize(calculation_sequence, categories=categories)
                entropy.entropy[scheme] = __calculate_entropy(calculation_sequence)
            elif scheme == "E6":
                calculation_sequence = __filter_aligned(alignment_list=aligned_with_minimum_score, allowed=list(chain(*categories)))
                calculation_sequence = __recategorize(calculation_sequence, categories=categories)
                entropy.entropy[scheme] = __calculate_entropy(calculation_sequence)

        entropies.append(entropy)

    return entropies

def write_entropies(entropies, filename):
    with open(filename, "w") as f:
        # header
        f.write(f"index,residue,{','.join(entropies[0].entropy.keys())}\n")
        for entropy in entropies:
            f.write(f"{entropy.index+1},{entropy.residue},{','.join([str(value) for value in entropy.entropy.values()])}\n")

if __name__ == '__main__':
    main()