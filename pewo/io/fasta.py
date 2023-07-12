"""
A module to work with .fasta-formatted files.
"""

__author__ = "Nikolai Romashchenko"
__license__ = "MIT"


import os
from string import Formatter
from Bio import SeqIO
from typing import List


def _seq_id_filter(id: str) -> str:
    """
    Replaces underscores and semicolons with dashes in the sequence IDs.
    This is needed to have nice output filename templates with underscore
    as a delimiter for parameters
    """
    result = id.replace("_", "-")
    return result.replace(";", "-")


def get_sequence_ids(input_file: str) -> List[str]:
    """
    Retrieves sequence IDs from the input .fasta file.
    """
    return [_seq_id_filter(record.id) for record in SeqIO.parse(input_file, "fasta")]


def _write_fasta(records: List[SeqIO.SeqRecord], filename: str) -> None:
    """
    Writes a list of fasta records to file.
    """
    with open(filename, "w") as output:
        SeqIO.write(records, output, "fasta")


def split_fasta(input_file: str, nametemplate: str) -> List[str]:
    """
    Splits the input .fasta file into multiple .fasta files,
    one sequence per file. Returns the list of resulting files.
    """
    files = []
    for record in SeqIO.parse(input_file, "fasta"):
        #FIXME:
        # By the convention, _r0 means "variable read length". Adding this
        # makes implicit dependency on the read file name convention in ALL rules
        # looking for read files: alignment_hmm, placement_rappas_dbinram etc.

            #Fill the in the name template with record id and 0 read length,
            #no matter what the keyworkd arguments are:
        keywords = [x[1] for x in Formatter().parse(nametemplate) if x[1]]
        replace = {key: val for key, val in zip(keywords,
                                                (_seq_id_filter(record.id), 0))}
        filename = nametemplate.format(**replace)

        _write_fasta([record], filename)
        files.append(filename)

    return files
