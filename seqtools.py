"""
Sequencing tools

This module provides utility functions for basic bioinformatics tasks.
Use run_dna_rna_tools() for various transformations of nucleotide sequences. 
Use filter_fastq() to filter reads from FASTQ-files by length, GC content, and quality.

Requirements:
- typing (standard library)
- statistics (standard library)
- os (standard library)
"""

from typing import Dict, Tuple, Union, List
import os
from seqtools_modules.dna_rna_tools import *
from seqtools_modules.fastq_tools import *


def run_dna_rna_tools(*args: str) -> Union[str, bool, List[Union[str, bool]]]:
    """
    This function returns complementary, reversed, transcribed and complementary
    to reversed sequences or checks if a given sequence is DNA/RNA for an unlimited
    number of sequences given as arguments.

    It accepts as input a list of arguments, the last of which should be the name
    of one of the 5 following functions: "is_nucleic_acid", "reverse", "transcribe",
    "complement", "reverse_complement".

    Arguments:
    *args: str - Variable number of arguments.
    The last argument must be the command.

    Returns:
    str or bool or list - if multiple sequences are provided, returns a list of results.

    Raises:
    ValueError - If the command is unknown or if the sequences are invalid.
    """

    sequences = args[:-1]
    procedure = args[-1]
    functions = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    if procedure not in functions:
        raise ValueError(f"Unknown procedure: {procedure}. Please try again.")

    if len(sequences) > 1:
        return [functions[procedure](seq) for seq in sequences]
    else:
        return functions[procedure](args[0])


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2 ** 32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters FASTQ sequences in a on-the-fly manner from an input file based on GC-content,
    length, and quality score, and saves the filtered records to a new output file.

    Arguments:
    input_fastq: str - Path to the input FASTQ file.
    output_fastq: str - Name of the output file (will be saved in 'filtered/' directory).
    gc_bounds: tuple or float - Interval (min, max) or upper bound (max) for GC content (default (0, 100)).
    length_bounds: tuple or int - Interval or upper bound (max) for sequence length (default (0, 2**32)).
    quality_threshold: float - Minimum average Phred33 quality score (default 0).
    """
    if not os.path.isdir("filtered"):
        os.mkdir("filtered")

    with open(input_fastq, "r") as fastq_file, open(
        os.path.join('filtered', output_fastq), "w"
    ) as filtered_output:
        while True:
            head_line = fastq_file.readline()
            if not head_line:
                break
            seq_line = fastq_file.readline()
            sep_line = fastq_file.readline()
            qual_line = fastq_file.readline()

            sequence = seq_line.strip()
            quality = qual_line.strip()

            if check_filters(
                sequence, quality, quality_threshold, gc_bounds, length_bounds
            ):
                filtered_output.writelines([head_line, seq_line, sep_line, qual_line])