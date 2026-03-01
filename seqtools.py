"""
Sequencing tools

This module provides utility functions for basic bioinformatics tasks.
Use run_dna_rna_tools() for various transformations of nucleotide sequences. 
Use filter_fastq() to filter FASTQ reads by length, GC content, and quality.

Requirements:
- typing (standard library)
- statistics (standard library)
"""

from typing import Dict, Tuple, Union, List
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
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> Dict[str, Tuple[str, str]]:
    """
    Filters FASTQ sequences based on GC-content, length, and quality score.

    Arguments:
    seqs: dict - Dictionary of  FASTQ sequences {name: (sequence, quality_string)}.
    gc_bounds: tuple or float - Interval or upper bound for GC content (default (0, 100)).
    length_bounds: tuple or int - Interval or upper bound for sequence length (default (0, 2**32)).
    quality_threshold: float - Minimum average Phred33 quality score (default 0).

    Returns:
    dict - Filtered dictionary of sequences.
    """
    output = {}
    for read_id, (seq, qual) in seqs.items():
        if (
            quality_filter(qual, quality_threshold)
            and gc_filter(seq, gc_bounds)
            and len_filter(seq, length_bounds)
        ):
            output[read_id] = (seq, qual)

    return output
