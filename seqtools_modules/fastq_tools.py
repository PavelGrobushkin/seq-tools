from typing import Union, Tuple
from statistics import mean


def quality_filter(phred33: str, quality_threshold: float) -> bool:
    """
    Checks if the average quality of the read is above the threshold.

    Arguments:
    phred33: str - Quality string encoded in Phred33 scale.
    quality_threshold: float - Minimum average quality score required.

    Returns:
    bool - True if average quality >= threshold, False otherwise.
    """
    return quality_threshold <= mean(ord(char) - 33 for char in phred33)


def gc_filter(seq: str, gc_bounds: Union[Tuple[float, float], float]) -> bool:
    """
    Checks if the GC-content of the sequence is within the specified bounds.

    Arguments:
    seq: str - DNA sequence.
    gc_bounds: tuple or float - Interval (min, max) or upper bound (max) for GC content (%).

    Returns:
    bool - True if GC content is within bounds, False otherwise.
    """
    gc_content = 100 * (seq.count("C") + seq.count("G")) / len(seq)
    if isinstance(gc_bounds, tuple):
        return gc_bounds[0] <= gc_content <= gc_bounds[1]
    else:
        return gc_content <= gc_bounds


def len_filter(seq: str, length_bounds: Union[Tuple[int, int], int]) -> bool:
    """
    Checks if the length of the sequence is within the specified bounds.

    Arguments:
    seq: str - DNA sequence.
    length_bounds: tuple or int - Interval (min, max) or upper bound (max) for length.

    Returns:
    bool - True if length is within bounds, False otherwise.
    """
    length = len(seq)
    if isinstance(length_bounds, tuple):
        return length_bounds[0] <= length <= length_bounds[1]
    else:
        return length <= length_bounds
