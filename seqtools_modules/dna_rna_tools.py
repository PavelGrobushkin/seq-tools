rna_nucleotides = {"A", "a", "C", "c", "G", "g", "U", "u"}
dna_nucleotides = {"A", "a", "C", "c", "G", "g", "T", "t"}

complement_dna = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
}


def is_dna(seq: str):
    """
    Helper function that checks if sequence is valid DNA.
    Raises ValueError if not.
    """
    if not set(seq) <= dna_nucleotides:
        raise ValueError(f"Sequence {seq} is not a valid DNA sequence.")


def is_nucleic_acid(seq: str) -> bool:
    """
    Checks if the sequence consists only of valid DNA or RNA nucleotides.

    Arguments:
    seq: str - Input sequence

    Returns:
    bool - True if sequence is valid DNA or RNA, False otherwise.
    """
    return set(seq) <= rna_nucleotides or set(seq) <= dna_nucleotides


def transcribe(seq: str) -> str:
    """
    Transcribes a DNA sequence into RNA (replaces T with U).

    Arguments:
    seq: str - Input DNA sequence

    Returns:
    str - Transcribed RNA sequence

    Raises:
    ValueError - If the input sequence contains invalid characters or is RNA.
    """
    is_dna(seq)
    return seq.replace("T", "U").replace("t", "u")


def reverse(seq: str) -> str:
    """
    Reverses the input nucleic acid sequence.

    Arguments:
    seq: str - Input DNA or RNA sequence

    Returns:
    str - Reversed sequence

    Raises:
    ValueError - If the input sequence is not a nucleic acid.
    """
    if not is_nucleic_acid(seq):
        raise ValueError(f"Sequence {seq} is not a nucleic acid.")
    return seq[::-1]


def complement(seq: str) -> str:
    """
    Returns the complementary DNA sequence.

    Arguments:
    seq: str - Input DNA sequence

    Returns:
    str - Complementary DNA sequence

    Raises:
    ValueError - If the input sequence contains invalid characters or is RNA.
    """
    is_dna(seq)
    return "".join(complement_dna[char] for char in seq)


def reverse_complement(seq: str) -> str:
    """
    Returns the reverse complementary DNA sequence.

    Arguments:
    seq: str - Input DNA sequence

    Returns:
    str - Reverse complementary DNA sequence

    Raises:
    ValueError - If the input sequence contains invalid characters or is RNA.
    """
    return complement(reverse(seq))
