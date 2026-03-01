"""
Toolkit for biological sequences transformations and FASTQ filtering.

- Biological sequences slasses: Objects for transformations and proteolytic analysis of molecules.
- FASTQ processing: Filtering by GC-content, length, and quality.
"""

import os

from statistics import mean
from typing import Tuple, Union, List
from abc import ABC, abstractmethod
from collections.abc import Sequence, Number

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


class BiologicalSequence(ABC):
    """
    A base class representing a generic biological sequence. This class is abstract.
    """

    def __init__(self, sequence: str):
        self.sequence = sequence
        if not self._is_valid_alphabet():
            raise ValueError(
                f"Sequence {self.sequence} is not a valid {self.__class__.__name__} sequence."
            )

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        item = self.sequence[index]
        return self.__class__(item)

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod  # This was unnecessary, since it creates code duplication in subclasses, but I am obliged to do it by the task
    def _is_valid_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    A class representing a nucleic acid sequence (DNA or RNA) with various utility methods.
    """

    _ALPHABET = set()
    _COMPLEMENT_MAP = {}

    def _is_valid_alphabet(self) -> bool:
        if not self._ALPHABET:
            raise NotImplementedError("Subclass must define alphabet.")
        return set(self.sequence) <= self._ALPHABET

    def reverse(self) -> str:
        """
        Reverses the input nucleic acid sequence.

        Arguments:
        seq: str - Input DNA or RNA sequence

        Returns:
        str - Reversed sequence
        """
        return self.__class__(self.sequence[::-1])

    def complement(self) -> str:
        """
        Returns the complementary nucleic acid sequence.

        Arguments:
        seq: str - Input sequence

        Returns:
        str - Complementary sequence
        """
        if not self._COMPLEMENT_MAP:
            raise NotImplementedError("Subclasses must define complement map.")

        complement_seq = "".join(self._COMPLEMENT_MAP[char] for char in self.sequence)

        return self.__class__(complement_seq)

    def reverse_complement(self) -> str:
        """
        Returns the reverse complementary nucleic acid sequence.

        Arguments:
        seq: str - Input nucleic acid sequence

        Returns:
        str - Reverse complementary nucleic acid sequence
        """
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """
    A class representing a DNA sequence.
    """

    _ALPHABET = set("ATCGatcg")
    _COMPLEMENT_MAP = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }

    def transcribe(self) -> str:
        """
        Transcribes a DNA sequence into RNA (replaces T with U).

        Arguments:
        seq: str - Input DNA sequence

        Returns:
        str - Transcribed RNA sequence
        """
        transcribed_seq = self.sequence.replace("T", "U").replace("t", "u")
        return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence):
    """
    A class representing an RNA sequence.
    """

    _ALPHABET = set("AUCGaucg")
    _COMPLEMENT_MAP = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }


class AminoAcidSequence(BiologicalSequence):
    """
    A class representing an amino acid sequence.
    """

    _ALPHABET = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    def _is_valid_alphabet(self) -> bool:
        return set(self.sequence) <= self._ALPHABET

    def trypsin_sites(self) -> List[int]:
        """
        Returns the indices of trypsin cleavage sites in the amino acid sequence.

        Returns:
        List[int] - A list of indices where trypsin would cleave the sequence (0-based).
        """
        cleavage_sites = []
        for i in range(len(self.sequence) - 1):
            if self.sequence[i] in "KRkr" and self.sequence[i + 1] != "Pp":
                cleavage_sites.append(i)
        return cleavage_sites


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filters FASTQ records based on GC-content, length, and quality score, and saves the filtered records to a new output file.
    Uses Biopython's SeqIO and SeqUtils for file processing.

    Arguments:
    input_fastq: str - Path to the input FASTQ file.
    output_fastq: str - Name of the output file (will be saved in 'filtered/' directory).
    gc_bounds: tuple or float - Interval (min, max) or upper bound (max) for GC content (default (0, 100)).
    length_bounds: tuple or int - Interval or upper bound (max) for sequence length (default (0, 2**32)).
    quality_threshold: float - Minimum average Phred33 quality score (default 0).
    """
    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)

    def is_in_range(value, bounds):
        if isinstance(bounds, Number):
            return value <= bounds
        elif isinstance(bounds, Sequence) and not isinstance(bounds, str):
            return bounds[0] <= value <= bounds[1]
        else:
            raise TypeError("Invalid type for bounds")

    def filtered_records():
        with open(input_fastq, "r") as handle:
            for record in SeqIO.QualityIO.FastqPhredIterator(handle):
                # Check length
                if not is_in_range(len(record.seq), length_bounds):
                    continue

                # Check GC content
                gc_val = gc_fraction(record.seq) * 100
                if not is_in_range(gc_val, gc_bounds):
                    continue

                # Check quality
                qualities = record.letter_annotations["phred_quality"]
                avg_quality = mean(qualities) if qualities else -1
                if avg_quality < quality_threshold:
                    continue

                # As recommended in the https://biopython.org/docs/latest/api/Bio.SeqIO.QualityIO.html#Bio.SeqIO.QualityIO.FastqPhredIterator documentation,
                # I made a generator using yield
                yield record

    with open(output_path, "w") as out_handle:
        SeqIO.write(filtered_records(), out_handle, "fastq")
