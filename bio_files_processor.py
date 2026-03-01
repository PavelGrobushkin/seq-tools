import os
from typing import Optional


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: Optional[str] = None
) -> None:
    """
    Converts a multi-line FASTA file (where sequence may be split across multiple lines)
    into a single-line FASTA file (where each sequence is stored on one line).

    Arguments:
    input_fasta: str - Path to the input multi-line FASTA file.
    output_fasta: str, optional - Name of the output FASTA file. If None, the filename
    will be automatically generated as <input_fasta>_onelined.fa

    Returns:
    None
    """
    if output_fasta is None:
        basename = os.path.splitext(input_fasta)[0]
        output_fasta = f"{basename}_onelined.fa"

    seq = ""

    with open(input_fasta, mode="r") as input_file, open(
        output_fasta, mode="w"
    ) as output_file:
        while True:
            line = input_file.readline()
            if not line:
                if seq:
                    output_file.write(seq + "\n")
                break
            elif line.startswith(">"):
                if seq:
                    output_file.write(seq + "\n")
                output_file.write(line)
                seq = ""
            else:
                seq += line.strip()


def parse_blast_output(input_file: str, output_file: Optional[str] = None) -> None:
    """
    Parses the standard BLAST output file (txt), extracting the description
    (Description) of the best match for each query.

    The extracted descriptions are sorted alphabetically and written
    to the specified output file.

    Arguments:
    input_file: str - Path to the BLAST results input file.
    output_file: str, optional - Path to the output file. If None, the filename
    will be automatically generated as <input_file>_best_sorted.txt.

    Returns:
    None
    """
    if output_file is None:
        basename = os.path.splitext(input_file)[0]
        output_file = f"{basename}_best_sorted.txt"
    best_prot_list = []
    with open(input_file, mode="r") as blast:
        while True:
            line = blast.readline()
            if not line:
                break
            elif line.startswith("Sequences producing significant alignments:\n"):
                line = blast.readline()
                position = line.find("S")
                for i in range(2):
                    line = blast.readline()
                best_prot_list.append(line[:position].strip() + "\n")
    best_prot_list.sort()
    with open(output_file, mode="w") as output:
        output.writelines(best_prot_list)
