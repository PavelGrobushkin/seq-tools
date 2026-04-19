# seq-tools 🧬

## What is it?

**seq-tools** is a Python package for biological sequence manipulation and NGS data filtration.

It can be used for:

-   DNA, RNA, and Protein sequence processing (create reverse, complement, or reverse-complement counterparts for nucleic acids, as well as identification of trypsin cleavage sites in a protein).

-   performing filtration of NGS reads in .fastq-files based on read length, GC content, and quality;

-   сonverting a multi-line FASTA file (where sequence may be split across multiple lines) into a single-line FASTA file (where each sequence is stored on one line);

-   parsing the standard BLAST output file (txt), extracting the description of the best match for each query.

It is an an educational project done as a part of Python course at the [Bioinformatics institute](https://bioinf.me/en) (2025-2026 cohort).

## Features

### 1. Sequence Classes

A class hierarchy with automated alphabet validation:

-   `DNASequence`: Methods for reverse complement and transcription.

-   `RNASequence`: Methods for RNA complementation.

-   `AminoAcidSequence`: Trypsin cleavage site identification.

### 2. FASTQ Filter (`filter_fastq` function)

Filters FASTQ files using Biopython's `FastqPhredIterator` for iterative record processing without loading entire files into memory.

| Argument | Default | Description |
|------------------------|------------------------|------------------------|
| `input_fastq` | `Required` | Path to input FASTQ file. |
| `output_fastq` | `Required` | Output filename (saved in `filtered/` directory). |
| `gc_bounds` | `(0, 100)` | Range `(min, max)` or upper threshold for GC content %. |
| `length_bounds` | `(0, 2**32)` | Range `(min, max)` or upper threshold for sequence length. |
| `quality_threshold` | `0` | Minimum mean **Phred33** score per read. |

### 3. Bio Files Processor

File format utilities:

-   `convert_multiline_fasta_to_oneline`: Converts a multi-line FASTA file (where sequences are split across lines) into a single-line format.

-   `parse_blast_output`: Parses a standard BLAST result file (.txt) to extract and sort the descriptions of the best matches (first hit) for each query.

## Usage Example

``` python
from seqtools import DNASequence, AminoAcidSequence, filter_fastq
import bio_files_processor as bfp

# Sequence manipulation
dna = DNASequence("ATGGC")
print(dna.reverse_complement()) # GCCAT
print(dna.transcribe())         # RNASequence('AUGGC')

prot = AminoAcidSequence("MKKRRPL")
print(prot.trypsin_sites())     # [1, 2, 3]

# FASTQ Filter - retain only reads with GC content between 40% and 60%, no longer than 90 nucleotides and with average quality > 20
filtered_seqs = st.filter_fastq(
    input_fastq="data.fastq",
    output_fastq="clean_data.fastq" 
    gc_bounds=(40, 60),
    length_bounds = 90,
    quality_threshold=20
)

# Bio Files Processor
# Convert FASTA
bfp.convert_multiline_fasta_to_oneline("input_multiline.fasta", "output_oneline.fasta")

# Parse BLAST results
bfp.parse_blast_output("blast_results.txt", "sorted_descriptions.txt")
```

## Project Structure

```         
.
├── README.md
├── seqtools.py              # Sequence classes and FASTQ filter
└── bio_files_processor.py   # FASTA and BLAST utilities
```