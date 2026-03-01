# seq-tools 🧬

## What is it?
**seq-tools** is a Python package containing bioinformatics tools that сan:
- check if a gived sequence is DNA or RNA and converts a DNA sequence into its reverse, complement, or reverse-complement counterpart;
- perform filtration of NGS reads in `.fastq`-files based on read length, GC content, and quality;
- сonvert a multi-line FASTA file (where sequence may be split across multiple lines) into a single-line FASTA file (where each sequence is stored on one line);
- parse the standard BLAST output file (txt), extracting the description of the best match for each query.

It is an an educational project done as a part of Python course in the [Bioinformatics Institute](https://bioinf.me/en) (2025-2026 cohort).

## Features

### 1. DNA/RNA Toolkit (`run_dna_rna_tools`)

This function executes standard molecular biology operations:

* `is_nucleic_acid`: Validates if a sequence is DNA or RNA.

* `reverse`: Reverses the sequence direction.

* `complement`: Finds the complementary strand (DNA only).

* `reverse_complement`: Finds the reverse complementary strand (DNA only).

* `transcribe`: Converts DNA to RNA (T → U).

### 2. FASTQ Filter (`filter_fastq`)

This tool processes FASTQ data directly from files. It filters reads "on-the-fly" (line-by-line) to handle large files efficiently without loading everything into RAM.

| Argument | Default Value | Description | 
| ----- | ----- | ----- | 
| `input_fastq` | `Required` | Path to the input FASTQ file.| 
| `output_fastq` | `Required` | Name of the output file. It will be created inside the `filtered/` directory. | 
| `gc_bounds` | `(0, 100)` | Range (tuple) or upper limit (float) for percentage of G/C content. | 
| `length_bounds` | `(0, 2**32)` | Range (tuple) or upper limit (int) for sequence length. | 
| `quality_threshold` | `0` | Minimum average **Phred33 Q-score** [(from 0 to 40)](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/QualityScoreEncoding_swBS.htm) required across the entire read. | 

### 3. Bio Files Processor (`bio_files_processor`)

A set of utilities for handling common bioinformatics file formats:

- `convert_multiline_fasta_to_oneline`: Converts a multi-line FASTA file (where sequences are split across lines) into a single-line format.
- `parse_blast_output`: Parses a standard BLAST result file (.txt) to extract and sort the descriptions of the best matches (first hit) for each query.

## Usage Example

```python
import seqtools as st
import bio_files_processor as bfp

# DNA/RNA Toolkit - perform reverse complement on two DNA sequences
result = st.run_dna_rna_tools("ATGGC", "GCt", "reverse_complement")
print(result) 
# Output:
#['GCCAT', 'aGC']

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
```text
.
├── README.md                # You are here
├── seqtools.py              # Main script for sequence manipulation and FASTQ filtering
├── bio_files_processor.py   # Script for FASTA and BLAST file processing
└── seqtools_modules/        
    ├── dna_rna_tools.py     # Modules for DNA/RNA sequence manipulation
    └── fastq_tools.py       # Modules for FASTQ read filtering
```
