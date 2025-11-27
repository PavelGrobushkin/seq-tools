# seq-tools 🧬

## What is it?
**seq-tools** is a Python package containing bioinformatics tools that:
- checks if a gived sequence is DNA or RNA and converts a DNA sequence into its reverse, complement, or reverse-complement counterpart;
- performs filtration of FASTQ reads based on read length, GC content, and quality. 

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

This a tool for processing FASTQ data (passed as a dictionary of `{'read_id': ('sequence', 'quality')}`).

| Filter | Argument | Default Value | Description | 
| ----- | ----- | ----- | ----- | 
| **GC Content** | `gc_bounds` | `(0, 100)` | Range (tuple) or upper limit (float) for percentage of G/C content. | 
| **Length** | `length_bounds` | `(0, 2**32)` | Range (tuple) or upper limit (int) for sequence length. | 
| **Quality** | `quality_threshold` | `0` | Minimum average **Phred33 Q-score** [(from 0 to 40)](https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/QualityScoreEncoding_swBS.htm) required across the entire read. | 

## Usage Example

```python
import seqtools as sq

# DNA/RNA Toolkit - perform reverse complement on two DNA sequences
result = sq.run_dna_rna_tools("ATGGC", "GCt", "reverse_complement")
print(result) 
# Output:
#['GCCAT', 'aGC']

# FASTQ Filter - retain only reads with GC content between 40% and 60%, no longer than 90 nucleotides and with average quality > 20
filtered_seqs = sq.filter_fastq(
    seqs=EXAMPLE_FASTQ, # let's imagine that we have previously imported example data as ```from example_data import EXAMPLE_FASTQ``` 
    gc_bounds=(40, 60),
    length_bounds = 90,
    quality_threshold=20
)

print(len(filtered_seqs)})
# Possible output:
# 7
```

## Project Structure
```text
.
├── README.md                   # You are here
├── seqtools.py                 # Main entry point script (contains wrapping functions run_dna_rna_tools and filter_fastq)
└── seqtools_modules/           
    ├── dna_rna_tools.py        # Modules for DNA/RNA sequence manipulation
    └── fastq_tools.py          # Modules for FASTQ read filtering
```
