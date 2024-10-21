# Specimux
Dual barcode and primer demultiplexing for MinION sequenced reads

Specimux is a substantial revision and enhancement of minibar.py, originally developed by the
California Academy of Sciences, Institute for Biodiversity Science & Sustainability. 

Specimux is designed  to improve the accuracy and throughput of DNA barcode identification for
multiplexed MinION sequencing data, with a primary focus on serving the fungal sequencing
community.  Whereas minibar.py includes several processing methods supporting a variety of barcode
designs and matching regimes, specimux focuses specifically on high precision for demultiplexing
dual-indexed sequences.  

Developed as a volunteer contribution, Specimux attempts to fill a small niche in processing large-scale, community-driven
fungal DNA barcoding projects. 

The tool was developed and tested using the Mycomap ONT037 dataset, which comprises 768 specimens
and approximately 765,000 nanopore reads in FastQ format. This real-world dataset provided a robust
testing ground, ensuring Specimux's capabilities align closely with the needs of contemporary
fungal biodiversity research. By enhancing the speed and accuracy of sequence demultiplexing,
Specimux aims to allow human sequence validators to focus their efforts more effectively.

## Requirements
**specimux.py** is written in Python and is compatible with Python 3. It requires the following dependencies:
- edlib (install using `pip install edlib`)
- BioPython

You can install all required dependencies using:
```
pip install edlib biopython
```

Specimux has been tested on MacOS, and Linux machines.

## Usage
Here's a basic example of how to use Specimux:

```
python specimux.py barcode_file sequence_file
```

For a full list of options, use the `-h` or `--help` flag:

```
python specimux.py -h
```

Here are the main options available:

```
usage: specimux.py [-h] [--min-length MIN_LENGTH] [--max-length MAX_LENGTH]
                   [-n NUM_SEQS] [-p PERCENT_MATCH] [-e INDEX_EDIT_DISTANCE]
                   [-E PRIMER_EDIT_DISTANCE] [-l SEARCH_LEN]
                   [-A AMBIGUITY_THRESHOLD] [-F]
                   [-P OUTPUT_FILE_PREFIX] [-D OUTPUT_DIR] [--color]
                   [--trim {none,tails,primers}] [--diagnostics] [--debug]
                   [-t THREADS] [-v]
                   barcode_file sequence_file

Specimux: Demultiplex MinION sequence by dual barcode indexes and primers.

positional arguments:
  barcode_file          File containing barcode information
  sequence_file         Sequence file in Fasta or Fastq format, gzipped or
                        plain text

options:
  -h, --help            show this help message and exit
  --min-length MIN_LENGTH
                        Minimum sequence length. Shorter sequences will be
                        skipped (default: 400)
  --max-length MAX_LENGTH
                        Maximum sequence length. Longer sequences will be
                        skipped (default: 1400)
  -n NUM_SEQS, --num-seqs NUM_SEQS
                        Number of sequences to read from file (e.g., -n 100 or
                        -n 102,3)
  -p PERCENT_MATCH, --percent-match PERCENT_MATCH
                        Percentage match (default: 0.75)
  -e INDEX_EDIT_DISTANCE, --index-edit-distance INDEX_EDIT_DISTANCE
                        Barcode edit distance value, overrides -p
  -E PRIMER_EDIT_DISTANCE, --primer-edit-distance PRIMER_EDIT_DISTANCE
                        Primer edit distance value, overrides -p
  -l SEARCH_LEN, --search-len SEARCH_LEN
                        Length to search for index and primer at start and end
                        of sequence (default: 80)
  -A AMBIGUITY_THRESHOLD, --ambiguity-threshold AMBIGUITY_THRESHOLD
                        Threshold for considering the edit distance between
                        two barcodes to be different
  -F, --output-to-files
                        Create individual sample files for sequences
  -P OUTPUT_FILE_PREFIX, --output-file-prefix OUTPUT_FILE_PREFIX
                        Prefix for individual files when using -F (default:
                        sample_)
  -D OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory for individual files when using -F (default:
                        .)
  --color               Highlight barcode matches in blue, primer matches in
                        green
  --trim {none,tails,primers}
                        trimming to apply
  --diagnostics         Output extra diagnostics
  --debug               Enable debug logging
  -t THREADS, --threads THREADS
                        Number of worker threads to use
  -v, --version         show program's version number and exit
```

This provides a more comprehensive set of options compared to the original minibar.py, allowing for greater flexibility and control in the demultiplexing process.

## Major Enhancements in Specimux

Specimux introduces several changes compared to the original minibar.py program, focused on
maximizing the precision of specimen assignment, minimizing wrongly-assigned specimens.

### 1. Middle-Out Barcode and Primer Matching

Specimux first searches for the primer sequences, beginning at the specified search length away from each end 
of the sequence.  Once a primer is found, specimux searches the immediately adjacent portion of the sequence
for the appropriate barcodes.  This effectively anchors one edge of the barcode against the primer. This approach:
- Provides a more reliable starting point for barcode identification
- Reduces misidentification due to partial matches or errors at the barcode-primer junction
- Increases the effective distance between valid codewords, improving error correction capabilities

### 3. Enhanced Ambiguity Detection

Specimux aggressively detects multiple matches for different barcodes in the same sequence:
- Addresses a limitation in the original minibar.py program
- Introduces the `--ambiguity-threshold` parameter for user-controlled match stringency

### 4. Sequence Orientation Detection and Normalization

The program automatically detects and normalizes sequence orientation:
- Checks for primer matches to determine forward or reverse orientation
- Standardizes all sequences to forward orientation before barcode matching
- Ensures consistent processing regardless of original sequence orientation

### 5. Sequence Length Filtering

Specimux provides options to filter sequences based on length:
- `--min-length` parameter sets the minimum acceptable sequence length
- `--max-length` parameter sets the maximum acceptable sequence length
- Helps exclude potentially problematic sequences (e.g., truncated reads, artifacts, chimeras)

## Multiprocessing

Specimux utilizes multiprocessing to improve performance, particularly to offset the additional computational overhead introduced by the enhanced alignment checks. This parallel processing approach allows Specimux to leverage multi-core systems effectively.

### How it works

1. The input sequences are divided into batches.
2. Each batch is processed independently by a separate worker process.
3. The results from all workers are combined to produce the final output.

### Performance impact

On a typical 8-core system, the multiprocessing implementation in Specimux allows it to maintain performance comparable to the original minibar.py, despite the additional alignment checks. The actual performance gain may vary depending on your specific hardware and the nature of your input data.

### Controlling multiprocessing

By default, Specimux will use all available CPU cores. You can control the number of processes used with the `-t` or `--threads` option:
```
python specimux.py -t 4 barcode_file sequence_file
```
This example would limit Specimux to using 4 processes. If you want to disable multiprocessing entirely and run in a single process, you can set `-t 1`.

### Considerations

- Memory usage will increase with the number of processes, as each process needs its own memory space.
- If you're running on a shared system, be considerate of other users and avoid using all available cores.

By leveraging multiprocessing, Specimux aims to provide improved performance on multi-core systems while maintaining the accuracy and flexibility of its enhanced demultiplexing algorithm.

## Optimizing Edit Distance Parameters

The selection of appropriate edit distance parameters is crucial for balancing precision and recall in barcode assignment. Specimux allows users to fine-tune these parameters to optimize performance for their specific experimental setup.

### Importance of Parameter Selection

- **Too low edit distance**: May result in failure to recall and assign valid sequences, leading to data loss.
- **Too high edit distance**: Can reduce precision, potentially causing incorrect assignment of sequences to specimens.

Careful selection of a barcode set with high error-correcting power can simplify this optimization process. The paper by Buschmann and Bystrykh [1] provides valuable insights into the rationale and design of such barcode sets.

### Diagnostic Information

Specimux provides detailed information about edit distances when run in `--diagnostics` mode. This information can be crucial for understanding the error-correcting capabilities of your barcode set. For example:
```
2024-09-16 15:52:37,784 - INFO - Maximum edit distances set to 3 (barcode) and 6 (primer)
2024-09-16 15:52:37,785 - INFO - Minimum edit distance is 6 for Forward Barcodes
2024-09-16 15:52:37,796 - INFO - Minimum edit distance is 5 for Reverse Barcodes
2024-09-16 15:52:37,801 - INFO - Minimum edit distance is 4 for Forward Barcodes + Reverse Complement of Reverse Barcodes
2024-09-16 15:52:37,807 - INFO - Minimum edit distance is 4 for Forward Barcodes + Reverse Complement of Reverse Barcodes + Both Primers
```
This output helps users understand the robustness of their barcode set against various types of errors. By analyzing this information, researchers can make informed decisions about whether their chosen parameters strike the right balance between stringency and flexibility for their specific experimental needs.

We recommend users experiment with different parameter settings on a subset of their data to determine the optimal configuration for their particular sequencing setup and expected error rates.

## Sequence Output Options

By default, Specimux trims sequences to remove the forward and reverse primers, barcodes, and any other information between the 
primer and the respective end of the sequence.  If primers are not found, the sequence is not trimmed.  This behavior can be 
modified using the following options:

### No Trimming

You can disable all trimming with the `--trim none` option:
```
python specimux.py --trim none barcode_file sequence_file
```

This will output the full, unmodified sequences regardless of barcode and primer matches.

### Trimming tails only

You can leave the barcodes and primers intact, removing all other tails, with the `--trim tails` option:
```
python specimux.py --trim tails barcode_file sequence_file
```
Note that typically many sequences do not match one or more barcodes.  In this case, the sequence is
trimmed based on the expected primer and barcode lengths.

### Color Output

The `--color` option enhances the visual representation of barcode and primer matches:
```
python specimux.py --color barcode_file sequence_file
```
When using this option:
- Matched barcodes are highlighted in blue
- Matched primers are highlighted in green
- Bases with a quality score less than 10 are displayed in lowercase

Note: The color output uses ANSI escape codes and is intended for terminal viewing. Do not use this format for downstream processing as it may interfere with other tools.

## Diagnostic Output

Specimux offers two levels of diagnostic information:

### Basic Diagnostics

Use the `--diagnostics` option for additional summary information:
```
python specimux.py --diagnostics barcode_file sequence_file
```
This option provides:
- Detailed counts of different match types (e.g., full matches, partial matches, no matches)
- Statistics on barcode and primer edit distances
- Information about ambiguous matches and potential issues

### Debug Mode

For even more detailed information, use the `--debug` option:
```
python specimux.py --debug barcode_file sequence_file
```
In debug mode, Specimux will output:
- Detailed information for each processed sequence
- Step-by-step matching results for barcodes and primers
- Quality score information and its impact on matching (when applicable)
- Specific edit distances and locations for each match attempt

The debug output is verbose and is primarily intended for troubleshooting or understanding the internal workings of the algorithm.

## Barcode demultiplex file format

The program needs 5 pieces of information for each sample type. These are Sample ID, Forward Barcode index, Forward Primer, Reverse Barcode index, Reverse Primer. Even though the Forward Primer and Reverse Primer are the same for each sample in a run, this format requires them on every line describing a sample's indexes.

There can be a header line. However it and every sample line must have the same number of tab delimited fields.

Here's a simple example with the minimum of 5 tabbed delimited columns. 
```
SampleID	FwIndex	FwPrimer	RvIndex	RvPrimer
ONT01.01-A01-CM23-03848-Run25-iNat178308092	AGCAATCGCGCAC	CTTGGTCATTTAGAGGAAGTAA	AACCAGCGCCTAG	TCCTCCGCTTATTGATATGC
ONT01.02-B01-CM23-10261-Run25-iNat179083872	AGCAATCGCGCAC	CTTGGTCATTTAGAGGAAGTAA	ACTCGCGGTGCCA	TCCTCCGCTTATTGATATGC
ONT01.03-C01-CM23-10262-Run25-iNat179083983	AGCAATCGCGCAC	CTTGGTCATTTAGAGGAAGTAA	CCTGAATCATCTA	TCCTCCGCTTATTGATATGC
```
The header line is auto-detected by default.

## References

[1]: Buschmann T, Bystrykh LV. Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC Bioinformatics. 2013 Sep 11;14:272. doi: 10.1186/1471-2105-14-272. PMID: 24021088; PMCID: PMC3853030.

[2]: Henrik Krehenwinkel, Aaron Pomerantz, James B. Henderson, Susan R. Kennedy, Jun Ying Lim, Varun Swamy, Juan Diego Shoobridge, Nipam H. Patel, Rosemary G. Gillespie, Stefan Prost.  Nanopore sequencing of long ribosomal DNA amplicons enables portable and simple biodiversity assessments with high phylogenetic resolution across broad taxonomic scale. 
