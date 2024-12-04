# Specimux

Dual barcode and primer demultiplexing for MinION sequenced reads

Specimux is a substantial revision and enhancement of minibar.py, originally developed by the
California Academy of Sciences, Institute for Biodiversity Science & Sustainability.

Specimux is designed to improve the accuracy and throughput of DNA barcode identification for
multiplexed MinION sequencing data, with a primary focus on serving the fungal sequencing
community. Whereas minibar.py includes several processing methods supporting a variety of barcode
designs and matching regimes, specimux focuses specifically on high precision for demultiplexing
dual-indexed sequences.

Developed as a volunteer contribution, Specimux attempts to fill a small niche in processing large-scale,
community-driven fungal DNA barcoding projects.  Specimux was designed to work well with the Primary Data Analysis protocol
developed by Stephen Russell [1].

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

Specimux has been tested on MacOS and Linux machines.

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
                   [-n NUM_SEQS] [-e INDEX_EDIT_DISTANCE]
                   [-E PRIMER_EDIT_DISTANCE] [-l SEARCH_LEN]
                   [--group-unknowns] [-F] [-P OUTPUT_FILE_PREFIX]
                   [-O OUTPUT_DIR] [--color]
                   [--trim {none,tails,barcodes,primers}] [-d] [-D]
                   [--top-unmatched-barcodes TOP_UNMATCHED_BARCODES]
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
                        skipped (default: no filtering)
  --max-length MAX_LENGTH
                        Maximum sequence length. Longer sequences will be
                        skipped (default: no filtering)
  -n NUM_SEQS, --num-seqs NUM_SEQS
                        Number of sequences to read from file (e.g., -n 100 or
                        -n 102,3)
  -e INDEX_EDIT_DISTANCE, --index-edit-distance INDEX_EDIT_DISTANCE
                        Barcode edit distance value, default is half of min
                        distance between barcodes
  -E PRIMER_EDIT_DISTANCE, --primer-edit-distance PRIMER_EDIT_DISTANCE
                        Primer edit distance value, default is min distance
                        between primers
  -l SEARCH_LEN, --search-len SEARCH_LEN
                        Length to search for index and primer at start and end
                        of sequence (default: 80)
  --group-unknowns      Group unknown sequences based on partial matches and
                        classifications
  -F, --output-to-files
                        Create individual sample files for sequences
  -P OUTPUT_FILE_PREFIX, --output-file-prefix OUTPUT_FILE_PREFIX
                        Prefix for individual files when using -F (default:
                        sample_)
  -O OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory for individual files when using -F (default:
                        .)
  --color               Highlight barcode matches in blue, primer matches in
                        green
  --trim {none,tails,barcodes,primers}
                        trimming to apply
  -d, --diagnostics     Output extra diagnostics
  -D, --debug           Enable debug logging
  --top-unmatched-barcodes TOP_UNMATCHED_BARCODES
                        Display the top N unmatched barcode strings
  -t THREADS, --threads THREADS
                        Number of worker threads to use
  -v, --version         show program's version number and exit
```

## Major Features in Specimux

Specimux introduces several changes compared to the original minibar.py program, focused on
maximizing the precision of specimen assignment, minimizing wrongly-assigned specimens.

## Middle-Out Primer and Barcode Matching

Specimux uses a "middle-out" strategy to identify barcodes and primers in sequencing reads. Here's how it works:

1. Primer Detection:
   - The algorithm first searches for primer sequences within a specified search region at each end of the read
   - Primers are matched using "infix" alignment, allowing the match to float within the search region
   - The search region length can be adjusted using the --search-len parameter (default: 80bp)

2. Barcode Detection:
   - After finding a primer, Specimux looks for the corresponding barcode
   - For the forward (5') end:
     - Converts both the sequence and barcode to reverse complement
     - Searches in the region between the primer match and the sequence start
     - Must align immediately adjacent to the primer (prefix mode)
   - For the reverse (3') end:
     - Uses forward orientation
     - Searches in the region between the primer match and the sequence end
     - Must align immediately adjacent to the primer (prefix mode)
   - Barcode search example:
```
5'                                                                                    3'
[Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]
                 ^                                                    ^
      <-- Search |                                                    | Search -->
```          

3. Edit Distance Scoring:
   - Any gaps between primer and barcode count against the allowed edit distance
   - This anchoring effect helps prevent false matches from partial sequence similarities
   - Barcode matches must stay within the specified edit distance threshold
   - Default threshold is half the minimum edit distance between any barcode pair

Benefits of this approach:
- More reliable identification by using primers as anchor points
- Reduced false positives from partial matches or chimeric sequences 
- Better error correction by enforcing barcode-primer proximity
- Handles both perfect and imperfect matches in a principled way

## Sequence Orientation Detection and Normalization

The program automatically detects and normalizes sequence orientation:

- Checks for primer matches to determine forward or reverse orientation
- Standardizes all sequences to forward orientation before barcode matching
- Ensures consistent processing regardless of original sequence orientation

## Sequence Length Filtering

Specimux provides options to filter sequences based on length:

- `--min-length` parameter sets the minimum acceptable sequence length
- `--max-length` parameter sets the maximum acceptable sequence length
- Helps exclude potentially problematic sequences (e.g., truncated reads, artifacts, chimeras)

## Multiprocessing

Specimux utilizes multiprocessing to improve performance, particularly to offset the additional computational overhead
introduced by the enhanced alignment checks. This parallel processing approach allows Specimux to leverage multi-core
systems effectively.

### How it works

1. The input sequences are divided into batches by the parent process.
2. Each batch is processed independently by a separate worker process.
3. The results from all workers are combined in the parent process to produce the final output.

### Performance impact

On a typical 8-core system, the multiprocessing implementation in Specimux allows it to maintain performance comparable
to the original minibar.py, despite the additional alignment checks. The actual performance gain may vary depending on
your specific hardware and the nature of your input data.

#### Multiple Primer Pairs and Degenerate Primers
Specimux supports multiple primer pairs and degenerate primers containing ambiguous IUPAC codes.  However, these do
not come for free.  Each additional pair of primers will slow performance.  Degenerate primers may slow performance 
significantly.

### Controlling multiprocessing

By default, Specimux will use all available CPU cores. You can control the number of processes used with the `-t` or
`--threads` option:

```
python specimux.py -t 4 barcode_file sequence_file
```

This example would limit Specimux to using 4 processes. If you want to disable multiprocessing entirely and run in a
single process, you can set `-t 1`.

### Considerations

- Memory usage will increase with the number of processes, as each process needs its own memory space.
- If you're running on a shared system, be considerate of other users and avoid using all available cores.

By leveraging multiprocessing, Specimux aims to provide improved performance on multi-core systems while maintaining the
accuracy and flexibility of its enhanced demultiplexing algorithm.

## Optimizing Edit Distance Parameters

The selection of appropriate edit distance parameters is crucial for balancing precision and recall in barcode
assignment. Specimux allows users to fine-tune these parameters to optimize performance for their specific experimental
setup.

### Importance of Parameter Selection

- **Too low edit distance**: May result in failure to recall and assign valid sequences, leading to data loss.
- **Too high edit distance**: Can reduce precision, potentially causing incorrect assignment of sequences to specimens.

Careful selection of a barcode set with high error-correcting power can simplify this optimization process. The paper by
Buschmann and Bystrykh [2] provides valuable insights into the rationale and design of such barcode sets.

By default, specimux will assign the barcode edit distance threshold as half (round up) of the minimum edit distance between all
barcodes specified in the metadata file.  This can be overridden with the `-e` parameter.

By default, specimux will assign the primer edit distance threshold as the minimum edit distance between all primers specified
in the metadata file.  This can be overridden with the `-E` parameter.

## Multiple Primer Support

Specimux now supports demultiplexing of mixed-primer runs, where different specimens may use different primer pairs.
This is particularly useful for fungal DNA barcoding projects where some specimens may be amplified with ITS1+ITS2
primers while others use ITS2-only primers or primers for other loci.

For each specimen in the barcode file, you specify both the primers and barcodes used. Specimux will:

- Track which specimens use which primer pairs
- Attempt to match each sequence against all primer pairs in use
- Choose the best match based on a scoring hierarchy that considers both primers and barcodes
- Report any ambiguous cases where multiple primer pairs match equally well

When multiple primer pairs are in use, the orientation detection and matching logic becomes more complex. Specimux
handles this by:

1. First determining sequence orientation using all available primer pairs
2. Then attempting full matching with each primer pair
3. Selecting the best match based on a clear scoring hierarchy:
    - Full matches (both primers and both barcodes) score highest
    - Partial matches are scored progressively lower
    - When matches of equal quality occur with different primer pairs, the sequence is marked as ambiguous

## Understanding Classification Statistics

When processing sequences, Specimux follows a specific decision tree that determines how each sequence is classified.
Understanding these classifications can help diagnose issues with your data or parameters.

The matching process follows this hierarchy:

1. **Length Checks** (performed first)
    - "Sequence Too Short" - Sequence is shorter than --min-length
    - "Sequence Too Long" - Sequence is longer than --max-length

2. **Orientation Detection**
    - "Could Not Determine Orientation" - No clear primer matches to establish sequence direction

3. **Primer Matching**
    - "No Primer Matches" - Neither forward nor reverse primer found
    - "No Forward Primer Matches" - Forward primer not found but reverse primer present
    - "No Reverse Primer Matches" - Reverse primer not found but forward primer present

4. **Barcode Matching** (only attempted if primers found)
    - "No Forward Barcode Matches" - No matching forward barcode after forward primer
    - "No Reverse Barcode Matches" - No matching reverse barcode after reverse primer
    - "No Forward Barcode Matches (May be truncated)" - No forward barcode found but sequence is short at that end
    - "No Reverse Barcode Matches (May be truncated)" - No reverse barcode found but sequence is short at that end

5. **Multiple Primer Ambiguity**
   - "Multiple Primer Pair Full Matches" - Multiple primer pairs matched with both barcodes also matching
6. **Specimen Ambiguity Checks**
    - "Multiple Matches for Forward Barcode" - Multiple forward barcodes match within the ambiguity threshold
    - "Multiple Matches for Reverse Barcode" - Multiple reverse barcodes match within the ambiguity threshold
    - "Multiple Matches for Both Barcodes" - Both ends have multiple matching barcodes
    - "Multiple Primer Pair Matches" - Multiple primer pairs match equally well (only with mixed primer sets)
7. **Success**
    - "Matched" - Successful unambiguous match of primers and barcodes

Each classification represents a terminal state - once a sequence receives a classification, no further matching is
attempted. This helps both with performance (by avoiding unnecessary matching attempts) and with diagnostic clarity (by
identifying the specific point where matching failed).

The statistics output shows the count and percentage for each classification, ordered by frequency. This can help
identify systematic issues - for example, a high rate of orientation failures might suggest primer sequence problems,
while many truncated matches might indicate sequencing quality issues or upstream data pipeline issues.

### Diagnostic Information

Specimux provides detailed information about edit distances when run in `-d` / `--diagnostics` mode. This information can be
important for understanding the error-correcting capabilities of your barcode set. For example:

```
2024-09-16 15:52:37,784 - INFO - Maximum edit distances set to 3 (barcode) and 6 (primer)
2024-09-16 15:52:37,785 - INFO - Minimum edit distance is 6 for Forward Barcodes
2024-09-16 15:52:37,796 - INFO - Minimum edit distance is 5 for Reverse Barcodes
2024-09-16 15:52:37,801 - INFO - Minimum edit distance is 4 for Forward Barcodes + Reverse Complement of Reverse Barcodes
2024-09-16 15:52:37,807 - INFO - Minimum edit distance is 4 for Forward Barcodes + Reverse Complement of Reverse Barcodes + Both Primers
```

This output helps users understand the robustness of their barcode set against various types of errors. By analyzing
this information, users can make informed decisions about whether their chosen parameters strike the right balance
between stringency and flexibility for their specific experimental needs.

We recommend users experiment with different parameter settings on a subset of their data to determine the optimal
configuration for their particular sequencing setup and expected error rates.

## Sequence Output Options

Each sequence can have several distinct regions, trimmed differently by each option:

Raw sequence layout:
```
5'                                                                                                            3'
<---tail--->[Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]<---tail--->
```

--trim none:
```
<---tail--->[Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]<---tail--->
|<-----------------------------------------entire sequence unchanged----------------------------------------->|
```
--trim tails:
```
            [Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]
            |<----------------------------------output sequence---------------------------------->|
```
--trim barcodes (default):
```
                           [Forward Primer]---target sequence---[Reverse Primer]
                           |<-----------------output sequence----------------->|
```
--trim primers:
```
                                         ---target sequence---
                                         |<-output sequence->|
```

Note: When barcodes cannot be located, their positions are estimated based on expected barcode length
when using ```--trim tails```. Trimming of primers only occurs when primers are explicitly found.

By default, Specimux trims sequences to remove the forward and reverse barcodes, and any other information
between the primer and the respective end of the sequence. Primers are left intact.  If primers are not found, 
the sequence is not trimmed. This behavior can be modified using the following options.

### No Trimming

You can disable all trimming with the `--trim none` option:

```
python specimux.py --trim none barcode_file sequence_file
```

This will output the full, unmodified sequences regardless of barcode and primer matches.  This trimming mode is especially
useful for debugging, in combination with the `--color` output.

### Trimming tails only

You can leave the barcodes and primers intact, removing all other tails, with the `--trim tails` option:

```
python specimux.py --trim tails barcode_file sequence_file
```

Note that typically many sequences do not match one or more barcodes. In this case, the sequence is
trimmed based on the expected primer and barcode lengths.

### Full Trimming, including primers

You can remove everything except the target sequence with the `--trim primers` option:

```
python specimux.py --trim primers barcode_file sequence_file
```

Note, however, that downstream processing such as consensus generation may benefit from leaving the primers
intact as anchor points for alignment.

### Color Output

The `--color` option enhances the visual representation of barcode and primer matches:

```
python specimux.py --color barcode_file sequence_file
```

When using this option:

- Matched barcodes are highlighted in blue
- Matched primers are highlighted in green
- Bases with a quality score less than 10 are displayed in lowercase

Note: The color output uses ANSI escape codes and is intended for terminal viewing. Do not use this format for
downstream processing as it may interfere with other tools.

## Output File Organization

When using `-F` / `--output-to-files`:
- Each specimen gets its own output file named {prefix}{specimen_id}.{extension}
- Invalid characters in specimen IDs are replaced with underscores
- Files are created in the specified output directory (default: current directory)
- Unknown sequences go to {prefix}unknown.{extension}
- Ambiguous matches go to {prefix}ambiguous.{extension}

### Experimental Grouping of Unknown Sequences

The `--group-unknowns` option provides an experimental feature for debugging matching problems by grouping sequences 
that share partial matches. Instead of putting all unmatched sequences into a single "unknown" file, it creates separate 
files based on partial match patterns:

- {prefix}barcode_fwd_{barcode}.{extension} - Sequences where only the forward barcode matched
- {prefix}barcode_rev_{barcode}.{extension} - Sequences where only the reverse barcode matched

This grouping can help identify:
- Truncated sequences where one end is missing
- Barcode synthesis or primer binding issues
- Potential contamination or cross-talk between specimens
- Novel barcodes that might indicate sample tracking errors

Note that this is a diagnostic tool and should not be used for production analysis. The groupings are based on partial 
matches that did not meet the full matching criteria and may contain incorrectly assigned sequences.  This feature may
be removed or enhanced in future versions.  Feedback is welcome.

## Diagnostic Output

Specimux offers two levels of diagnostic information:

### Basic Diagnostics

Use the `-d` or `--diagnostics` option for additional summary information:

```
python specimux.py --diagnostics barcode_file sequence_file
```

This option provides:

- Detailed counts of different match types (e.g., full matches, partial matches, no matches)
- Statistics on barcode and primer edit distances
- Information about ambiguous matches and potential issues

### Debug Mode

For even more detailed information, use the `-D` or `--debug` option:

```
python specimux.py --debug barcode_file sequence_file
```

In debug mode, Specimux will output:

- Detailed information for each processed sequence
- Step-by-step matching results for barcodes and primers
- Quality score information and its impact on matching (when applicable)
- Specific edit distances and locations for each match attempt

The debug output is verbose and is primarily intended for troubleshooting or understanding the internal workings of the
algorithm.

## Barcode demultiplex file format

The program needs 5 pieces of information for each sample type. These are Sample ID, Forward Barcode index, Forward
Primer, Reverse Barcode index, Reverse Primer. Even though the Forward Primer and Reverse Primer are often the same for each
sample in a run, this format requires them on every line describing a sample's indexes.

There can be a header line. However it and every sample line must have the same number of tab delimited fields.

Here's a simple example with the minimum of 5 tabbed delimited columns.

```
SampleID	FwIndex	FwPrimer	RvIndex	RvPrimer
ONT01.01-A01-CM23-03848-Run25-iNat178308092	AGCAATCGCGCAC	CTTGGTCATTTAGAGGAAGTAA	AACCAGCGCCTAG	TCCTCCGCTTATTGATATGC
ONT01.02-B01-CM23-10261-Run25-iNat179083872	AGCAATCGCGCAC	CTTGGTCATTTAGAGGAAGTAA	ACTCGCGGTGCCA	TCCTCCGCTTATTGATATGC
ONT01.03-C01-CM23-10262-Run25-iNat179083983	AGCAATCGCGCAC	CTTGGTCATTTAGAGGAAGTAA	CCTGAATCATCTA	TCCTCCGCTTATTGATATGC
```

The header line is auto-detected (and ignored) by default.

## References
[1]: Stephen Douglas Russell 2023. Primary Data Analysis - Basecalling, Demultiplexing, and Consensus Building for ONT Fungal Barcodes. 
protocols.io https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v3

[2]: Buschmann T, Bystrykh LV. Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC Bioinformatics.
2013 Sep 11;14:272. doi: 10.1186/1471-2105-14-272. PMID: 24021088; PMCID: PMC3853030.

[3]: Henrik Krehenwinkel, Aaron Pomerantz, James B. Henderson, Susan R. Kennedy, Jun Ying Lim, Varun Swamy, Juan Diego
Shoobridge, Nipam H. Patel, Rosemary G. Gillespie, Stefan Prost. Nanopore sequencing of long ribosomal DNA amplicons
enables portable and simple biodiversity assessments with high phylogenetic resolution across broad taxonomic scale. 
