# Specimux

Dual barcode and primer demultiplexing for MinION sequenced reads

Specimux is a substantial revision and enhancement of minibar.py, originally developed by the
California Academy of Sciences, Institute for Biodiversity Science & Sustainability.

Specimux is designed to improve the accuracy and throughput of DNA barcode identification for
multiplexed MinION sequencing data, with a primary focus on serving the fungal sequencing
community. Whereas minibar.py includes several processing methods supporting a variety of barcode
designs and matching regimes, specimux focuses specifically on high precision for demultiplexing
dual-indexed sequences.

The tool was developed and tested using the Mycomap ONT037 dataset, which comprises 768 specimens
and approximately 765,000 nanopore reads in FastQ format. This real-world dataset provided a robust
testing ground, ensuring Specimux's capabilities align closely with the needs of contemporary
fungal biodiversity research. Specimux was designed to work seamlessly with the Primary Data Analysis protocol developed by Stephen Russell [1], serving the needs of community-driven fungal DNA barcoding projects.

## Requirements

**specimux.py** is written in Python and is compatible with Python 3. It requires several Python packages which can be installed using the provided `requirements.txt` file:

```bash 
pip install -r requirements.txt
```

Specimux has been tested on MacOS and Linux machines.

## Basic Usage

Specimux uses primer pools to organize specimens and their associated primers. Here's a basic example:

1. Define primers and their pools (primers.fasta):
```
>ITS1F pool=ITS position=forward
CTTGGTCATTTAGAGGAAGTAA
>ITS4 pool=ITS position=reverse
TCCTCCGCTTATTGATATGC
```

2. Create specimen file mapping barcodes to pools (specimens.txt):
```
SampleID    PrimerPool    FwIndex    FwPrimer    RvIndex    RvPrimer
specimen1   ITS           ACGTACGT   ITS1F       TGCATGCA   ITS4
specimen2   ITS           GTACGTAC   ITS1F       CATGCATG   ITS4
```

3. Run specimux:
```bash
python specimux.py primers.fasta specimens.txt sequences.fastq -F -d
```

For a full list of options:
```bash
python specimux.py -h
```

## Primer Pool Organization

Primer pools are a core organizing principle in Specimux, allowing logical grouping of primers and specimens. A pool defines:
- Which primers can be used together
- Which specimens belong to which primer sets
- How output files are organized

### Pool Design Benefits
- Organize specimens by target region (e.g., ITS, RPB2)
- Support shared primers between pools
- Improve performance by limiting primer search space
- Provide logical output organization

### Primer File Format

Primers are specified in a text file in FASTA format with metadata in the description line:

```
>primer_name pool=pool1,pool2 position=forward
PRIMER_SEQUENCE
```

Required metadata:
- `pool=` - Comma/semicolon separated list of pool names
- `position=` - Either "forward" or "reverse"

Example for fungal ITS and RPB2 regions:
```
>ITS1F pool=ITS,Mixed position=forward
CTTGGTCATTTAGAGGAAGTAA
>ITS4 pool=ITS position=reverse
TCCTCCGCTTATTGATATGC
>fRPB2-5F pool=RPB2 position=forward
GAYGAYMGWGATCAYTTYGG
>RPB2-7.1R pool=RPB2 position=reverse
CCCATRGCYTGYTTMCCCATDGC
```

Although the file is technically in FASTA format, you can name it primers.fasta, primers.txt, or anything that makes sense for your workflow.

### Specimen File Format 

Tab-separated file with columns:
- SampleID - Unique identifier for specimen
- PrimerPool - Which pool the specimen belongs to
- FwIndex - Forward barcode sequence
- FwPrimer - Forward primer name or wildcard (*/-) 
- RvIndex - Reverse barcode sequence
- RvPrimer - Reverse primer name or wildcard (*/-) 

Example:
```
SampleID         PrimerPool  FwIndex         FwPrimer  RvIndex         RvPrimer
specimen1        ITS         ACGTACGT        ITS1F     TGCATGCA        ITS4
specimen2        RPB2        GTACGTAC        *         CATGCATG        *
```

### Output Organization

Specimux organizes output with match quality at the top level, making it easy to access your primary data (full matches) while keeping partial matches and unknowns organized separately:

```
output_dir/
  full/                            # All complete matches (PRIMARY DATA)
    ITS/                           # Pool-level aggregation
      sample_specimen1.fastq       # All ITS full matches collected here
      sample_specimen2.fastq
      primers.fasta                # All primers in the ITS pool
      ITS1F-ITS4/                  # Primer-pair specific matches
        sample_specimen1.fastq
        sample_specimen2.fastq
        primers.fasta              # Just this primer pair
    RPB2/
      sample_specimen3.fastq
      primers.fasta
      fRPB2-5F-RPB2-7.1R/
        sample_specimen3.fastq
        primers.fasta
  
  partial/                         # One barcode matched (RECOVERY CANDIDATES)
    ITS/
      ITS1F-ITS4/
        sample_barcode_fwd_ACGTACGT.fastq
      ITS1F-unknown/               # Forward primer only detected
        sample_barcode_fwd_ACGTACGT.fastq
      unknown-ITS4/                # Reverse primer only detected
        sample_barcode_rev_TGCATGCA.fastq
  
  ambiguous/                       # Multiple specimens matched
    ITS/
      ITS1F-ITS4/
        sample_ambiguous.fastq
  
  unknown/                         # No barcodes matched
    ITS/
      ITS1F-ITS4/                  # Primers detected but no barcodes
        sample_unknown.fastq
    unknown/
      unknown-unknown/             # No primers detected at all
        sample_unknown.fastq
  
  log.txt                          # Run statistics and classification report
```

#### Directory Structure Benefits

The match-type-first organization provides several advantages:

1. **Primary Data Access**: Full matches are immediately accessible in the `full/` directory without navigating through multiple subdirectories
2. **Clean Separation**: Partial matches and unknowns are segregated, reducing clutter when accessing your primary demultiplexed data
3. **Convenient Aggregation**: Pool-level directories (e.g., `full/ITS/`) collect all successful matches for that target region
4. **Recovery Options**: The `partial/` directory contains sequences that may be recoverable using tools like `speciharvest.py`
5. **Automatic Cleanup**: Empty directories are automatically removed after processing to keep the output clean

#### Key Directories

- **full/[pool]/**: Contains ALL full matches for that pool, regardless of primer pair. Sequences appear both here and in their specific primer-pair subdirectory for maximum flexibility
- **full/[pool]/[primer1-primer2]/**: Contains full matches for this specific primer pair only
- **partial/[pool]/[primer1-unknown]/**: Contains sequences where only one primer was detected (potential recovery candidates)
- **unknown/unknown/unknown-unknown/**: Contains sequences where no primers could be identified

### Pool Management

Specimux automatically:
- Validates pool configurations (minimum one forward/reverse primer)
- Tracks which pools are used by specimens
- Prunes unused pools
- Reports pool statistics

Example output:
```
INFO - Loaded 4 primers in 3 pools
INFO - Pool RPB2: 1 forward, 1 reverse primers
INFO - Pool TEF1: 1 forward, 1 reverse primers
INFO - Removing unused pools: MultiLocus
```

## Sequence Matching Strategy

Specimux uses a "middle-out" strategy to identify primers and barcodes:

1. Primer Detection:
   - Search within specified region at each end (--search-len, default: 80bp)
   - Use "infix" alignment allowing float within search region
   - Match against primers from assigned pool

2. Barcode Detection:
   - After finding primer, look for corresponding barcode
   - Must align immediately adjacent to primer
   - Forward (5') end: search between primer and sequence start
   - Reverse (3') end: search between primer and sequence end

```
5'                                                                                    3'
[Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]
                 ^                                                    ^
      <-- Search |                                                    | Search -->
```          

3. Match Scoring:
   - Full matches (both primers + barcodes) score highest
   - Partial matches scored progressively lower
   - Pool consistency considered in scoring
   - Ambiguity reported when multiple high-scoring matches exist

All sequences are automatically normalized to forward orientation after matching, ensuring consistent output regardless of input orientation.

## Edit Distance Parameters

The selection of appropriate edit distance parameters is crucial for balancing precision and recall in sequence assignment. Specimux provides separate controls for barcode and primer matching:

### Barcode Edit Distance
- Default is half (rounded up) of the minimum edit distance between any barcode pair
- Provides good balance between error tolerance and uniqueness
- Can override with -e/--index-edit-distance parameter
- Too low: Fails to recall valid sequences with errors
- Too high: May cause incorrect assignments

### Primer Edit Distance
- Default based on primer complexity and IUPAC codes
- Accounts for degenerate bases in calculation
- Can override with -E/--primer-edit-distance parameter
- More tolerant than barcode matching due to primer length

### Optimization Tips
- Use high-quality barcode sets with good error-correcting properties
- Consider the paper by Buschmann and Bystrykh [2] for barcode design
- Test on subset of data to find optimal parameters
- Use -d/--diagnostics to view edit distance statistics

## Sequence Processing Options

### Trimming Modes

Each mode trims a different portion of the sequence:

```
Raw sequence:
5' <---tail--->[Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]<---tail---> 3'

--trim none: (entire sequence unchanged)
5' <---tail--->[Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]<---tail---> 3'

--trim tails: (remove external regions)
5'             [Forward Barcode][Forward Primer]---target sequence---[Reverse Primer][Reverse Barcode]             3'

--trim barcodes: (default, remove barcodes)
5'                              [Forward Primer]---target sequence---[Reverse Primer]                              3'

--trim primers: (remove primers)
5'                                              ---target sequence---                                              3'
```

### Multiprocessing

- Enabled with -F/--output-to-files option
- Uses all cores by default, controllable with --threads
- Improves performance through parallel processing
- Memory usage increases with thread count

### Performance Optimizations

#### Bloom Filter Prefiltering (v0.4)
- Uses hashing before sequence alignment to speed up barcode matching
- Best for barcodes ≤13nt and edit distances ≤3
- Can disable with --disable-prefilter

#### Sequence Pre-orientation
- Heuristic orientation detection
- Reduces alignment operations
- Can disable with --disable-preorient

#### Pool-Based Optimization
- Limits primer search space to active pools
- Automatic pruning of unused pools
- Efficient file organization and buffering

## Diagnostic Features

### Run Log

At the end of each run, Specimux creates a log.txt file in the output directory containing:
- Date and time of the run
- Command line used
- Input files (primer, specimen, and sequence files)
- Run statistics (total sequences, match rate, processing time, etc.)
- Detailed classification statistics

This log is generated regardless of whether the -d/--diagnostics flag is used and provides a permanent record of each processing run.

### Classification Statistics (-d)

Shows sequence assignment outcomes:

1. Length Checks
   - "Sequence Too Short/Long"

2. Primer Matching
   - "No Primer Matches"
   - "No Forward/Reverse Primer Matches"

3. Barcode Matching
   - "No Forward/Reverse Barcode Matches"
   - "No Barcode Matches (May be truncated)"

4. Ambiguity Checks
   - "Multiple Primer/Orientation Full Matches"
   - "Multiple Matches for Forward/Reverse Barcode"
   - "No Specimen for Barcodes"

5. Success
   - "Matched"

### Debug Output (-D)

Provides detailed matching information:
- Step-by-step alignment results
- Quality score impacts
- Edit distances and locations
- Pool assignment decisions

### Quality Visualization (--color)

Highlights sequence components:
- Barcodes in blue
- Primers in green
- Low quality bases (<Q10) in lowercase

## Processing Flow Visualization

Specimux generates a `stats.json` file containing detailed statistics about sequence processing flow, compatible with Sankey diagram visualization tools. This provides a visual representation of how sequences move through the processing pipeline: primer detection → outcome classification → pool assignment.

### Visualization Tool

The included `visualize_stats.py` script creates interactive Sankey diagrams from the stats.json file:

```bash
# Basic usage
python visualize_stats.py stats.json

# Custom output file
python visualize_stats.py stats.json my_flow_diagram.html

# Custom dimensions
python visualize_stats.py stats.json --width 1600 --height 800
```

### Dependencies

The visualization tool requires plotly:

```bash
pip install plotly
```

### Features

- **Interactive Diagrams**: Hover over nodes and flows to see exact counts
- **Color Coding**: Different colors for primer types, outcomes, and pools
- **Customizable**: Adjustable dimensions for different display needs
- **Self-Contained**: Generated HTML files work offline and can be shared easily

### Flow Stages

The diagram shows sequence processing through these stages:

1. **Total Sequences**: Starting point showing all input sequences
2. **First Primer Detection**: Sequences grouped by detected forward primer
3. **Primer Pair Formation**: Complete primer pair identification
4. **Outcome Classification**: Matched, Partial, Ambiguous, or Unknown outcomes
5. **Pool Assignment**: Final assignment to primer pools

This visualization helps identify:
- Which primers are most/least effective
- Where sequences are lost in the pipeline
- Pool assignment patterns and potential issues
- Overall processing efficiency

## Specimen File Converter

If you're upgrading from an earlier version of specimux (or from minibar.py), a converter tool is included to help migrate your specimen files to the new format.

### Legacy Format

Earlier versions used a different specimen file format that included primer sequences directly:

```
SampleID         FwIndex         FwPrimer                    RvIndex         RvPrimer
ONT01.01-A01     AGCAATCGCGCAC   CTTGGTCATTTAGAGGAAGTAA      AACCAGCGCCTAG   TCCTCCGCTTATTGATATGC
ONT01.02-B01     AGCAATCGCGCAC   CTTGGTCATTTAGAGGAAGTAA      ACTCGCGGTGCCA   TCCTCCGCTTATTGATATGC
```

### Converter Tool

The `v04_specimen_converter.py` script automatically:

1. Extracts all unique primer sequences
2. Generates a `primers.fasta` file with proper pool annotations
3. Creates a new specimen file with the required `PrimerPool` column
4. Replaces primer sequences with primer names

### Usage

```bash
python v04_specimen_converter.py Index.txt --output-specimen=IndexPP.txt --output-primers=primers.fasta --pool-name=ITS
```

### Arguments

- `input_file`: The old format specimen file (required)
- `--output-specimen`: Path for the new format specimen file (default: specimen_new.txt)
- `--output-primers`: Path for the primers FASTA file (default: primers.fasta)
- `--pool-name`: Name to use for the primer pool (default: pool1)

## Version History
- 0.5.1 (August 2025): Reorganized output with match-type-first directory structure for easier access to primary data, added automatic cleanup of empty directories, added processing flow statistics (stats.json) and interactive Sankey diagram visualization
- 0.5 (March 2025): Added Primer Pools, Hierarchical Output with pool-level full match collections, and detailed run log
- 0.4 (February 2025): Added Bloom filter optimization
- 0.3 (December 2024): Code cleanup and write pooling improvements
- 0.2 (November 2024): Multiple primer pair support
- 0.1 (September 2024): Initial release

## Troubleshooting

### File Organization
Specimux maintains several directory structures for efficient operation:

~/.specimux/cache/
- Stores cached Bloom filters for barcode matching
- Safe to delete if issues arise
- Will be recreated as needed

output_dir/.specimux_locks/
- Temporary lock files for multiprocessing
- Automatically cleaned up after successful runs
- Can be safely deleted if program crashes

### Common Issues
- If multiprocessing hangs: Remove .specimux_locks directory
- If matching seems incorrect: Clear ~/.specimux/cache
- If file output fails: Check directory permissions

## References

[1]: Stephen Douglas Russell 2023. Primary Data Analysis - Basecalling, Demultiplexing, and Consensus Building for ONT Fungal Barcodes. 
protocols.io https://dx.doi.org/10.17504/protocols.io.dm6gpbm88lzp/v3

[2]: Buschmann T, Bystrykh LV. Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC Bioinformatics.
2013 Sep 11;14:272. doi: 10.1186/1471-2105-14-272.