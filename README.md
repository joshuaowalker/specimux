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

Primers are specified in FASTA format with metadata in the description line:

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

Specimux organizes output by pool and primer pair:

```
output_dir/
  ITS/
    ITS1F-ITS4/
      sample_specimen1.fastq
    unknown/
      sample_unknown.fastq
  RPB2/
    fRPB2-5F-RPB2-7.1R/
      sample_specimen2.fastq
    unknown/
      sample_unknown.fastq
  unknown/
    sample_unknown.fastq
```

Each pool directory contains:
- Subdirectories for each primer pair combination
- "unknown" directory for sequences matching primers but not barcodes
- Separate files for ambiguous matches

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
   |<-----------------------------------------------full sequence----------------------------------------------->|

--trim tails: (remove external regions)
               |<----------------------------------trimmed sequence--------------------------------->|

--trim barcodes: (default, remove barcodes)
                                |<-----------------trimmed sequence---------------->|

--trim primers: (remove primers)
                                                |<trimmed  sequence>|
```

### Multiprocessing

- Enabled with -F/--output-to-files option
- Uses all cores by default, controllable with --threads
- Improves performance through parallel processing
- Memory usage increases with thread count

### Performance Optimizations

#### Bloom Filter Prefiltering (v0.4)
- Quick elimination of impossible barcode matches
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

## Version History
- 0.5 (February 2025): Added Primer Pools and Hierarchical Output
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