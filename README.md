# Specimux

Dual barcode and primer demultiplexing for MinION sequenced reads

Specimux is an independent project inspired by minibar.py (originally developed by the
California Academy of Sciences). While building upon core demultiplexing concepts from minibar,
Specimux represents a complete reimplementation with substantial algorithmic enhancements and
architectural improvements.

Specimux is designed to improve the accuracy and throughput of DNA barcode identification for
multiplexed MinION sequencing data, with a primary focus on serving the fungal sequencing
community. Whereas minibar.py includes several processing methods supporting a variety of barcode
designs and matching regimes, specimux focuses specifically on high precision for demultiplexing
dual-indexed sequences.

The tool was developed and tested using the Mycomap ONT037 dataset, which comprises 768 specimens
and approximately 765,000 nanopore reads in FastQ format. This real-world dataset provided a robust
testing ground, ensuring Specimux's capabilities align closely with the needs of contemporary
fungal biodiversity research. Specimux was designed to work seamlessly with the Primary Data Analysis protocol developed by Stephen Russell [1], serving the needs of community-driven fungal DNA barcoding projects.

## Installation

### Option 1: Install from GitHub (Recommended)

**Virtual Environment Recommended**: It's strongly recommended to use a virtual environment to avoid dependency conflicts:

```bash
# Create and activate virtual environment
python3 -m venv specimux-env
source specimux-env/bin/activate  # On Windows: specimux-env\Scripts\activate

# Install latest version (includes visualization support)
pip install git+https://github.com/joshuaowalker/specimux.git

# Install with development tools
pip install "git+https://github.com/joshuaowalker/specimux.git#egg=specimux[dev]"
```

After installation, specimux commands are available:
```bash
specimux --version
specimux primers.fasta specimens.txt sequences.fastq -F -d
```

**Note**: Remember to activate your virtual environment (`source specimux-env/bin/activate`) each time you want to use specimux.

### Option 2: Local Development Installation

For development or testing modifications:

```bash
# Clone the repository
git clone https://github.com/joshuaowalker/specimux.git
cd specimux

# Create virtual environment (Python 3.10+ required)
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e .

# Install with development tools
pip install -e ".[dev]"
```

### Requirements

**Python Version**: Specimux requires Python 3.10 or newer, with full support for Python 3.10-3.13.

Specimux automatically installs these dependencies:
- edlib>=1.1.2 (sequence alignment)
- biopython>=1.81 (sequence handling)
- pybloomfilter3>=0.7.3 (performance optimization)
- cachetools>=5.3.0 (file handle caching)
- tqdm>=4.65.0 (progress bars)
- plotly>=5.0.0 (visualization support)
- watchdog>=3.0.0 (file system monitoring for specimux-watch)

Specimux has been tested on MacOS and Linux machines.

## Available Commands

After installation, specimux provides several command-line tools:

- **`specimux`** - Main demultiplexer for dual barcode and primer matching
- **`specimux-watch`** - Automatic file watcher for live MinKNOW sequencing workflows
- **`specimine`** - Mine additional sequences from partial barcode matches
- **`specimux-convert`** - Convert legacy specimen files to current format
- **`specimux-stats`** - Analyze trace files to generate statistics
- **`specimux-visualize`** - Create interactive Sankey diagrams from statistics

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
specimux primers.fasta specimens.txt sequences.fastq -F -d
```

### Dereplication

When the same read matches multiple primer pairs (common with complex primer pools), specimux deduplicates the output by default:

```bash
# Default: Select best match per specimen/barcode group (recommended)
specimux primers.fasta specimens.txt sequences.fastq --dereplicate best

# Output all equivalent matches (may cause read amplification)
specimux primers.fasta specimens.txt sequences.fastq --dereplicate none
```

The `best` strategy selects the optimal match per specimen using tiebreakers: barcode distance, primer count, primer distance, and file order. This prevents artificial read amplification that can occur when overlapping primers cause the same read to match through multiple permutations.

Note: If a read legitimately matches multiple different specimens (e.g., through different barcode combinations), it will still appear in each specimen's output file. Dereplication only prevents duplicate output when the same read matches the same specimen through different primer pair routes.

For a full list of options:
```bash
specimux -h
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
      specimen1.fastq              # All ITS full matches collected here
      specimen2.fastq
      primers.fasta                # All primers in the ITS pool
      ITS1F-ITS4/                  # Primer-pair specific matches
        specimen1.fastq
        specimen2.fastq
        primers.fasta              # Just this primer pair
    RPB2/
      specimen3.fastq
      primers.fasta
      fRPB2-5F-RPB2-7.1R/
        specimen3.fastq
        primers.fasta
  
  partial/                         # One barcode matched (RECOVERY CANDIDATES)
    ITS/
      ITS1F-ITS4/
        barcode_fwd_ACGTACGT.fastq
      ITS1F-unknown/               # Forward primer only detected
        barcode_fwd_ACGTACGT.fastq
      unknown-ITS4/                # Reverse primer only detected
        barcode_rev_TGCATGCA.fastq
  
  # Note: ambiguous/ directory removed in v0.5+
  # Multiple equivalent matches now output to their respective specimen files
  
  unknown/                         # No barcodes matched
    ITS/
      ITS1F-ITS4/                  # Primers detected but no barcodes
        unknown.fastq
    unknown/
      unknown-unknown/             # No primers detected at all
        unknown.fastq
  
  trace/                           # Diagnostic trace files (with -d flag)
    specimux_trace_TIMESTAMP_WORKER.tsv  # Detailed processing events per worker
  
  log.txt                          # Complete console output log with run parameters and results
```

#### Directory Structure Benefits

The match-type-first organization provides several advantages:

1. **Primary Data Access**: Full matches are immediately accessible in the `full/` directory without navigating through multiple subdirectories
2. **Clean Separation**: Partial matches and unknowns are segregated, reducing clutter when accessing your primary demultiplexed data
3. **Convenient Aggregation**: Pool-level directories (e.g., `full/ITS/`) collect all successful matches for that target region
4. **Recovery Options**: The `partial/` directory contains sequences that may be recoverable using tools like `specimine.py`
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
   - Multiple equivalent matches handled with configurable strategies

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

### Quality-Based Subsampling

The `--sample-topq N` option creates subsampled datasets containing only the highest-quality sequences:

```bash
specimux primers.fasta specimens.txt sequences.fastq -F -O output --sample-topq 500
```

Features:
- Creates `subsample/` directory mirroring the structure of `full/`
- Sorts sequences by average Phred quality score
- Retains only the top N sequences from each file
- Preserves primer files (primers.fasta and primers.txt) in subsample directories
- Runs as post-processing step after demultiplexing completes

Use cases:
- Generate high-confidence datasets for downstream analysis
- Create smaller representative datasets for testing
- Focus computational resources on highest-quality reads
- Compatibility with tools like NGSpeciesID that benefit from quality filtering

Example output structure:
```
output/
├── full/           # Complete demultiplexed sequences
│   └── ITS/
│       └── specimen_001.fastq (1000 sequences)
└── subsample/      # Top 500 highest-quality sequences
    └── ITS/
        └── specimen_001.fastq (500 sequences)
```

### Live Sequencing with specimux-watch

For live MinKNOW sequencing workflows, `specimux-watch` automatically monitors a directory and processes new FASTQ files as they are written:

```bash
specimux-watch primers.fasta specimens.txt /path/to/minknow/output -F -O demux_output/ -d
```

**Key features:**

- **Automatic detection**: Monitors directory for new `.fastq` files as MinKNOW writes them
- **File stability checking**: Waits for files to finish writing before processing (default: 30s settle time)
- **Sequential processing**: Ensures only one file is processed at a time to avoid resource conflicts
- **Real-time output**: Progress bars and logs display in real-time during processing
- **State persistence**: Tracks processed files to avoid reprocessing if restarted
- **Safe restarts**: Ignores pre-existing files on startup, only processes new arrivals

**Common options:**

```bash
# Basic live sequencing
specimux-watch primers.fasta specimens.txt watch_dir/ -F -O output/

# Custom settle time for large files
specimux-watch primers.fasta specimens.txt watch_dir/ -F -O output/ --settle-time 60

# With diagnostics and specific file pattern
specimux-watch primers.fasta specimens.txt watch_dir/ -F -O output/ -d --pattern "*.fastq"

# Run as background daemon (logs to file)
specimux-watch primers.fasta specimens.txt watch_dir/ -F -O output/ --daemon
```

**Behavior:**

- On startup, all existing `.fastq` files in the watch directory are marked as "ignored" and not processed
- Only files that arrive **after** `specimux-watch` starts are automatically processed
- Each successfully processed file is recorded in a state file (`.specimux-watch-state.json`)
- If you need to (re)process an existing file, run `specimux` on it directly

**Use cases:**

- Live demultiplexing during long sequencing runs
- Processing files as they complete writing
- Automated pipeline integration
- Continuous monitoring of sequencing output

All standard `specimux` arguments (edit distances, trimming modes, diagnostics, etc.) are supported and passed through to the demultiplexer.

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

### Run Logging

Specimux provides comprehensive logging during processing:
- Real-time progress updates with sequence counts and match rates
- Pool configuration validation and statistics
- Final processing summary (total sequences, match rate, processing time)
- Command line parameters used for the run
- Error reporting and diagnostic information

When using file output (`-F` flag), all console output is automatically duplicated to `log.txt` in the output directory, providing a permanent record of each processing run.

### Trace Logging System (-d)

The diagnostic mode now provides comprehensive trace logging for detailed pipeline analysis:

**Verbosity Levels:**
- `-d` or `-d1`: Standard events (match results, decisions, outputs)
- `-d2`: Detailed events including successful search attempts  
- `-d3`: Verbose events including all search attempts (successful and failed)

**Trace Files:**
- Created in `output_dir/trace/` directory
- One TSV file per worker process: `specimux_trace_TIMESTAMP_WORKER.tsv`
- Contains timestamped events tracking each sequence through the pipeline

**Key Events Logged:**
- Sequence received/filtered/output decisions
- Primer/barcode search attempts and matches
- Multiple match detection and resolution
- Specimen identification and pool assignment
- Match scoring and selection logic

This system enables detailed analysis of processing efficiency, match patterns, and troubleshooting of specific sequences.

**For complete trace event documentation, see [trace_event_schema.md](trace_event_schema.md).**

### Trace-Based Statistics and Visualization

The trace system enables comprehensive post-processing analysis through two complementary tools:

#### specimux-stats - Flexible Statistics Engine

Converts trace events into statistical summaries with any combination of analysis dimensions:

```bash
# Hierarchical text analysis
specimux-stats trace/ --hierarchical pool primer_pair outcome
specimux-stats trace/ --hierarchical orientation match_type --count-by sequences

# Export data for visualization  
specimux-stats trace/ --sankey-data pool outcome --output flow.json

# List all available dimensions
specimux-stats trace/ --list-dimensions

# Classification diagnostics (similar to v0.5 classification system)
specimux-stats trace/ --hierarchical pool primer_pair match_type --count-by sequences
```

To obtain similar diagnostic information as the v0.5 classification system, use the last command above. This provides a biologically meaningful breakdown showing exactly which primers and barcodes were detected for each sequence, organized by pool and primer pair.

**Available dimensions:** `orientation`, `pool`, `primer_pair`, `forward_primer`, `reverse_primer`, `forward_barcode`, `reverse_barcode`, `barcode_count`, `match_type`, `outcome`, `selection_strategy`, `discard_reason`, `outcome_detailed`, and more.

**Counting modes:**
- `candidate_matches`: Count every primer-pair match attempt (detailed pipeline analysis)  
- `sequences`: Count unique sequences only (overall success rates)

#### specimux-visualize - Interactive Flow Diagrams

Creates interactive Sankey diagrams from trace statistics:

```bash  
# Basic flow diagram
specimux-visualize flow.json diagram.html

# Custom styling
specimux-visualize flow.json diagram.html --theme dark --width 1600 --height 800
```

**Features:**
- Semantic coloring based on processing pipeline stages
- Interactive hover details with flow counts
- Automatic layout adaptation to data structure
- Support for arbitrary user-defined pools and dimensions

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

The trace-based statistics system can generate data for Sankey flow diagrams showing how sequences move through the processing pipeline.

### Visualization Tool

Use the `specimux-visualize` command to create interactive Sankey diagrams:

```bash
# Generate flow data first
specimux-stats trace/ --sankey-data pool outcome --output flow.json

# Create visualization
specimux-visualize flow.json my_flow_diagram.html

# Custom styling
specimux-visualize flow.json diagram.html --theme dark --width 1600 --height 800
```

### Dependencies

Visualization support is included by default with plotly>=5.0.0 dependency.

### Features

- **Interactive Diagrams**: Hover over nodes and flows to see exact counts
- **Semantic Coloring**: Automatic colors based on processing pipeline stages
- **Customizable**: Adjustable dimensions and themes for different display needs
- **Self-Contained**: Generated HTML files work offline and can be shared easily
- **Flexible Data**: Works with any combination of trace analysis dimensions

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

The `specimux-convert` command automatically:

1. Extracts all unique primer sequences
2. Generates a `primers.fasta` file with proper pool annotations
3. Creates a new specimen file with the required `PrimerPool` column
4. Replaces primer sequences with primer names

### Usage

```bash
specimux-convert Index.txt --output-specimen=IndexPP.txt --output-primers=primers.fasta --pool-name=ITS
```

### Arguments

- `input_file`: The old format specimen file (required)
- `--output-specimen`: Path for the new format specimen file (default: specimen_new.txt)
- `--output-primers`: Path for the primers FASTA file (default: primers.fasta)
- `--pool-name`: Name to use for the primer pool (default: pool1)

## Version History
- 0.7.0 (January 2026): Add dereplication to prevent artificial read amplification. When complex primer pools cause the same read to match multiple primer permutations, specimux now selects the best match per specimen/barcode group by default. New `--dereplicate` option replaces `--resolve-multiple-matches` with values `best` (default) and `none`. The `downgrade-full` strategy has been removed. Also fixes a bug where ~180 reads per run were silently dropped when trimming would produce empty sequences (now routed to unknown output).
- 0.6.9 (December 2025): Fix specimine path derivation for current output structure. The tool now correctly finds partial match files after the specimux output reorganization. Also fixes --partial-reverse flag (changed to --no-partial-reverse) which was broken due to argparse store_true with default=True
- 0.6.8 (December 2025): Add validation for empty barcodes in specimen file. Single-indexed demultiplexing (where FwIndex or RvIndex is empty) is not supported and now produces a clear error message listing affected specimens instead of crashing during alignment
- 0.6.7 (November 2025): Fix crash when using --trim with --sample-topq on sequences where trimming produces empty result. Very short sequences with overlapping primers now skip output instead of writing empty records that caused division by zero during subsampling
- 0.6.6 (October 2025): Fix pool assignment bug for sequences matching primers shared across multiple pools. Full matches now correctly use the specimen's declared pool from Index.txt rather than ambiguous primer-based pool selection. Pool selection is now deterministic (alphabetical) for partial matches and edge cases. This fixes incorrect routing where specimens declared in one pool were being output to a different pool when primers belonged to multiple pools
- 0.6.5 (October 2025): Fix bug in --disable-prefilter flag where code attempted to call .match() on None prefilter object, causing AttributeError. Updated type hints to Optional[BarcodePrefilter] and added None check before prefilter usage
- 0.6.4 (October 2025): Change default output file prefix from "sample_" to empty string for cleaner filenames. All tools (specimux, specimux-watch, specimine) now produce files like "specimen_001.fastq" instead of "sample_specimen_001.fastq". Backward compatible with legacy "sample_" prefixed files. Users can still specify custom prefix with -P flag
- 0.6.3 (October 2025): Add specimux-watch for live MinKNOW sequencing workflows with automatic file monitoring and processing. Fix duplicate output bug when primers belong to multiple pools. Pool assignment now uses the attempted primer pair context rather than matched primers. Fix empty directory pruning to ignore primer metadata files when determining if a directory should be removed. Update validation script to handle sequences appearing in multiple locations
- 0.6.2 (August 2025): Fix primer orientation detection bug introduced in commit f0209a3 (January 29, 2025). The determine_orientation function now correctly searches for reverse primers at the beginning of the reverse complement sequence, properly detecting sequence orientation for pre-filtering
- 0.6.1 (August 2025): Fix sequence orientation normalization bug introduced on August 11, 2025. Sequences are now properly normalized to canonical orientation regardless of input orientation, ensuring consistent output for the same biological sequences
- 0.6.0-dev (August 2025): Modern Python packaging with pip installation support, Python 3.10-3.13 compatibility with maintained bloom filter dependency (pybloomfilter3), dedicated CLI commands for all tools, major code refactoring with modular architecture, multiple match processing (replacing "ambiguity" concept), reorganized output with match-type-first directory structure for easier access to primary data, comprehensive trace event system with 3 verbosity levels, trace-based statistics framework with hierarchical analysis capabilities, interactive Sankey flow diagrams, automatic cleanup of empty directories
- 0.5.1 (March 2025): Primer Pools implementation with hierarchical output and pool-level full match collections, detailed run logging with log.txt files
- 0.4 (February 2025): Added Bloom filter optimization for performance improvements
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