# Specimux Trace Event Schema Documentation

## Overview

This document defines the schema for trace events logged by Specimux when diagnostic mode (`-d`) is enabled. Each worker process writes events to a separate TSV file, capturing the complete journey of sequences through the demultiplexing pipeline.

## Event Verbosity Levels

Trace events support three verbosity levels controlled by the `-d` flag:
- `-d` or `-d1`: **Standard** - Key events only (match results, decisions, outputs)
- `-d2`: **Detailed** - Includes successful search attempts
- `-d3`: **Verbose** - Includes all search attempts (successful and failed)

## Base Fields (All Events)

Every event contains these base fields:
```
timestamp | worker_id | event_seq | sequence_id | event_type | ...event-specific-fields...
```

### Base Field Definitions
- `timestamp`: ISO 8601 format (YYYY-MM-DDTHH:MM:SS.mmm)
- `worker_id`: String, unique worker identifier (e.g., "worker_0", "worker_1", "main")
- `event_seq`: Integer, incrementing counter per worker to ensure event ordering
- `sequence_id`: String, globally unique sequence identifier (see below)
- `event_type`: String, event type identifier (see Event Types)

### Candidate Match Tracking

Some events track individual primer pair match attempts (candidate matches) against a sequence. These events include an additional `candidate_match_id` field to enable aggregation by match location rather than just by sequence:

- `candidate_match_id`: String, unique identifier for a specific primer pair match attempt within a sequence (format: "{sequence_id}_match_{counter}", e.g., "seq123#00042#worker_1_match_0")

### Sequence ID Format
The `sequence_id` is constructed to be both globally unique and consistent across runs:
```
{original_seq_id}#{file_record_num}#{worker_id}
```
Example: `@M00123:456:000000000-ABCDE:1:1101:15847:1234#00001523#worker_2`

## Standard Events (Verbosity Level 1)

### SEQUENCE_RECEIVED
**Purpose**: Log when a sequence enters the pipeline
```
timestamp | worker_id | event_seq | sequence_id | SEQUENCE_RECEIVED | sequence_length | sequence_name
```
- `sequence_length`: Integer, length of input sequence
- `sequence_name`: String, original sequence identifier from FASTQ/FASTA

### SEQUENCE_FILTERED  
**Purpose**: Log when a sequence is filtered out before processing
```
timestamp | worker_id | event_seq | sequence_id | SEQUENCE_FILTERED | sequence_length | filter_reason
```
- `sequence_length`: Integer, length of filtered sequence
- `filter_reason`: Enum: `too_short`, `too_long`

### ORIENTATION_DETECTED
**Purpose**: Log orientation detection result
```
timestamp | worker_id | event_seq | sequence_id | ORIENTATION_DETECTED | orientation | forward_score | reverse_score | confidence
```
- `orientation`: Enum: `forward`, `reverse`, `unknown`
- `forward_score`: Integer, number of forward primer matches found
- `reverse_score`: Integer, number of reverse primer matches found
- `confidence`: Float, confidence score (e.g., ratio of scores)

### PRIMER_MATCHED
**Purpose**: Log successful primer match (any configuration) for a specific candidate match
```
timestamp | worker_id | event_seq | sequence_id | PRIMER_MATCHED | candidate_match_id | match_type | forward_primer | reverse_primer | forward_distance | reverse_distance | pool | orientation_used
```
- `candidate_match_id`: String, unique identifier for this primer pair match attempt
- `match_type`: Enum: `both`, `forward_only`, `reverse_only`
- `forward_primer`: String, forward primer name (or "none")
- `reverse_primer`: String, reverse primer name (or "none")
- `forward_distance`: Integer, edit distance for forward primer (-1 if not found)
- `reverse_distance`: Integer, edit distance for reverse primer (-1 if not found)
- `pool`: String, primer pool name
- `orientation_used`: Enum: `as_is`, `reverse_complement`

### BARCODE_MATCHED
**Purpose**: Log barcode match result for a specific candidate match
```
timestamp | worker_id | event_seq | sequence_id | BARCODE_MATCHED | candidate_match_id | match_type | forward_barcode | reverse_barcode | forward_distance | reverse_distance | forward_primer | reverse_primer
```
- `candidate_match_id`: String, unique identifier for this primer pair match attempt
- `match_type`: Enum: `both`, `forward_only`, `reverse_only`, `none`
- `forward_barcode`: String, forward barcode name (or "none")
- `reverse_barcode`: String, reverse barcode name (or "none")
- `forward_distance`: Integer, edit distance for forward barcode (-1 if not found)
- `reverse_distance`: Integer, edit distance for reverse barcode (-1 if not found)
- `forward_primer`: String, adjacent forward primer
- `reverse_primer`: String, adjacent reverse primer

### MATCH_SCORED
**Purpose**: Log match scoring for a specific candidate match
```
timestamp | worker_id | event_seq | sequence_id | MATCH_SCORED | candidate_match_id | forward_primer | reverse_primer | forward_barcode | reverse_barcode | total_edit_distance | barcode_presence | score
```
- `candidate_match_id`: String, unique identifier for this primer pair match attempt
- `forward_primer`: String, forward primer name
- `reverse_primer`: String, reverse primer name
- `forward_barcode`: String, forward barcode (or "none")
- `reverse_barcode`: String, reverse barcode (or "none")
- `total_edit_distance`: Integer, sum of all edit distances
- `barcode_presence`: Enum: `both`, `forward_only`, `reverse_only`, `none`
- `score`: Float, calculated match score

### MATCH_SELECTED
**Purpose**: Log final match selection decision
```
timestamp | worker_id | event_seq | sequence_id | MATCH_SELECTED | selection_strategy | forward_primer | reverse_primer | forward_barcode | reverse_barcode | pool | is_unique
```
- `selection_strategy`: Enum: `unique`, `first`, `strict`, `all`
- `forward_primer`: String, selected forward primer
- `reverse_primer`: String, selected reverse primer
- `forward_barcode`: String, selected forward barcode (or "none")
- `reverse_barcode`: String, selected reverse barcode (or "none")
- `pool`: String, selected pool name
- `is_unique`: Boolean, whether selection was unique

### MATCH_DISCARDED
**Purpose**: Log when a candidate match is discarded during selection
```
timestamp | worker_id | event_seq | sequence_id | MATCH_DISCARDED | candidate_match_id | forward_primer | reverse_primer | forward_barcode | reverse_barcode | score | discard_reason
```
- `candidate_match_id`: String, unique identifier for the discarded match
- `forward_primer`: String, forward primer name (or "none")
- `reverse_primer`: String, reverse primer name (or "none")
- `forward_barcode`: String, forward barcode (or "none")
- `reverse_barcode`: String, reverse barcode (or "none")
- `score`: Float, the match score
- `discard_reason`: String, reason for discarding (e.g., "lower_score", "downgraded_multiple_full")

### SPECIMEN_RESOLVED
**Purpose**: Log specimen identification result
```
timestamp | worker_id | event_seq | sequence_id | SPECIMEN_RESOLVED | specimen_id | resolution_type | pool | forward_primer | reverse_primer | forward_barcode | reverse_barcode
```
- `specimen_id`: String, resolved specimen ID (or "UNKNOWN", "FWD_ONLY_XXX", "REV_ONLY_XXX")
- `resolution_type`: Enum: `full_match`, `partial_forward`, `partial_reverse`, `unknown`
- `pool`: String, pool name (or "unknown")
- `forward_primer`: String, forward primer used
- `reverse_primer`: String, reverse primer used
- `forward_barcode`: String, forward barcode matched (or "none")
- `reverse_barcode`: String, reverse barcode matched (or "none")

### SEQUENCE_OUTPUT
**Purpose**: Log final output decision
```
timestamp | worker_id | event_seq | sequence_id | SEQUENCE_OUTPUT | specimen_id | pool | primer_pair | file_path
```
- `specimen_id`: String, final specimen ID
- `pool`: String, pool name
- `primer_pair`: String, primer pair (format: "FWD-REV")
- `file_path`: String, relative path where sequence will be written (output type can be derived from path)

### NO_MATCH_FOUND
**Purpose**: Log when no viable matches found
```
timestamp | worker_id | event_seq | sequence_id | NO_MATCH_FOUND | stage_failed | reason
```
- `stage_failed`: Enum: `primer_search`, `barcode_search`, `specimen_lookup`
- `reason`: String, description of why matching failed

## Detailed Events (Verbosity Level 2+)

These events are only logged at higher verbosity levels.

### PRIMER_SEARCH (Level 2: successful only, Level 3: all attempts)
**Purpose**: Log individual primer search attempts
```
timestamp | worker_id | event_seq | sequence_id | PRIMER_SEARCH | primer_name | primer_direction | search_region_start | search_region_end | found | edit_distance | match_position
```
- `primer_name`: String, primer identifier (e.g., "ITS1F")
- `primer_direction`: Enum: `forward`, `reverse`
- `search_region_start`: Integer, start of search region in sequence
- `search_region_end`: Integer, end of search region in sequence
- `found`: Boolean, whether primer was found
- `edit_distance`: Integer, edit distance of match (-1 if not found)
- `match_position`: Integer, position where primer matched (-1 if not found)

### BARCODE_SEARCH (Level 3 only)
**Purpose**: Log individual barcode search attempts
```
timestamp | worker_id | event_seq | sequence_id | BARCODE_SEARCH | barcode_name | barcode_type | primer_adjacent | search_start | search_end | found | edit_distance | match_position
```
- `barcode_name`: String, barcode identifier (e.g., "BC01")
- `barcode_type`: Enum: `forward`, `reverse`
- `primer_adjacent`: String, name of adjacent primer
- `search_start`: Integer, start position of barcode search region
- `search_end`: Integer, end position of barcode search region
- `found`: Boolean, whether barcode was found
- `edit_distance`: Integer, edit distance of match (-1 if not found)
- `match_position`: Integer, position where barcode matched (-1 if not found)

## Field Conventions

### Special Values
- Use `-1` for numeric fields when value is not applicable/not found
- Use `none` for string fields when value is empty/not found
- Use `none` for enum fields when value cannot be determined

### Data Types
- **Timestamps**: ISO 8601 format (YYYY-MM-DDTHH:MM:SS.mmm)
- **Booleans**: `true` or `false` (lowercase)
- **Integers**: No quotes needed
- **Floats**: Decimal point required
- **Strings**: No quotes needed in TSV, tabs escaped as `\t`
- **Lists**: Comma-separated, no spaces (e.g., "ITS,RPB2,LSU")

## Example Event Sequences

### Successful Full Match (Level 1)
```
2025-01-10T10:15:23.456	worker_1	1	seq123#00042#worker_1	SEQUENCE_RECEIVED	1523	seq123
2025-01-10T10:15:23.457	worker_1	2	seq123#00042#worker_1	ORIENTATION_DETECTED	forward	2	0	1.0
2025-01-10T10:15:23.458	worker_1	3	seq123#00042#worker_1	PRIMER_MATCHED	seq123#00042#worker_1_match_0	both	ITS1F	ITS4	2	1	ITS	as_is
2025-01-10T10:15:23.459	worker_1	4	seq123#00042#worker_1	BARCODE_MATCHED	seq123#00042#worker_1_match_0	both	BC01	BC15	1	0	ITS1F	ITS4
2025-01-10T10:15:23.460	worker_1	5	seq123#00042#worker_1	MATCH_SCORED	seq123#00042#worker_1_match_0	ITS1F	ITS4	BC01	BC15	4	both	0.25
2025-01-10T10:15:23.461	worker_1	6	seq123#00042#worker_1	MATCH_SELECTED	unique	ITS1F	ITS4	BC01	BC15	ITS	true
2025-01-10T10:15:23.462	worker_1	7	seq123#00042#worker_1	SPECIMEN_RESOLVED	Sample_123	full_match	ITS	ITS1F	ITS4	BC01	BC15
2025-01-10T10:15:23.463	worker_1	8	seq123#00042#worker_1	SEQUENCE_OUTPUT	Sample_123	ITS	ITS1F-ITS4	full/ITS/ITS1F-ITS4/Sample_123.fastq
```

### Partial Match (Forward Barcode Only)
```
2025-01-10T10:15:24.100	worker_2	42	seq456#00100#worker_2	SEQUENCE_RECEIVED	1420	seq456
2025-01-10T10:15:24.101	worker_2	43	seq456#00100#worker_2	ORIENTATION_DETECTED	forward	1	1	0.5
2025-01-10T10:15:24.102	worker_2	44	seq456#00100#worker_2	PRIMER_MATCHED	seq456#00100#worker_2_match_0	both	RPB2-5F	RPB2-7R	1	2	RPB2	as_is
2025-01-10T10:15:24.103	worker_2	45	seq456#00100#worker_2	BARCODE_MATCHED	seq456#00100#worker_2_match_0	forward_only	BC05	none	0	-1	RPB2-5F	RPB2-7R
2025-01-10T10:15:24.104	worker_2	46	seq456#00100#worker_2	MATCH_SCORED	seq456#00100#worker_2_match_0	RPB2-5F	RPB2-7R	BC05	none	3	forward_only	0.5
2025-01-10T10:15:24.105	worker_2	47	seq456#00100#worker_2	MATCH_SELECTED	unique	RPB2-5F	RPB2-7R	BC05	none	RPB2	true
2025-01-10T10:15:24.106	worker_2	48	seq456#00100#worker_2	SPECIMEN_RESOLVED	FWD_ONLY_BC05	partial_forward	RPB2	RPB2-5F	RPB2-7R	BC05	none
2025-01-10T10:15:24.107	worker_2	49	seq456#00100#worker_2	SEQUENCE_OUTPUT	FWD_ONLY_BC05	RPB2	RPB2-5F-RPB2-7R	partial/RPB2/RPB2-5F-RPB2-7R/FWD_ONLY_BC05.fastq
```

## File Naming Convention

Trace files are named using the pattern:
```
trace_dir/specimux_trace_{timestamp}_{worker_id}.tsv
```

Example:
```
trace_dir/
  specimux_trace_20250110_101523_worker_0.tsv
  specimux_trace_20250110_101523_worker_1.tsv
  specimux_trace_20250110_101523_main.tsv
```

## Implementation Notes

1. **Buffering**: Events should be buffered and written in batches for performance
2. **File Rotation**: Consider rotating trace files if they exceed a size limit (e.g., 1GB)
3. **Compression**: Trace files compress well; consider gzip compression for archival
4. **Header**: Each trace file should start with a header row listing column names
5. **Thread Safety**: Each worker writes to its own file to avoid locking overhead

## Usage Examples

### Analyzing Match Rates by Pool
```sql
SELECT pool, 
       COUNT(*) as total_matches,
       SUM(CASE WHEN resolution_type = 'full_match' THEN 1 ELSE 0 END) as full_matches
FROM events
WHERE event_type = 'SPECIMEN_RESOLVED'
GROUP BY pool;
```

### Tracking Processing Time
```sql
SELECT sequence_id,
       MAX(timestamp) - MIN(timestamp) as processing_time_ms
FROM events
GROUP BY sequence_id
ORDER BY processing_time_ms DESC
LIMIT 10;
```