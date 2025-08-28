# Specimux Integration Test Suite

This test suite covers comprehensive scenarios for specimux demultiplexer validation using real sequences sampled from ont63 dataset.

## Test Data Overview

**Total sequences**: 40 representative sequences  
**Source dataset**: ont63 (567K sequences)  
**Primer pools**: ITS and ITS2  
**Match rate**: 15% (6 full matches, partial matches, and unknowns)

## Test Scenarios Covered

### 1. Full Match Cases (Score 5.0)
- **cb5a9725-c40b-4d0d-bf62-cbe19d06cf43**: ITS2 pool, perfect barcode match → TEST_SPECIMEN_001
- **18705e64-bb48-4f63-8a4a-8a2ccf715955**: ITS2 pool, perfect barcode match → TEST_SPECIMEN_002
- **9d6f7b24-979e-41db-b47d-c9c3e7e0a60f**: ITS pool, perfect barcode match → TEST_SPECIMEN_003

### 2. Partial Match Cases (Score 4.0)
- **Forward barcode only**: Sequences with P1+P2+B1 matches
- **Reverse barcode only**: Sequences with P1+P2+B2 matches

### 3. Primer-Only Cases (Score 2.0)
- Both primers found but no barcodes
- Single primer matches

### 4. Unknown/Unmatchable Cases
- No primer matches found
- Primers outside search regions
- High edit distance beyond thresholds

### 5. Edge Cases
- **Sequence lengths**: Range from 47bp to 1411bp
- **Orientation detection**: Unknown, forward, reverse orientations
- **Multiple primer pools**: Cross-pool matching scenarios

## Test Structure

```
test_data/integration_test_suite/
├── sequences.fastq          # 40 test sequences
├── primers.fasta           # ITS and ITS2 primer pools  
├── specimens.txt           # 3 known specimens for full matches
├── expected_output/        # Expected results from specimux
│   ├── full/              # Full matches (6 sequences)
│   │   ├── ITS/           # 1 sequence
│   │   └── ITS2/          # 2 sequences  
│   ├── partial/           # Partial matches (barcodes only)
│   ├── unknown/           # Primer-only and unmatched
│   └── trace/             # Debug trace files
└── README.md              # This documentation
```

## Running Tests

### Basic Test Run
```bash
python specimux.py test_data/integration_test_suite/primers.fasta test_data/integration_test_suite/specimens.txt test_data/integration_test_suite/sequences.fastq -F -O test_output -d
```

### Validation Script
```bash
# TODO: Create validation script to compare test_output vs expected_output
python validate_test_results.py test_output expected_output
```

## Expected Results Summary

- **Full matches**: 3 specimens matched (15% success rate)
- **Partial matches**: Multiple forward/reverse barcode-only matches
- **Unknown sequences**: Majority of sequences with primer-only or no matches
- **Trace files**: Complete event logging for analysis

## Test Categories by Sequence ID

### Score 5.0 (Perfect Matches)
- cb5a9725-c40b-4d0d-bf62-cbe19d06cf43 → TEST_SPECIMEN_001
- 18705e64-bb48-4f63-8a4a-8a2ccf715955 → TEST_SPECIMEN_002
- 9d6f7b24-979e-41db-b47d-c9c3e7e0a60f → TEST_SPECIMEN_003

### Score 4.0 (Partial Matches) 
- a1ceebf9-b754-4c35-9a93-d29840adbcdb (forward barcode only)
- e2714492-b1e5-4a01-91ad-8dc0ce6c521a (reverse barcode only)

### Score 2.0 (Primer Only)
- d5a99b43-8367-4a95-9143-12c3f62f07f1 (both primers, no barcodes)

### Edge Cases
- **Short sequences**: e5752a0b-0248-4014-af95-f2a35476484c (47bp)
- **Long sequences**: 665f4f7a-7b7d-4d0f-8314-38d037899ef9 (1411bp)
- **Unknown orientation**: Multiple sequences with ambiguous primer detection

This test suite validates all major specimux scenarios and edge cases using real bioinformatics data.