"""
Shared pytest fixtures for specimux tests.
"""

import pytest
import tempfile
from pathlib import Path
import sys


@pytest.fixture(scope="session")
def package_root():
    """Return the root directory of the specimux package."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def test_data_root():
    """Return the root directory containing test data."""
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def integration_test_data(test_data_root):
    """Return path to integration test suite data."""
    path = test_data_root / "integration_test_suite"
    if not path.exists():
        pytest.skip("Integration test data not found")
    return path


@pytest.fixture
def temp_dir():
    """Provide a temporary directory that's cleaned up after the test."""
    with tempfile.TemporaryDirectory(prefix="specimux_test_") as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_primers():
    """Provide sample primer data for testing."""
    return """>ITS1F pool=ITS position=forward
CTTGGTCATTTAGAGGAAGTAA
>ITS4 pool=ITS position=reverse
TCCTCCGCTTATTGATATGC
>gITS7 pool=ITS2 position=forward
GTGARTCATCGARTCTTTG
>ITS4 pool=ITS2 position=reverse
TCCTCCGCTTATTGATATGC"""


@pytest.fixture
def sample_specimens():
    """Provide sample specimen data for testing."""
    return """SampleID\tPrimerPool\tFwIndex\tFwPrimer\tRvIndex\tRvPrimer
TEST_SPECIMEN_001\tITS2\tTGCATACGGTGTC\tgITS7\tATAATATTCGGCA\tITS4
TEST_SPECIMEN_002\tITS2\tTCTTACCTCATTC\tgITS7\tCGCATCCATACTG\tITS4
TEST_SPECIMEN_003\tITS\tAGCATGACCTGTT\tITS1F\tGGTTCTCGCATCG\tITS4"""


@pytest.fixture
def sample_sequence_fastq():
    """Provide a sample FASTQ sequence for testing."""
    return """@test_seq_001
CTTGGTCATTTAGAGGAAGTAAACGTACGTACGTACGTACGTACGTACGTACGTTCCTCCGCTTATTGATATGC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"""


# Markers for test organization
def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as an integration test"
    )
    config.addinivalue_line(
        "markers", "unit: mark test as a unit test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "requires_data: mark test as requiring external data"
    )