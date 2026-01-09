"""Specimux: Demultiplexing tools for MinION sequenced reads."""

__version__ = "0.7.0"

# Import key classes and functions from the refactored modules
from .databases import PrimerDatabase, Specimens
from .models import MatchParameters
from .orchestration import specimux, specimux_mp
from .io_utils import read_primers_file, read_specimen_file
from .orchestration import setup_match_parameters
from .demultiplex import process_sequences
from .trace import TraceLogger

__all__ = [
    "PrimerDatabase",
    "Specimens",
    "MatchParameters",
    "specimux",
    "specimux_mp",
    "read_primers_file",
    "read_specimen_file",
    "setup_match_parameters",
    "process_sequences",
    "TraceLogger",
]
