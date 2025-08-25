"""Specimux: Demultiplexing tools for MinION sequenced reads."""

__version__ = "0.6.0-dev"

# Re-export key functions and classes that might be useful for programmatic access
from .core import (
    PrimerDatabase,
    Specimens,
    MatchParameters,
    specimux,
    specimux_mp,
)

__all__ = [
    "PrimerDatabase",
    "Specimens", 
    "MatchParameters",
    "specimux",
    "specimux_mp",
]