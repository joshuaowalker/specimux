#!/usr/bin/env python3

"""
Database and registry classes for specimux.

This module contains database classes for managing primers and specimens,
as well as the barcode prefilter protocol and simple implementations.
"""

import logging
from typing import Dict, List, Optional, Set, Protocol

from .constants import Primer
from .models import PrimerInfo


class PrimerDatabase:
    """Manages primers and their relationships to pools."""

    def __init__(self):
        self._primers: Dict[str, PrimerInfo] = {}  # name -> PrimerInfo
        self._pools: Dict[str, Set[str]] = {}  # pool -> set of primer names
        self._pool_primers: Dict[str, Dict[Primer, List[PrimerInfo]]] = {}  # pool -> {direction -> [PrimerInfo]}

    def add_primer(self, primer: PrimerInfo, pools: List[str]) -> None:
        """
        Add a primer to the registry and its pools.

        Args:
            primer: PrimerInfo object to add
            pools: List of pool names this primer belongs to

        Raises:
            ValueError: If primer name already exists
        """
        if primer.name in self._primers:
            raise ValueError(f"Duplicate primer name: {primer.name}")

        # Store primer
        self._primers[primer.name] = primer

        # Add to pools
        for pool in pools:
            if pool not in self._pools:
                self._pools[pool] = set()
                self._pool_primers[pool] = {
                    Primer.FWD: [],
                    Primer.REV: []
                }
            self._pools[pool].add(primer.name)
            self._pool_primers[pool][primer.direction].append(primer)

    def get_primer(self, name: str) -> Optional[PrimerInfo]:
        """Get a primer by name."""
        return self._primers.get(name)

    def get_primers_in_pool(self, pool: str) -> List[PrimerInfo]:
        """Get all primers in a pool."""
        if pool not in self._pools:
            return []
        primers = []
        for direction in [Primer.FWD, Primer.REV]:
            primers.extend(self._pool_primers[pool][direction])
        return primers

    def get_pools(self) -> List[str]:
        """Get list of all pool names."""
        return list(self._pools.keys())

    def get_pool_primers(self, pool: str, direction: Optional[Primer] = None) -> List[PrimerInfo]:
        """
        Get primers in a pool, optionally filtered by direction.

        Args:
            pool: Name of the pool
            direction: Optional primer direction to filter by

        Returns:
            List of PrimerInfo objects
        """
        if pool not in self._pools:
            return []
        if direction:
            return self._pool_primers[pool][direction]
        return self.get_primers_in_pool(pool)

    def primer_in_pool(self, primer_name: str, pool: str) -> bool:
        """Check if a primer is in a pool."""
        return pool in self._pools and primer_name in self._pools[pool]

    def validate_pools(self) -> None:
        """
        Validate pool configurations.

        Raises:
            ValueError: If validation fails
        """
        for pool in self._pools:
            # Each pool must have at least one forward and one reverse primer
            if not self._pool_primers[pool][Primer.FWD]:
                raise ValueError(f"Pool {pool} has no forward primers")
            if not self._pool_primers[pool][Primer.REV]:
                raise ValueError(f"Pool {pool} has no reverse primers")

    def get_pool_stats(self) -> Dict:
        """Get statistics about pools and primers."""
        stats = {
            'total_primers': len(self._primers),
            'total_pools': len(self._pools),
            'pools': {}
        }

        for pool in self._pools:
            stats['pools'][pool] = {
                'forward_primers': len(self._pool_primers[pool][Primer.FWD]),
                'reverse_primers': len(self._pool_primers[pool][Primer.REV]),
                'total_primers': len(self.get_primers_in_pool(pool))
            }

        return stats


class Specimens:
    """Manages specimen information including barcodes and primers."""
    
    def __init__(self, primer_registry: PrimerDatabase):
        self._specimens = []  # List of (id, pool, b1, p1s, b2, p2s) tuples for reference
        self._barcode_length = 0
        self._primers = {}
        self._specimen_ids = set()
        self._primer_pairings = {}
        self._primer_registry = primer_registry
        self._active_pools = set()  # Track pools actually used by specimens

    def add_specimen(self, specimen_id: str, pool: str, b1: str, p1: str, b2: str, p2: str):
        """Add a specimen with its barcodes and primers."""
        if specimen_id in self._specimen_ids:
            raise ValueError(format(f"Duplicate specimen id in index file: {specimen_id}"))
        self._specimen_ids.add(specimen_id)

        # Track active pools
        self._active_pools.add(pool)

        self._barcode_length = max(self._barcode_length, len(b1), len(b2))

        # Handle wildcards and get list of possible primers
        p1_list = self._resolve_primer_name(p1, pool, Primer.FWD)
        p2_list = self._resolve_primer_name(p2, pool, Primer.REV)

        # Register primers and barcodes
        for p1_info in p1_list:
            ps1 = p1_info.primer
            if ps1 not in self._primers:
                self._primers[ps1] = p1_info
            primer_info = self._primers[ps1]
            primer_info.barcodes.add(b1)
            primer_info.specimens.add(specimen_id)

        for p2_info in p2_list:
            ps2 = p2_info.primer
            if ps2 not in self._primers:
                self._primers[ps2] = p2_info
            primer_info = self._primers[ps2]
            primer_info.barcodes.add(b2)
            primer_info.specimens.add(specimen_id)

        self._specimens.append((specimen_id, pool, b1, p1_list, b2, p2_list))

    def prune_unused_pools(self):
        """Remove pools that aren't used by any specimens."""
        # Get list of all unused pools
        all_pools = set(self._primer_registry.get_pools())
        unused_pools = all_pools - self._active_pools

        if unused_pools:
            logging.info(f"Removing unused pools: {unused_pools}")

            # Remove unused pools from all primers
            for primer in self._primers.values():
                primer.pools = [p for p in primer.pools if p in self._active_pools]

            # Update primer registry
            for pool in unused_pools:
                # Remove pool from registry's internal data structures
                if pool in self._primer_registry._pools:
                    del self._primer_registry._pools[pool]
                if pool in self._primer_registry._pool_primers:
                    del self._primer_registry._pool_primers[pool]

            # Update pool stats in logs
            stats = self._primer_registry.get_pool_stats()
            logging.info(f"After pruning: {stats['total_primers']} primers in {stats['total_pools']} pools")
            for pool, pool_stats in stats['pools'].items():
                logging.info(f"Pool {pool}: {pool_stats['forward_primers']} forward, "
                             f"{pool_stats['reverse_primers']} reverse primers")

    def _resolve_primer_name(self, primer_name: str, pool: str, direction: Primer) -> List[PrimerInfo]:
        """Resolve a primer name (including wildcards) to a list of PrimerInfo objects."""
        if primer_name == '-' or primer_name == '*':  # Handle wildcards
            # Get all primers in the specified pool and direction
            primers = []
            for p in self._primer_registry.get_primers_in_pool(pool):
                if p.direction == direction:
                    primers.append(p)
            if not primers:
                raise ValueError(f"No {direction.name} primers found in pool {pool}")
            return primers
        else:
            # Get specific primer
            primer = self._primer_registry.get_primer(primer_name)
            if not primer:
                raise ValueError(f"Primer not found: {primer_name}")
            if primer.direction != direction:
                raise ValueError(f"Primer {primer_name} is not a {direction.name} primer")
            if not self._primer_registry.primer_in_pool(primer_name, pool):
                raise ValueError(f"Primer {primer_name} is not in pool {pool}")
            return [primer]

    def specimens_for_barcodes_and_primers(self, b1_list: List[str], b2_list: List[str],
                                           p1_matched: PrimerInfo, p2_matched: PrimerInfo) -> List[str]:
        """Find specimens matching given barcodes and primers."""
        matching_specimens = []
        for spec_id, pool, b1, p1s, b2, p2s in self._specimens:
            if (p1_matched in p1s and
                    p2_matched in p2s and
                    b1.upper() in b1_list and
                    b2.upper() in b2_list):
                matching_specimens.append(spec_id)

        return matching_specimens

    def specimen_for_exact_match(self, b1: str, b2: str,
                                 p1: PrimerInfo, p2: PrimerInfo) -> Optional[str]:
        """Find specimen matching exact barcode and primer combination.

        Unlike specimens_for_barcodes_and_primers() which accepts lists of barcodes,
        this returns the specimen for a single specific barcode+primer combination.
        Used by dereplication to map specific barcode choices to specimens.
        """
        for spec_id, pool, spec_b1, p1s, spec_b2, p2s in self._specimens:
            if (p1 in p1s and p2 in p2s and
                    spec_b1.upper() == b1.upper() and
                    spec_b2.upper() == b2.upper()):
                return spec_id
        return None

    def get_primers(self, direction: Primer) -> List[PrimerInfo]:
        """Get all primers in a given direction."""
        return [p for p in self._primers.values() if p.direction == direction]

    def get_paired_primers(self, primer: str) -> List[PrimerInfo]:
        """Get primers that are paired with the given primer."""
        if primer in self._primer_pairings:
            return self._primer_pairings[primer]

        specimens = self._primers[primer].specimens
        direction = self._primers[primer].direction
        rv = []
        for pi in self._primers.values():
            if direction != pi.direction and pi.specimens.intersection(specimens):
                rv.append(pi)

        self._primer_pairings[primer] = rv
        return rv

    def get_specimen_pool(self, specimen_id: str) -> Optional[str]:
        """Get the primer pool for a specimen."""
        for spec_id, pool, _, _, _, _ in self._specimens:
            if spec_id == specimen_id:
                return pool
        return None

    def b_length(self):
        """Get the maximum barcode length."""
        return self._barcode_length

    def validate(self):
        """Validate specimen database consistency."""
        self._validate_barcodes_globally_unique()
        self._validate_barcode_lengths()
        self.prune_unused_pools()

    def _validate_barcodes_globally_unique(self):
        """Check if barcodes are unique between forward and reverse."""
        all_b1s = set()
        all_b2s = set()
        for primer in self._primers.values():
            if primer.direction == Primer.FWD:
                all_b1s.update(primer.barcodes)
            else:
                all_b2s.update(primer.barcodes)
        dups = all_b1s.intersection(all_b2s)
        if len(dups) > 0:
            logging.warning(f"Duplicate Barcodes ({len(dups)}) in Fwd and Rev: {dups}")

    def _validate_barcode_lengths(self):
        """Check if all barcodes have consistent lengths."""
        all_b1s = set()
        all_b2s = set()
        for primer in self._primers.values():
            if primer.direction == Primer.FWD:
                all_b1s.update(primer.barcodes)
            else:
                all_b2s.update(primer.barcodes)
        if len(set(len(b) for b in all_b1s)) > 1:
            logging.warning("Forward barcodes have inconsistent lengths")
        if len(set(len(b) for b in all_b2s)) > 1:
            logging.warning("Reverse barcodes have inconsistent lengths")


class BarcodePrefilter(Protocol):
    """Protocol defining the interface for barcode prefilters."""
    
    def match(self, barcode: str, sequence: str) -> bool:
        """Check if barcode potentially matches sequence."""
        ...


class PassthroughPrefilter:
    """A prefilter that always returns True (no filtering)."""

    def match(self, barcode: str, sequence: str) -> bool:
        # always call edlib.align
        return True