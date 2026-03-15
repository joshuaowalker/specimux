"""Optional JSONL progress reporting for orchestration tools."""

import json
import time
from pathlib import Path
from typing import Optional


class ProgressReporter:
    """Write JSONL progress lines to a file for external consumption.

    When path is None, all methods are no-ops (zero overhead).
    Progress lines are throttled to at most one per second.
    """

    def __init__(self, path: Optional[Path] = None):
        self._path = path
        self._file = None
        self._last_write = 0.0
        self._throttle = 1.0  # seconds between writes
        if path:
            self._file = open(path, 'w')

    def update(self, processed: int, matched: int, total_est: int) -> None:
        """Write a throttled progress line."""
        if self._file is None:
            return
        now = time.monotonic()
        if now - self._last_write < self._throttle:
            return
        self._last_write = now
        rate = matched / processed if processed > 0 else 0.0
        line = json.dumps({
            "type": "progress",
            "processed": processed,
            "matched": matched,
            "total_est": total_est,
            "rate": round(rate, 4),
        })
        self._file.write(line + '\n')
        self._file.flush()

    def complete(self, processed: int, matched: int) -> None:
        """Write the final completion line."""
        if self._file is None:
            return
        rate = matched / processed if processed > 0 else 0.0
        line = json.dumps({
            "type": "complete",
            "processed": processed,
            "matched": matched,
            "rate": round(rate, 4),
        })
        self._file.write(line + '\n')
        self._file.flush()

    def close(self) -> None:
        """Close the progress file."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
