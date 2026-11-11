"""Independent statistical cross-check of the frozen exosearch v1.0.0 reliability products.

This package does not replace the interim threshold or the star-quality config. It reads the same
injection and null products (or the sample CSVs) and writes its own JSON and figures.
"""

from .crosscheck import run_crosscheck

__all__ = ["run_crosscheck"]
