"""
Fluxomics I/O operations package.

Contains parsers and writers for various fluxomics data formats:
- FluxML: XML-based format used by 13CFlux/13CFlux2
- MTF: Multi-file text format used by influx_si
- Freeflux: Tabular format (tsv/csv/xlsx) used by Freeflux
"""

from .fluxml_parser import FluxMLParser, parse_fluxml_file
from .fluxml_writer import FluxMLWriter, write_fluxml
from .mtf_parser import MTFParser, parse_mtf
from .mtf_writer import MTFWriter, write_mtf
from .freeflux_parser import FreefluxParser, parse_freeflux
from .freeflux_writer import FreefluxWriter, write_freeflux

__all__ = [
    # FluxML format
    "FluxMLParser",
    "parse_fluxml_file",
    "FluxMLWriter",
    "write_fluxml",
    # MTF format (influx_si)
    "MTFParser",
    "parse_mtf",
    "MTFWriter",
    "write_mtf",
    # Freeflux format (tabular)
    "FreefluxParser",
    "parse_freeflux",
    "FreefluxWriter",
    "write_freeflux",
]
