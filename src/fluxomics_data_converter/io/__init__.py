"""Fluxomics I/O operations package.

Contains parsers and writers for three fluxomics data formats:

- **FluxML** — hierarchical XML format used by 13CFlux/13CFlux2
- **MTF** — multi-file text format used by influx_si
- **Freeflux** — tabular format (tsv/csv/xlsx) used by Freeflux

Quick start
-----------
Auto-detect format from file extension and parse::

    from fluxomics_data_converter.io import parse
    model = parse("models/ecoli.fml")          # FluxML
    model = parse("models/ecoli.netw")         # MTF  (pass the .netw file)
    model = parse("models/reactions.tsv")      # Freeflux (pass reactions file)

Or use format-specific convenience functions::

    from fluxomics_data_converter.io import parse_fluxml_file, parse_mtf, parse_freeflux
    model = parse_fluxml_file("models/ecoli.fml")
    model = parse_mtf("models/ecoli")          # base path without extension
    model = parse_freeflux("models/ecoli/")    # directory containing reactions.*

Format protocols
----------------
All parsers satisfy :class:`~fluxomics_data_converter.io.base.FluxomicsParser` and
all writers satisfy :class:`~fluxomics_data_converter.io.base.FluxomicsWriter`.
See those protocols for the expected interface when writing format-agnostic code.
"""

from pathlib import Path

from .base import FluxomicsParser, FluxomicsWriter
from .fluxml_parser import FluxMLParser, parse_fluxml_file
from .fluxml_writer import FluxMLWriter, write_fluxml
from .mtf_parser import MTFParser, parse_mtf
from .mtf_writer import MTFWriter, write_mtf
from .freeflux_parser import FreefluxParser, parse_freeflux
from .freeflux_writer import FreefluxWriter, write_freeflux

# Extension → parser class mapping used by parse()
_PARSERS: dict[str, type] = {
    ".fml": FluxMLParser,
    ".xml": FluxMLParser,
    ".netw": MTFParser,
    ".tsv": FreefluxParser,
    ".csv": FreefluxParser,
    ".xlsx": FreefluxParser,
}


def parse(path: str | Path):
    """Auto-detect format from file extension and parse.

    This is the single entry point for reading any supported fluxomics
    format.  The format is inferred from the file extension:

    - ``.fml`` / ``.xml``             → :class:`FluxMLParser`
    - ``.netw``                        → :class:`MTFParser`
      (pass the ``.netw`` file; all associated ``.linp``, ``.miso`` etc.
      files are resolved automatically from the same directory)
    - ``.tsv`` / ``.csv`` / ``.xlsx`` → :class:`FreefluxParser`
      (pass the ``reactions.*`` file or the directory containing it)

    Args:
        path: Path to the model file (or base path for multi-file formats).

    Returns:
        A fully validated :class:`~fluxomics_data_converter.FluxomicsData`.

    Raises:
        ValueError: If the file extension is not recognised.
        FileNotFoundError: If *path* does not exist.

    Examples::

        from fluxomics_data_converter.io import parse

        model = parse("data/ecoli.fml")
        model = parse("data/ecoli.netw")
        model = parse("data/reactions.tsv")
    """
    path = Path(path)
    parser_cls = _PARSERS.get(path.suffix.lower())
    if parser_cls is None:
        supported = sorted(_PARSERS.keys())
        raise ValueError(
            f"Unrecognised file extension '{path.suffix}'. "
            f"Supported extensions: {supported}"
        )
    return parser_cls().parse(path)


__all__ = [
    # Protocols
    "FluxomicsParser",
    "FluxomicsWriter",
    # Unified entry point
    "parse",
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
