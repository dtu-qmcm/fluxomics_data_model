# IO

The IO module provides parsers and writers for different formats.

## Convenience functions

::: fluxomics_data_model.io.parse_fluxml_file
    options:
      show_root_heading: true

::: fluxomics_data_model.io.write_fluxml
    options:
      show_root_heading: true

::: fluxomics_data_model.io.parse_mtf
    options:
      show_root_heading: true

::: fluxomics_data_model.io.write_mtf
    options:
      show_root_heading: true

::: fluxomics_data_model.io.parse_freeflux
    options:
      show_root_heading: true

## Parsers

::: fluxomics_data_model.io.fluxml_parser.FluxMLParser
    options:
      show_root_heading: true
      members:
        - parse_file

::: fluxomics_data_model.io.mtf_parser.MTFParser
    options:
      show_root_heading: true
      members:
        - parse

::: fluxomics_data_model.io.freeflux_parser.FreefluxParser
    options:
      show_root_heading: true
      members:
        - parse

## Writers

::: fluxomics_data_model.io.fluxml_writer.FluxMLWriter
    options:
      show_root_heading: true
      members:
        - write

::: fluxomics_data_model.io.mtf_writer.MTFWriter
    options:
      show_root_heading: true
      members:
        - write
