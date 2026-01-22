# FluxML Format Specification

## Overview

FluxML is an XML-based standard for 13C metabolic flux analysis models, developed by the ModSim group.

| Property | Value |
|----------|-------|
| **File Extension** | `.fml`, `.xml` |
| **Encoding** | UTF-8 |
| **Namespace** | `http://www.13cflux.net/fluxml` |
| **Schema** | https://13cflux.net/fluxml |
| **Documentation** | https://github.com/modsim/FluxML |

## Version History

| Version | Features |
|---------|----------|
| v1.0 | Isotopically stationary 13C-MFA |
| v1.1 | XML schema validation, flux ratios, error models, metadata attributes |
| v2.0 | Non-stationary 13C-MFA, pool sizes, time-based measurements, substrate profiles |
| v3.0 | Multi-tracer support (C, N, H, O, S), extended atom mapping |

## Document Structure

```xml
<?xml version="1.0" encoding="UTF-8"?>
<fluxml xmlns="http://www.13cflux.net/fluxml">
  <info>...</info>                          <!-- Optional: metadata -->
  <reactionnetwork>...</reactionnetwork>    <!-- Required: network definition -->
  <constraints>...</constraints>            <!-- Optional: global constraints -->
  <configuration name="exp1">...</configuration>  <!-- Optional: experiments -->
</fluxml>
```

## Namespace Handling

FluxML documents may use either default or prefixed namespaces:

```xml
<!-- Default namespace (recommended) -->
<fluxml xmlns="http://www.13cflux.net/fluxml">

<!-- Prefixed namespace (also valid) -->
<fml:fluxml xmlns:fml="http://www.13cflux.net/fluxml">
```

The parser must handle both forms equivalently. When parsing:

1. First detect namespace from root element
2. Register namespace prefix mapping
3. Use namespace-aware XPath queries

**Example parser setup:**

```python
from lxml import etree

tree = etree.parse(path)
root = tree.getroot()
nsmap = root.nsmap

# Handle default namespace
if None in nsmap:
    ns = nsmap[None]
    nsmap = {'fml': ns}
else:
    ns = next(iter(nsmap.values()))

# Use namespace in queries
pools = root.xpath('//fml:pool', namespaces=nsmap)
```

## Elements

### info (Optional)

Model metadata and documentation.

```xml
<info>
  <name>E. coli central carbon metabolism</name>
  <version>2.0</version>
  <date>2024-01-15 10:30:00</date>
  <modeler>John Doe</modeler>
  <strain>E. coli K-12 MG1655</strain>
  <comment>Model for glucose metabolism study</comment>
</info>
```

| Element | Description |
|---------|-------------|
| `name` | Model identifier |
| `version` | FluxML version number |
| `date` | Timestamp (format: `YYYY-MM-DD HH:MM:SS`) |
| `modeler` | Creator information |
| `strain` | Organism/strain specification |
| `comment` | Descriptive text |
| `signature` | Base64-encoded validation signature |

### reactionnetwork (Required)

Contains the metabolic network structure.

#### metabolitepools

```xml
<metabolitepools>
  <pool id="G6P" atoms="6"/>
  <pool id="Pyr" atoms="3" size="0.5"/>
  <pool id="AcCoA" atoms="2" cfg="C2N0">
    <annotation name="KEGG">C00024</annotation>
  </pool>
</metabolitepools>
```

| Attribute | Description | Default |
|-----------|-------------|---------|
| `id` | Unique identifier (required) | - |
| `atoms` | Number of traceable atoms (0-1024) | 0 |
| `size` | Pool size value | 1.0 |
| `cfg` | Atom configuration (e.g., `C6N0`) | - |

#### reaction

```xml
<reaction id="pgi" bidirectional="true">
  <reduct id="G6P" cfg="abcdef"/>
  <rproduct id="F6P" cfg="abcdef"/>
  <annotation name="EC">5.3.1.9</annotation>
</reaction>

<!-- Multiple substrates/products -->
<reaction id="ald" bidirectional="true">
  <reduct id="FBP" cfg="abcdef"/>
  <rproduct id="DHAP" cfg="cba"/>
  <rproduct id="GAP" cfg="def"/>
</reaction>

<!-- With variants for ambiguous mappings -->
<reaction id="tk1" bidirectional="true">
  <reduct id="X5P" cfg="abcde"/>
  <reduct id="R5P" cfg="fghij"/>
  <rproduct id="S7P" cfg="abfghij"/>
  <rproduct id="GAP" cfg="cde"/>
  <variant ratio="0.5">
    <!-- Alternative atom mapping -->
  </variant>
</reaction>
```

| Attribute | Description | Default |
|-----------|-------------|---------|
| `id` | Unique reaction identifier | - |
| `bidirectional` | Reversibility flag | `true` |

**Atom Mapping Notations:**

1. **Letter-based** (standard):
   ```xml
   <reduct id="G6P" cfg="abcdef"/>
   ```

2. **InChI-like** (multi-tracer, v3.0):
   ```xml
   <reduct id="G6P" cfg="C#1@1 C#2@1 C#3@1 C#4@1 C#5@1 C#6@1"/>
   ```
   Format: `AtomType#Position@MoleculeIndex`

### constraints (Optional)

Global stoichiometric constraints on fluxes and pool sizes.

```xml
<constraints>
  <net>
    <textual>v1 - v2 = 0</textual>
    <!-- Or MathML format -->
  </net>
  <xch>
    <textual>xch_pgi >= 0.1</textual>
  </xch>
  <psize>
    <textual>Pyr + OAA = 1</textual>
  </psize>
</constraints>
```

| Element | Description |
|---------|-------------|
| `<net>` | Net flux constraints |
| `<xch>` | Exchange flux constraints |
| `<psize>` | Pool size constraints |

### configuration (Multiple allowed)

Experiment-specific settings. Multiple configurations enable comparative studies.

```xml
<configuration name="glucose_labeling" stationary="true">
  <input>...</input>
  <measurement>...</measurement>
  <simulation>...</simulation>
</configuration>
```

| Attribute | Description | Default |
|-----------|-------------|---------|
| `name` | Unique configuration name | - |
| `stationary` | Steady-state flag | `true` |
| `time` | Time point (non-stationary) | - |

#### input

Substrate labeling specification.

```xml
<input pool="Gluc_ext" type="isotopomer">
  <label cfg="111111" purity="0.99">0.25</label>
  <label cfg="100000" purity="0.99">0.75</label>
</input>
```

| Attribute | Description |
|-----------|-------------|
| `pool` | Reference to metabolite pool |
| `type` | Input type: `isotopomer`, `cumomer`, `emu` |

**Label element:**

| Attribute | Description |
|-----------|-------------|
| `cfg` | Isotope pattern (`0`=unlabeled, `1`=labeled, `x`=unknown) |
| `purity` | Isotope purity (0-1) |
| `cost` | Label cost (for experimental design) |

#### measurement

Contains measurement model and experimental data.

```xml
<measurement>
  <mlabel>
    <date>2024-01-15</date>
    <version>1.0</version>
  </mlabel>
  <model>
    <labelingmeasurement>...</labelingmeasurement>
    <fluxmeasurement>...</fluxmeasurement>
    <poolsizemeasurement>...</poolsizemeasurement>
  </model>
  <data>
    <dlabel>...</dlabel>
    <datum>...</datum>
  </data>
</measurement>
```

**labelingmeasurement** (MS/NMR data):

```xml
<labelingmeasurement>
  <group id="ms_Ala" scale="auto">
    <textual>Ala[2,3]#M0,1,2</textual>
  </group>
  <group id="ms_Glu" scale="one" times="0,0.5,1,2">
    <textual>Glu[1-5]#M0,1,2,3,4,5</textual>
  </group>
</labelingmeasurement>
```

Format: `Metabolite[atoms]#M0,1,2,...`
- `[1,2,3]` or `[1-3]` - atom positions
- `#M0,1,2` - mass isotopomers to measure

**fluxmeasurement**:

```xml
<fluxmeasurement>
  <netflux id="v_upt">
    <textual>v_upt</textual>
  </netflux>
  <xchflux id="xch_pgi">
    <textual>pgi</textual>
  </xchflux>
</fluxmeasurement>
```

**poolsizemeasurement**:

```xml
<poolsizemeasurement>
  <poolsize id="ps_Pyr">
    <textual>Pyr</textual>
  </poolsize>
</poolsizemeasurement>
```

#### data

Experimental measurement values.

```xml
<data>
  <dlabel>
    <experiment>Exp001</experiment>
    <operator>Jane Doe</operator>
  </dlabel>

  <!-- Stationary MS data -->
  <datum id="ms_Ala" stddev="0.01" weight="0">0.35</datum>
  <datum id="ms_Ala" stddev="0.01" weight="1">0.42</datum>
  <datum id="ms_Ala" stddev="0.01" weight="2">0.23</datum>

  <!-- Non-stationary MS data -->
  <datum id="ms_Glu" stddev="0.01" weight="0" time="0">0.95</datum>
  <datum id="ms_Glu" stddev="0.01" weight="0" time="0.5">0.72</datum>
  <datum id="ms_Glu" stddev="0.01" weight="0" time="1">0.45</datum>

  <!-- Flux measurement -->
  <datum id="v_upt" stddev="0.05">1.0</datum>

  <!-- Pool size measurement -->
  <datum id="ps_Pyr" stddev="0.02">0.5</datum>
</data>
```

| Attribute | Description |
|-----------|-------------|
| `id` | Reference to measurement group |
| `stddev` | Standard deviation |
| `weight` | Mass isotopomer index (M+0, M+1, ...) |
| `time` | Time point (non-stationary) |
| `row` | Row index (1-256) |
| `pos` | Atom position |
| `type` | Measurement type: `S`, `DL`, `DR`, `DD`, `T` |

#### simulation

Simulation parameters and initial values.

```xml
<simulation type="auto" method="emu">
  <variables>
    <fluxvalue flux="v1" type="net" lo="0" hi="100" inc="1">10.5</fluxvalue>
    <fluxvalue flux="pgi" type="xch" lo="0" hi="1">0.5</fluxvalue>
    <poolsizevalue pool="Pyr" lo="0.01" hi="10">0.5</poolsizevalue>
  </variables>
</simulation>
```

| Attribute | Description |
|-----------|-------------|
| `type` | `auto`, `explicit`, `full` |
| `method` | `auto`, `cumomer`, `emu` |

## Measurement Formats

### MS Labeling Notation

```
Metabolite[Atoms]#MassIsotopomers
```

Examples:
- `Ala[2,3]#M0,1,2` - Alanine carbons 2,3; measure M+0, M+1, M+2
- `Glu[1-5]#M0,1,2,3,4,5` - Glutamate all 5 carbons; full MDV
- `Phe[1,2]#M0,1,2` - Phenylalanine fragment

### NMR Labeling Notation

Uses cumomer patterns with `0`, `1`, `x`:
- `0` - unlabeled
- `1` - labeled
- `x` - unknown/either

Example: `Ala:x1x` - measure carbon 2 labeling regardless of positions 1 and 3

## Data Files in Repository

| Path | Description |
|------|-------------|
| `data/convert/FluxML/` | Conversion test files |
| `data/tests/fluxml_models/` | Unit test models |

## Example File Structure

```
EC.fml
├── <info> - E. coli model metadata
├── <reactionnetwork>
│   ├── <metabolitepools> - 60+ metabolites
│   └── <reaction> - 70+ reactions
├── <constraints> - Mass balance constraints
└── <configuration name="exp1">
    ├── <input> - 13C glucose labeling
    ├── <measurement>
    │   ├── <labelingmeasurement> - MS groups
    │   ├── <fluxmeasurement> - Uptake flux
    │   └── <poolsizemeasurement>
    └── <data> - Experimental values
```

## Related Tools

- **libFluxML**: C/C++ library for parsing FluxML
- **FluxML Validator**: Online validation tool
- **13CFLUX2**: Simulation software using FluxML
