# Format Compatibility Analysis

This document provides a comprehensive analysis of compatibility issues between FluxML, FreeFlux, influx_si, and INCA formats, with implementation strategies for each issue.

---

## Executive Summary

| Issue Category | Count | Impact |
|----------------|-------|--------|
| **Blocking (Cannot Convert)** | 3 | Multi-tracer, flux ratios, error models |
| **Lossy (Information Loss)** | 8 | Experiment organization, metadata, grouping |
| **Representational (Transformable)** | 9 | Atom mapping, symmetry, external metabolites |

**Key Finding**: FluxML is the most feature-rich format. Converting FROM FluxML to others will always lose information. Converting TO FluxML preserves all information.

---

## 1. BLOCKING Incompatibilities

These features exist in one format but CANNOT be represented in others. Conversion must either **fail** or **drop the feature with a warning**.

### 1.1 Multi-Tracer Atom Mapping (FluxML v3.0 Only)

| Aspect | Details |
|--------|---------|
| **FluxML** | `cfg="C#1@1 C#2@1 N#1@2"` - Tracks 13C, 15N, 2H, 18O, 34S simultaneously |
| **FreeFlux** | Carbon only - `(abcdef)` |
| **influx_si** | Carbon only - `(ABCDEF)` |
| **INCA** | Carbon only - `(abcdef)` |

**Problem**: Multi-tracer models track multiple isotopes through the same reaction network. For example, tracking both 13C and 15N through amino acid biosynthesis. This information structure cannot be represented in single-tracer formats.

**Impact**:
- Nitrogen flux in amino acid synthesis pathways
- Hydrogen/deuterium exchange studies
- 18O incorporation in oxidation reactions

**Implementation Strategy**:

```python
class MultiTracerConversionError(ConversionError):
    """Raised when multi-tracer model cannot be converted."""
    pass

def convert_fluxml_to_freeflux(model: FluxomicsModel) -> None:
    if model.has_multi_tracer_mapping():
        raise MultiTracerConversionError(
            "Model contains multi-tracer atom mappings (13C + 15N/2H/18O). "
            "FreeFlux only supports 13C. Options:\n"
            "  1. Use --extract-carbon to export 13C-only subset\n"
            "  2. Convert to FluxML format instead"
        )

# Option: Extract carbon-only subset
def extract_carbon_mapping(mapping: ReactionAtomMapping) -> ReactionAtomMapping:
    """Extract only 13C transitions from multi-tracer mapping."""
    carbon_transitions = [t for t in mapping.transitions if t.atom_type == "C"]
    return ReactionAtomMapping(
        transitions=carbon_transitions,
        notation_type="lowercase",
        _warning="Multi-tracer mapping reduced to 13C only"
    )
```

### 1.2 Flux Ratio Measurements (FluxML Only)

| Aspect | Details |
|--------|---------|
| **FluxML** | `<fluxratios>` element for direct ratio fitting |
| **Others** | Not supported |

**Problem**: FluxML can fit ratios like `v1/v2 = 0.3 ± 0.05` directly. Other formats require deriving this from individual flux estimates.

**Example**:
```xml
<!-- FluxML -->
<fluxratios>
  <ratio id="ppp_split">
    <textual>zwf / (pgi + zwf)</textual>
    <value>0.3</value>
    <stddev>0.05</stddev>
  </ratio>
</fluxratios>
```

**Implementation Strategy**:

```python
def convert_flux_ratios(ratios: List[FluxRatio], target_format: str) -> ConversionResult:
    if target_format != "fluxml":
        warnings.warn(
            DataLossWarning(
                f"Flux ratio measurements dropped: {[r.id for r in ratios]}. "
                f"{target_format} does not support direct ratio fitting."
            )
        )
        return ConversionResult(converted=[], dropped=ratios)
```

### 1.3 Custom Error Models (FluxML Only)

| Aspect | Details |
|--------|---------|
| **FluxML** | Custom uncertainty formulas based on measurement value |
| **Others** | Simple standard deviation only |

**Problem**: FluxML can express heteroscedastic errors (e.g., `stddev = 0.01 * value + 0.005`) that scale with measurement magnitude.

**Implementation Strategy**:

```python
def convert_error_model(error_model: ErrorModel, target_format: str) -> float:
    if error_model.is_custom():
        # Evaluate at representative value to get approximate stddev
        representative_stddev = error_model.evaluate(value=0.5)
        warnings.warn(
            DataLossWarning(
                f"Custom error model '{error_model.formula}' converted to "
                f"constant stddev={representative_stddev:.4f}"
            )
        )
        return representative_stddev
    return error_model.stddev
```

---

## 2. LOSSY Conversions

These features can be partially represented but with information loss. Conversion succeeds with warnings.

### 2.1 Multiple Experiments in One File

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Structure** | Multiple `<configuration>` | One experiment per directory | One experiment per file set | Multiple via `m.expts` array |
| **Max experiments** | Unlimited | 1 | 1 (or .vmtf batch) | Unlimited |

**Problem**: FluxML stores parallel labeling experiments (e.g., U-13C vs 1-13C glucose) in one file. FreeFlux/influx_si require separate directories.

**Example**:
```xml
<!-- FluxML: Two experiments in one file -->
<fluxml>
  <reactionnetwork>...</reactionnetwork>  <!-- Shared network -->
  <configuration name="U13C_glucose">...</configuration>
  <configuration name="1_13C_glucose">...</configuration>
</fluxml>
```

**Implementation Strategy**:

```python
def convert_multi_experiment_fluxml(model: FluxomicsModel, output_dir: Path) -> List[Path]:
    """Split multi-experiment FluxML into multiple output files."""
    output_paths = []

    for exp_name, experiment in model.experiments.items():
        # Create experiment-specific directory
        exp_dir = output_dir / sanitize_filename(exp_name)
        exp_dir.mkdir(exist_ok=True)

        # Create single-experiment model
        single_exp_model = model.copy_with_experiment(exp_name)

        # Write to format-specific files
        writer.write(single_exp_model, exp_dir)
        output_paths.append(exp_dir)

    return output_paths
```

**Output structure for FreeFlux**:
```
output/
├── U13C_glucose/
│   ├── reactions.tsv        # Same network
│   ├── labeling_input.tsv   # U-13C specific
│   └── measured_MDVs.tsv    # Experiment-specific data
└── 1_13C_glucose/
    ├── reactions.tsv        # Same network (duplicate)
    ├── labeling_input.tsv   # 1-13C specific
    └── measured_MDVs.tsv    # Experiment-specific data
```

### 2.2 Measurement Grouping with Scaling

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Groups** | `<group scale="auto">` | No | Dataset column | No |
| **Systematic error** | Group-level scale factors | No | No | No |

**Problem**: FluxML can group measurements (e.g., all MS fragments from same GC-MS run) and estimate a common scaling factor for systematic error correction.

**Implementation Strategy**:

```python
def flatten_measurement_groups(groups: List[MeasurementGroup]) -> List[Measurement]:
    """Flatten grouped measurements to individual measurements."""
    measurements = []
    for group in groups:
        if group.scale_type != "one":
            warnings.warn(
                DataLossWarning(
                    f"Group '{group.id}' scale='{group.scale_type}' converted to "
                    f"individual measurements (systematic error handling lost)"
                )
            )
        for measurement in group.measurements:
            measurements.append(measurement.copy(group_id=None))
    return measurements
```

### 2.3 Rich Metadata

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Model name** | Yes | Partial (filename) | Partial (filename) | Comments |
| **Version** | Yes | No | No | No |
| **Date** | Yes | No | No | No |
| **Modeler** | Yes | No | No | No |
| **Strain** | Yes | No | No | No |
| **Signature** | Yes (validation) | No | No | No |

**Implementation Strategy**:

```python
def convert_metadata(metadata: ModelMetadata, target_format: str) -> dict:
    """Convert metadata, preserving what's possible in target format."""
    lost_fields = []

    result = {}
    if target_format == "freeflux":
        # Only filename available
        if metadata.name:
            result["filename_base"] = sanitize_filename(metadata.name)
        lost_fields = ["version", "date", "modeler", "strain", "signature"]

    if lost_fields:
        warnings.warn(
            DataLossWarning(
                f"Metadata fields lost in {target_format}: {lost_fields}"
            )
        )

    return result
```

### 2.4 Pool Size Ratio Measurements (FluxML Only)

Similar to flux ratios - ratio measurements like `[Pyr]/[PEP] = 2.3 ± 0.2` cannot be represented.

### 2.5 Time-Varying Substrate Profiles

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Constant substrate** | Yes | Yes | Yes | Yes |
| **Step change** | Yes | Limited | Partial | Yes |
| **Continuous profile** | Yes | No | No | Partial |

**Problem**: FluxML supports arbitrary substrate labeling profiles over time. Others assume constant labeling input.

### 2.6 NMR Measurement Notation

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Notation** | Cumomer `x1x` | Not documented | `2->1,3` doublet | Not documented |
| **Support** | Full | Limited | Full | Limited |

**Problem**: Different NMR notations for coupling patterns. Some formats lack documentation.

### 2.7 Constraint Expressiveness

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Linear equality** | Yes (textual/MathML) | No | Yes (.cnstr) | Implicit (bounds) |
| **Linear inequality** | Yes | No | Yes | Implicit |
| **Flux bounds** | Yes (lo/hi) | Code-only (`set_flux_bounds`) | No | Yes (lb/ub) |
| **Measured flux constraint** | Yes | Code-only (`set_measured_flux`) | Yes (.mflux) | Yes |
| **Non-linear** | Via MathML | No | No | No |

**Problem**: FreeFlux has **no file-based constraints** and **no linear equality/inequality constraints**. It only supports:
- Flux bounds via `set_flux_bounds()` Python API
- Soft equality via `set_measured_flux()` (measured values with SD)

**Note**: FreeFlux's `set_concentration_bounds()` does NOT constrain optimization - it only sets initial guess sampling range.

### 2.8 Variable Classification (influx_si Exclusive)

| Aspect | FluxML | FreeFlux | influx_si | INCA |
|--------|--------|----------|-----------|------|
| **Free/Dependent/Constrained** | Not explicit | Not explicit | `.tvar` file | Not explicit |

**Problem**: influx_si explicitly classifies variables as Free (to estimate), Dependent (calculated), or Constrained (fixed). This classification is lost when converting to other formats.

---

## 3. REPRESENTATIONAL Differences

These are equivalent representations that require transformation but no information loss.

### 3.1 Atom Mapping Notation

| Format | Case | Example | Multi-substrate |
|--------|------|---------|-----------------|
| FluxML | lowercase | `cfg="abcdef"` | Separate `<reduct>` elements |
| FreeFlux | lowercase | `OAA(abcd)+AcCoA(ef)` | `+` concatenation |
| influx_si | UPPERCASE | `OAA (ABCD) + AcCoA (EF)` | `+` concatenation |
| INCA | lowercase | `OAA (abcd) + AcCoA (ef)` | `+` concatenation |

**Transformation**:

```python
def convert_atom_mapping_notation(
    mapping: str,
    source_format: str,
    target_format: str
) -> str:
    """Convert between atom mapping notations."""
    if source_format == "influx_si":
        # UPPERCASE -> lowercase
        mapping = mapping.lower()
    elif target_format == "influx_si":
        # lowercase -> UPPERCASE
        mapping = mapping.upper()
    return mapping
```

### 3.2 Symmetric Metabolite Handling

| Format | Representation | Example |
|--------|---------------|---------|
| FluxML | `<variant ratio="0.5">` | Reaction variants |
| FreeFlux | Comma notation | `Suc(abcd,dcba)` |
| influx_si | Split reactions + constraint | `fum_a - fum_b == 0` |
| INCA | `.sym` property | `m.mets{'Suc'}.sym = atommap('1:4 2:3...')` |

**Transformation Logic**:

```python
class SymmetryConverter:
    """Convert symmetric metabolite representation between formats."""

    def to_freeflux(self, metabolite_id: str, symmetry: SymmetryGroup) -> str:
        """Generate comma notation: Suc(abcd,dcba)"""
        variants = []
        for transform in symmetry.transformations:
            variant = self._apply_transform(transform)
            variants.append(variant)
        return f"{metabolite_id}({','.join(variants)})"

    def to_influx_si(
        self,
        reaction: Reaction,
        symmetry: SymmetryGroup
    ) -> Tuple[List[Reaction], LinearConstraint]:
        """Generate split reactions + equality constraint."""
        reactions = []
        for i, transform in enumerate(symmetry.transformations):
            new_reaction = reaction.copy()
            new_reaction.id = f"{reaction.id}_{chr(97+i)}"  # fum_a, fum_b
            new_reaction.atom_mapping = self._apply_transform(
                reaction.atom_mapping, transform
            )
            reactions.append(new_reaction)

        # Add constraint: fum_a - fum_b == 0
        constraint = LinearConstraint(
            coefficients={r.id: (1 if i == 0 else -1)
                         for i, r in enumerate(reactions)},
            operator=ConstraintOperator.EQ,
            value=0.0,
            kind="net"
        )

        return reactions, constraint

    def to_inca(self, metabolite_id: str, symmetry: SymmetryGroup) -> str:
        """Generate INCA .sym notation."""
        mappings = []
        for transform in symmetry.transformations[0].items():
            mappings.append(f"{transform[0]}:{transform[1]}")
        return f"m.mets{{'{metabolite_id}'}}.sym = list('{symmetry.name}', atommap('{' '.join(mappings)}'));"
```

### 3.3 External Metabolite Conventions

| Format | Convention | Example |
|--------|------------|---------|
| FluxML | Pool attribute or naming | `<pool id="Gluc_ext" external="true"/>` |
| FreeFlux | Implicit (substrate inputs) | Not explicit in network |
| influx_si | `_ext` suffix | `Gluc_ext` |
| INCA | `.ext` or `.ex` suffix | `Gluc.ext` |

**Transformation**:

```python
def normalize_external_metabolite_id(
    metabolite_id: str,
    source_format: str,
    target_format: str,
    is_external: bool
) -> str:
    """Convert external metabolite naming conventions."""
    # Remove source-specific suffix
    base_id = metabolite_id
    if source_format == "influx_si" and base_id.endswith("_ext"):
        base_id = base_id[:-4]
    elif source_format == "inca" and (base_id.endswith(".ext") or base_id.endswith(".ex")):
        base_id = base_id.rsplit(".", 1)[0]

    # Add target-specific suffix
    if is_external:
        if target_format == "influx_si":
            return f"{base_id}_ext"
        elif target_format == "inca":
            return f"{base_id}.ext"

    return base_id
```

### 3.4 Labeling Input Representation

| Format | Notation | Example |
|--------|----------|---------|
| FluxML | Binary pattern + purity | `<label cfg="111111" purity="0.99">0.25</label>` |
| FreeFlux | Binary in API | `labeling_patterns=['111111']` |
| influx_si | Binary in column | `Isotopomer: 111111` |
| INCA | Position list | `@ 1 2 3 4 5 6` |

**Transformation**:

```python
def convert_labeling_notation(
    pattern: str,
    source_format: str,
    target_format: str
) -> str:
    """Convert between labeling notations."""
    # Normalize to binary
    if source_format == "inca":
        # "@ 1 2 3" -> "111000" (for 6-carbon)
        positions = [int(p) for p in pattern.replace("@", "").split()]
        n_atoms = max(positions)  # Need actual atom count from metabolite
        binary = ''.join('1' if i in positions else '0'
                        for i in range(1, n_atoms + 1))
    else:
        binary = pattern

    # Convert to target
    if target_format == "inca":
        # "111000" -> "@ 1 2 3"
        positions = [str(i + 1) for i, c in enumerate(binary) if c == '1']
        return "@ " + " ".join(positions) if positions else "@"

    return binary
```

### 3.5 Reversibility Representation

| Format | Reversible | Irreversible | Special |
|--------|------------|--------------|---------|
| FluxML | `bidirectional="true"` | `bidirectional="false"` | - |
| FreeFlux | `reversibility=1` | `reversibility=0` | - |
| influx_si | `<->` or `<->>` | `->` | `->>` (non-negative net) |
| INCA | `<->` | `->` | - |

**Note**: influx_si's `->>` (non-negative net flux) has no direct equivalent in other formats. Convert to constraint:

```python
def convert_arrow_type(arrow: str, target_format: str) -> Tuple[bool, Optional[LinearConstraint]]:
    """Convert influx_si arrow types."""
    if arrow == "->>":  # Irreversible, non-negative net
        if target_format == "fluxml":
            # Add constraint: net_flux >= 0
            return False, LinearConstraint(
                coefficients={"reaction_id": 1},
                operator=ConstraintOperator.GE,
                value=0.0
            )
        return False, None
    elif arrow == "<->>":  # Reversible, non-negative net
        return True, LinearConstraint(
            coefficients={"reaction_id": 1},
            operator=ConstraintOperator.GE,
            value=0.0
        )
    elif arrow == "<->":
        return True, None
    else:  # "->"
        return False, None
```

### 3.6 Flux Value Representation

| Format | Representation |
|--------|---------------|
| FluxML | Net + exchange OR forward + backward |
| FreeFlux | `v1` (net) OR `v1_f`/`v1_b` (fwd/bwd) |
| influx_si | NET and XCH columns |
| INCA | Array values |

**Transformation**:

```python
def convert_flux_values(
    net: float,
    exchange: float,
    target_format: str
) -> dict:
    """Convert between flux representations."""
    forward = exchange + max(0, net)
    backward = exchange + max(0, -net)

    if target_format == "freeflux":
        if backward > 0:
            return {"_f": forward, "_b": backward}
        return {"": net}
    elif target_format == "influx_si":
        return {"NET": net, "XCH": exchange}

    return {"net": net, "exchange": exchange}
```

### 3.7 MS Measurement Fragment Notation

| Format | Notation | Example |
|--------|----------|---------|
| FluxML | `Metabolite[atoms]#M0,1,2` | `Ala[2,3]#M0,1,2` |
| FreeFlux | `Metabolite_positions` | `Ala_23` |
| influx_si | Columns | `Specie: Ala, Fragment: 2,3` |
| INCA | `name: metabolite @ positions` | `Ala2: Ala @ 2 3` |

**Transformation**:

```python
class MSFragmentNotation:
    """Convert between MS fragment notations."""

    @staticmethod
    def parse(notation: str, format: str) -> Tuple[str, List[int]]:
        """Parse fragment notation to (metabolite, positions)."""
        if format == "fluxml":
            # Ala[2,3]#M0,1,2 -> Ala, [2, 3]
            match = re.match(r'(\w+)\[([^\]]+)\]', notation)
            metabolite = match.group(1)
            positions = [int(p) for p in match.group(2).replace('-', ',').split(',')]
        elif format == "freeflux":
            # Ala_23 -> Ala, [2, 3]
            metabolite, pos_str = notation.rsplit('_', 1)
            positions = [int(p) for p in pos_str]
        # ... etc

        return metabolite, positions

    @staticmethod
    def format(metabolite: str, positions: List[int], format: str) -> str:
        """Format fragment to target notation."""
        if format == "fluxml":
            return f"{metabolite}[{','.join(map(str, positions))}]"
        elif format == "freeflux":
            return f"{metabolite}_{''.join(map(str, positions))}"
        elif format == "influx_si":
            return f"{','.join(map(str, positions))}"  # Fragment column only
        elif format == "inca":
            return f"{metabolite}{len(positions)}: {metabolite} @ {' '.join(map(str, positions))}"
```

### 3.8 Biomass Reaction Representation

| Format | Representation |
|--------|---------------|
| FluxML | Regular reaction with stoichiometry |
| FreeFlux | Implicit or simple reaction |
| influx_si | `bs_*` prefix reactions + constraints |
| INCA | Stoichiometry in reaction string |

**Transformation**:

```python
def convert_biomass_reaction(
    biomass: BiomassComposition,
    target_format: str
) -> Union[Reaction, Tuple[List[Reaction], List[LinearConstraint]]]:
    """Convert biomass representation."""

    if target_format == "influx_si":
        # Create individual bs_* reactions + constraints
        reactions = []
        constraints = []

        for precursor_id, coefficient in biomass.precursors.items():
            rxn_id = f"bs_{precursor_id.lower()}"
            reactions.append(Reaction(
                id=rxn_id,
                reactants=[ReactionParticipant(precursor_id, 1.0)],
                products=[ReactionParticipant("Biomass", 1.0)],
                reaction_type=ReactionType.BIOMASS
            ))

            # Constraint: bs_ala - 0.488*BM == 0
            constraints.append(LinearConstraint(
                coefficients={rxn_id: 1.0, "BM": -coefficient},
                operator=ConstraintOperator.EQ,
                value=0.0
            ))

        return reactions, constraints

    else:
        # Create single reaction with stoichiometry
        reactants = [
            ReactionParticipant(met_id, coef)
            for met_id, coef in biomass.precursors.items()
        ]
        return Reaction(
            id="biomass",
            reactants=reactants,
            products=[ReactionParticipant("Biomass", 1.0)],
            reaction_type=ReactionType.BIOMASS
        )
```

### 3.9 Constraint Syntax

| Format | Syntax | Example |
|--------|--------|---------|
| FluxML | Textual or MathML | `v1 - v2 = 0` |
| influx_si | Column-based | `Formula: v1-v2, Operator: ==, Value: 0` |
| FreeFlux | Not supported | - |
| INCA | Implicit in bounds | - |

**Transformation**:

```python
def parse_constraint_formula(formula: str) -> Dict[str, float]:
    """Parse constraint formula to coefficients."""
    # "v1 + 0.5*v2 - v3" -> {"v1": 1.0, "v2": 0.5, "v3": -1.0}
    coefficients = {}
    # Parse using tokenizer
    return coefficients

def format_constraint(
    constraint: LinearConstraint,
    target_format: str
) -> str:
    """Format constraint for target format."""
    if target_format == "fluxml":
        terms = [f"{coef}*{var}" if coef != 1 else var
                for var, coef in constraint.coefficients.items()]
        return f"{' + '.join(terms)} {constraint.operator.value} {constraint.value}"
    # ... etc
```

---

## 4. Conversion Strategy Matrix

### 4.1 Feature Support Matrix

| Feature | FluxML | FreeFlux | influx_si | INCA |
|---------|:------:|:--------:|:---------:|:----:|
| **Basic Features** |
| Reactions with atom mapping | ✓ | ✓ | ✓ | ✓ |
| Reversible reactions | ✓ | ✓ | ✓ | ✓ |
| External metabolites | ✓ | ✓ | ✓ | ✓ |
| MS measurements | ✓ | ✓ | ✓ | ✓ |
| Flux measurements | ✓ | ✓ | ✓ | ✓ |
| Pool size measurements | ✓ | ✓ | ✓ | ✓ |
| **Advanced Features** |
| Symmetric metabolites | ✓ | ✓ | ✓ | ✓ |
| Non-stationary experiments | ✓ | ✓ | ✓ | ✓ |
| Linear constraints | ✓ | ✗ | ✓ | ~ |
| Multiple experiments | ✓ | ✗ | ✗ | ✓ |
| **FluxML-Exclusive** |
| Multi-tracer (C+N+H) | ✓ | ✗ | ✗ | ✗ |
| Flux ratio measurements | ✓ | ✗ | ✗ | ✗ |
| Custom error models | ✓ | ✗ | ✗ | ✗ |
| Measurement grouping | ✓ | ✗ | ~ | ✗ |
| Pool size ratios | ✓ | ✗ | ✗ | ✗ |

Legend: ✓ = Full support, ~ = Partial support, ✗ = Not supported

### 4.2 Conversion Loss Matrix

This matrix shows what is LOST when converting FROM row TO column:

| FROM → TO | FluxML | FreeFlux | influx_si | INCA |
|-----------|:------:|:--------:|:---------:|:----:|
| **FluxML** | - | High | Medium | Medium |
| **FreeFlux** | None | - | Low | Low |
| **influx_si** | None | Low | - | Low |
| **INCA** | None | Low | Low | - |

**Details**:

- **FluxML → FreeFlux** (High loss): Multi-tracer, flux ratios, error models, constraints, multi-experiment, measurement groups, pool size ratios
- **FluxML → influx_si** (Medium loss): Multi-tracer, flux ratios, error models, multi-experiment, pool size ratios
- **FluxML → INCA** (Medium loss): Multi-tracer, flux ratios, error models, measurement groups, pool size ratios
- **FreeFlux → influx_si** (Low loss): Symmetric notation converted; labeling input restructured
- **influx_si → FreeFlux** (Low loss): Constraints dropped; variable classification lost

### 4.3 Recommended Conversion Paths

```
                    ┌─────────────────┐
                    │     FluxML      │
                    │ (Most Complete) │
                    └────────┬────────┘
                             │
         ┌───────────────────┼───────────────────┐
         │                   │                   │
         ▼                   ▼                   ▼
    ┌─────────┐        ┌──────────┐        ┌─────────┐
    │ FreeFlux │ ◄────► │ influx_si│ ◄────► │  INCA   │
    └─────────┘        └──────────┘        └─────────┘
         │                   │                   │
         └───────────────────┼───────────────────┘
                             │
                    ┌────────▼────────┐
                    │   JAX Export    │
                    │ (Computation)   │
                    └─────────────────┘
```

**Recommendations**:

1. **For archival**: Always save in FluxML format (lossless)
2. **For computation**: Convert to tool-specific format as needed
3. **For round-trip**: FluxML → Other → FluxML may lose information; track original FluxML

---

## 5. Implementation Strategies

### 5.1 Strict vs. Permissive Mode

```python
class ConversionMode(Enum):
    STRICT = "strict"       # Fail on any data loss
    PERMISSIVE = "permissive"  # Warn but continue
    EXTRACT = "extract"     # Extract compatible subset

def convert(
    model: FluxomicsModel,
    target_format: str,
    mode: ConversionMode = ConversionMode.PERMISSIVE
) -> ConversionResult:

    issues = []

    # Check for blocking issues
    if model.has_multi_tracer_mapping():
        if mode == ConversionMode.STRICT:
            raise MultiTracerConversionError("...")
        elif mode == ConversionMode.EXTRACT:
            model = extract_carbon_only(model)
            issues.append(DataLossWarning("Extracted 13C-only mapping"))
        else:
            raise MultiTracerConversionError("...")

    # Check for lossy features
    if model.has_flux_ratios():
        if mode == ConversionMode.STRICT:
            raise ConversionError("Flux ratios cannot be converted")
        issues.append(DataLossWarning("Flux ratios dropped"))

    return ConversionResult(model=converted_model, warnings=issues)
```

### 5.2 Conversion Report

Generate detailed report of what was converted, modified, or lost:

```python
@dataclass
class ConversionReport:
    source_format: str
    target_format: str
    timestamp: str

    # Statistics
    reactions_converted: int
    metabolites_converted: int
    measurements_converted: int

    # Modifications
    notation_changes: List[str]      # e.g., "Atom mapping: uppercase -> lowercase"
    structure_changes: List[str]     # e.g., "Symmetric reactions split"

    # Losses
    features_dropped: List[str]      # e.g., "Flux ratios"
    data_simplified: List[str]       # e.g., "Error model -> constant stddev"

    # Warnings
    warnings: List[str]

    def to_markdown(self) -> str:
        """Generate human-readable report."""
        ...
```

### 5.3 Round-Trip Fidelity

Track original format information for better round-trip conversion:

```python
@dataclass
class OriginalFormatInfo:
    """Preserved information for round-trip conversion."""
    format: str
    original_notation: Dict[str, str]  # {element_id: original_string}
    original_ids: Dict[str, str]       # {canonical_id: original_id}
    format_specific: Dict[str, Any]    # Format-specific metadata

@dataclass
class FluxomicsModel:
    # ... existing fields ...
    _original_format_info: Optional[OriginalFormatInfo] = None

    def can_round_trip_to(self, format: str) -> bool:
        """Check if model can round-trip to original format."""
        if self._original_format_info is None:
            return False
        return self._original_format_info.format == format
```

---

## 6. Test Cases for Compatibility

### 6.1 Atom Mapping Round-Trip

```python
def test_atom_mapping_round_trip():
    """Test atom mapping conversion between all formats."""
    original = "OAA(abcd)+AcCoA(ef) -> Cit(dcbfea)"

    # Parse as FreeFlux
    mapping = FreeFluxParser.parse_reaction(original)

    # Convert to each format and back
    for format in ["fluxml", "influx_si", "inca"]:
        converted = mapping.to_format(format)
        back = parse_from_format(converted, format)
        assert back.to_freeflux() == original
```

### 6.2 Symmetric Metabolite Conversion

```python
def test_symmetric_conversion_freeflux_to_influx():
    """FreeFlux comma notation -> influx_si split reactions."""
    freeflux_reaction = "v1: Suc(abcd,dcba) -> Fum(abcd,dcba)"

    model = FreeFluxReader.read_reaction(freeflux_reaction)
    influx_model = convert_to_influx_si(model)

    # Should have two reactions
    assert "v1_a" in influx_model.reactions
    assert "v1_b" in influx_model.reactions

    # Should have constraint
    assert any(
        c.formula == "v1_a - v1_b" and c.value == 0
        for c in influx_model.constraints.net_flux
    )
```

### 6.3 Multi-Experiment Split

```python
def test_multi_experiment_split():
    """FluxML multi-experiment -> multiple FreeFlux directories."""
    fluxml_model = FluxMLReader.read("multi_exp.fml")
    assert len(fluxml_model.experiments) == 3

    output_dirs = convert_to_freeflux(fluxml_model, Path("output/"))

    assert len(output_dirs) == 3
    for exp_name in fluxml_model.experiments:
        assert (Path("output") / exp_name / "reactions.tsv").exists()
```

---

## 7. Summary of Implementation Requirements

### 7.1 Core Converters Needed

| Converter | Complexity | Key Challenges |
|-----------|------------|----------------|
| Atom Mapping | Medium | Case conversion, InChI parsing |
| Symmetric Metabolites | High | Four different representations |
| Labeling Input | Medium | Position ↔ binary conversion |
| MS Fragments | Medium | Four different notations |
| Constraints | High | FreeFlux lacks support |
| Biomass Reactions | High | Detect and normalize |
| External Metabolites | Low | Suffix conventions |
| Flux Values | Low | Net/exchange ↔ fwd/bwd |

---

## 8. Flux Variable Representation Comparison

This section compares how flux variables are classified and represented across formats.

### 8.1 influx_si .tvar Format (Most Explicit)

influx_si is the ONLY format with explicit variable classification:

| Column | Description |
|--------|-------------|
| `Name` | Reaction ID |
| `Kind` | `NET` (net flux), `XCH` (exchange flux), `METAB` (pool size) |
| `Type` | `F` (free/to estimate), `D` (dependent/calculated), `C` (constrained/fixed) |
| `Value` | Initial value (required for F and C types) |

**Example .tvar:**
```
Name    Kind    Type    Value
v_upt   NET     F       1.0       # Free parameter - will be estimated
pgi     NET     D                 # Dependent - calculated from constraints
pgi     XCH     F       0.75      # Free exchange flux
pfk     NET     C       0.8       # Constrained - fixed value
Pyr     METAB   F       0.5       # Pool size - free parameter
```

### 8.2 FluxML Flux Representation

FluxML uses `<simulation>` element with bounds and initial values:

```xml
<simulation type="auto" method="emu">
  <variables>
    <!-- Net flux with bounds -->
    <fluxvalue flux="pgi" type="net" lo="0" hi="100" inc="1" ed="0">10.5</fluxvalue>

    <!-- Exchange flux -->
    <fluxvalue flux="pgi" type="xch" lo="0" hi="1">0.5</fluxvalue>

    <!-- Pool size -->
    <poolsizevalue pool="Pyr" lo="0.01" hi="10">0.5</poolsizevalue>
  </variables>
</simulation>
```

**Key Differences from influx_si:**
- No explicit Free/Dependent/Constrained classification
- Uses bounds (`lo`, `hi`) instead of type
- `ed="0"` means "do not estimate" (equivalent to `C`)
- Dependent fluxes are implicit from stoichiometry

### 8.3 FreeFlux Flux Representation

FreeFlux uses `fluxes.tsv` for file-based values, but **bounds are set in Python code**:

**File-based (fluxes.tsv):**
```
#flux_ID    value
v1          10
v6_f        12.5    # Forward flux
v6_b        7.5     # Backward flux
```

**Code-based bounds (Python API):**
```python
# Set bounds for all fluxes
fit.set_flux_bounds('all', bounds=[-100, 100])

# Set bounds for specific flux
fit.set_flux_bounds('v1', bounds=[0, 50])

# Measured flux (acts as equality constraint)
fit.set_measured_flux('v1', mean=10, sd=1)
```

**Key Differences:**
- **NO file-based constraints** - bounds only via Python API
- **NO linear equality/inequality constraints** (like `v1 - v2 = 0`)
- **NO explicit variable classification** (Free/Dependent/Constrained)
- Uses `_f`/`_b` suffix for forward/backward instead of NET/XCH
- `set_measured_flux()` acts as soft equality constraint via fitting

### 8.4 INCA Flux Representation

INCA uses MATLAB arrays:

```matlab
% Bounds (apply to all)
m.rates.flx.lb = 1e-7;    % Lower bound
m.rates.flx.ub = 1e7;     % Upper bound

% Values (indexed by reaction order)
m.rates.flx.val = [10.5, 8.2, ...];
```

**Key Differences:**
- **Array-based** - values ordered by reaction position
- **Global bounds** - same bounds for all fluxes by default
- **NO explicit variable classification**
- Individual bounds can be set via `m.rates.flx{i}.lb`

### 8.5 Flux Variable Conversion Matrix

| Feature | influx_si | FluxML | FreeFlux | INCA |
|---------|:---------:|:------:|:--------:|:----:|
| Net flux values | ✓ | ✓ | ✓ (implicit) | ✓ |
| Exchange flux values | ✓ | ✓ | ✓ (`_f`/`_b`) | ~ |
| Free/Dependent/Constrained | ✓ | ~ (via bounds) | ✗ | ✗ |
| Bounds (lo, hi) | ✗ | ✓ (file) | ✓ (code only) | ✓ |
| Bounds in file | ✗ | ✓ | ✗ | ✓ |
| Per-flux uncertainty | ✗ | ~ | ✗ | ✗ |
| Linear constraints | ✓ (.cnstr) | ✓ (textual) | ✗ | ~ |

### 8.6 Conversion Strategies

**influx_si → Others:**
```python
def convert_tvar_to_fluxml(tvar_data):
    """Convert .tvar classification to FluxML bounds."""
    for var in tvar_data:
        if var.type == 'C':  # Constrained
            yield FluxValue(
                flux=var.name,
                type=var.kind.lower(),
                value=var.value,
                lo=var.value,  # Fixed: lo = hi = value
                hi=var.value,
                ed=0  # Do not estimate
            )
        elif var.type == 'F':  # Free
            yield FluxValue(
                flux=var.name,
                type=var.kind.lower(),
                value=var.value,
                lo=0,
                hi=1000,  # Default range
                ed=1  # Estimate
            )
        # Dependent (D) - omit from explicit list
```

**Others → influx_si:**
```python
def convert_fluxml_to_tvar(simulation):
    """Convert FluxML to .tvar - LOSSY for dependent classification."""
    warnings.warn("Dependent flux classification cannot be determined")

    for fluxval in simulation.fluxvalues:
        if fluxval.lo == fluxval.hi:  # Fixed value
            yield TVar(name=fluxval.flux, kind=fluxval.type.upper(),
                       type='C', value=fluxval.value)
        else:  # Free parameter
            yield TVar(name=fluxval.flux, kind=fluxval.type.upper(),
                       type='F', value=fluxval.value)
    # NOTE: Dependent variables must be inferred from stoichiometry
```

---

## 9. Metabolite Concentration Representation Comparison

### 9.1 influx_si .mmet Format

Measured metabolite concentrations with uncertainty:

```
Id    Comment    Specie    Value    SD
                 Pyr       0.426    0.05
                 PEP       0.182    0.02
                 AKG       0.41e-5  0.01e-5
```

**Features:**
- Measurements with standard deviation (for fitting)
- Sparse: only measured metabolites listed
- Units typically mM or mmol/gDW

### 9.2 FluxML Pool Size Representation

Two locations for pool sizes:

**1. Initial values in `<metabolitepools>`:**
```xml
<pool id="Pyr" atoms="3" size="0.5"/>
```

**2. Measurements in `<poolsizemeasurement>`:**
```xml
<poolsizemeasurement>
  <poolsize id="ps_Pyr">
    <textual>Pyr</textual>
  </poolsize>
</poolsizemeasurement>

<data>
  <datum id="ps_Pyr" stddev="0.05">0.426</datum>
</data>
```

**3. Estimation parameters in `<simulation>`:**
```xml
<poolsizevalue pool="Pyr" lo="0.01" hi="10">0.5</poolsizevalue>
```

### 9.3 FreeFlux Concentration Representation

**File-based (concentrations.tsv):**
```
#metab_ID    value
OAA          0.1
Cit          5
Pyr          0.5
```

**Code-based bounds (Python API):**
```python
# Set concentration bounds - NOTE: only for initial guess sampling!
fit.set_concentration_bounds('all', bounds=[0.001, 10])

# Per-metabolite bounds
fit.set_concentration_bounds('Pyr', bounds=[0.1, 1.0])
```

**IMPORTANT:** From FreeFlux docs: "Bounds specified by `set_concentration_bounds` are not used as constraints in the optimization, it just defines the sampling range of the initial guess of concentrations."

**Limitations:**
- **NO uncertainty/SD in files** - just values
- **NO real optimization bounds** - only initial guess sampling
- No measured vs initial distinction in files
- Units: μmol gCDW⁻¹

### 9.4 INCA Pool Size Representation

Array-based with global bounds:

```matlab
% Bounds
m.states.lb = 1e-6;
m.states.ub = 1;

% Values (ordered by metabolite index)
m.states.val = [0.5, 0.3, 1.6e-3, ...];
```

### 9.5 Pool Size Conversion Matrix

| Feature | influx_si | FluxML | FreeFlux | INCA |
|---------|:---------:|:------:|:--------:|:----:|
| Pool values | ✓ | ✓ | ✓ | ✓ |
| Uncertainty (SD) | ✓ | ✓ | ✗ | ✗ |
| Optimization bounds | ✗ | ✓ | ✗* | ✓ |
| Measured vs Initial | ✓ (separate files) | ✓ (separate elements) | ✗ | ~ |
| Per-metabolite | ✓ | ✓ | ✓ | Array-indexed |
| File-based | ✓ | ✓ | ✓ | ✓ |
| Code-based bounds | N/A | N/A | ✓ (initial guess only) | N/A |

*FreeFlux's `set_concentration_bounds()` only defines initial guess sampling, NOT optimization constraints

---

## 10. Stationary vs Non-Stationary Experiment Comparison

This is a critical difference affecting data structure and file organization.

### 10.1 Conceptual Difference

| Aspect | Stationary (Isotopic Steady-State) | Non-Stationary (INST/Instationary) |
|--------|-----------------------------------|-----------------------------------|
| **Assumption** | Labeling patterns constant | Labeling patterns change over time |
| **Data** | Single MDV per fragment | Time-course MDVs |
| **Pool sizes** | Not needed | Required for fitting |
| **Information** | Relative fluxes | Absolute fluxes + pool sizes |
| **Experiment** | Hours to reach steady-state | Minutes of labeling dynamics |

### 10.2 FluxML Non-Stationary Support

Most explicit and feature-rich:

```xml
<!-- Configuration-level flag -->
<configuration name="inst_exp" stationary="false">

  <!-- Measurement model with time points -->
  <labelingmeasurement>
    <group id="ms_Glu" scale="one" times="0,0.5,1,2,5">
      <textual>Glu[1-5]#M0,1,2,3,4,5</textual>
    </group>
  </labelingmeasurement>

  <!-- Time-indexed data -->
  <data>
    <datum id="ms_Glu" stddev="0.01" weight="0" time="0">0.95</datum>
    <datum id="ms_Glu" stddev="0.01" weight="0" time="0.5">0.72</datum>
    <datum id="ms_Glu" stddev="0.01" weight="0" time="1">0.45</datum>
    <datum id="ms_Glu" stddev="0.01" weight="0" time="2">0.32</datum>
  </data>

  <!-- Pool size estimation -->
  <simulation>
    <variables>
      <poolsizevalue pool="Glu" lo="0.01" hi="10">0.5</poolsizevalue>
    </variables>
  </simulation>
</configuration>
```

**Key Features:**
- `stationary="false"` attribute
- `times="0,0.5,1,2,5"` in measurement group
- `time="X"` attribute on datum elements

### 10.3 influx_si Non-Stationary Support

Uses Time column in .miso file:

```
# Stationary measurement (empty Time)
Id    Comment    Specie    Fragment    Dataset    Isospecies    Value    SD    Time
                 Ala       2,3         MS-1       M0            0.35     0.01

# Non-stationary measurement (Time specified)
                 AKG       1-5         inst       M0            0.95     0.01   0
                 AKG       1-5         inst       M0            0.72     0.01   0.5
                 AKG       1-5         inst       M0            0.45     0.01   1.0
                 AKG       1-5         inst       M0            0.32     0.01   2.0
```

**Features:**
- Same file format, Time column populated
- Dataset often named "inst" for instationary
- Pool sizes in .mmet file (required)
- .tvar includes METAB kind for pool size variables

### 10.4 FreeFlux Non-Stationary Support

Uses separate file `measured_inst_MDVs.tsv`:

```
# Stationary file: measured_MDVs.tsv
#fragment_ID    mean                              sd
Glu_12345       0.328,0.276,0.274,0.088,0.03,0.004    0.01,0.01,...

# Non-stationary file: measured_inst_MDVs.tsv
#fragment_ID    time    mean                     sd
Glu_12345       0       0.948,0.051,0.001,...    0.01,0.01,...
Glu_12345       0.1     0.928,0.06,0.012,...     0.01,0.01,...
Glu_12345       0.2     0.873,0.085,0.041,...    0.01,0.01,...
Glu_12345       1       0.487,0.233,0.212,...    0.01,0.01,...
```

**Features:**
- Separate file types for stationary vs non-stationary
- Time column added
- `concentrations.tsv` required for pool sizes

### 10.5 INCA Non-Stationary Support

Uses simulation options:

```matlab
% Simulation mode
m.options.sim_ss = false;  % false = non-stationary

% Time points
m.options.int_tspan = [0 0.1 0.2 0.5 1 2];

% Pool size initial values
m.states.val = [0.5, 0.3, 1.6e-3, ...];
m.states.lb = 1e-6;
m.states.ub = 1;
```

### 10.6 Non-Stationary Conversion Matrix

| Feature | FluxML | influx_si | FreeFlux | INCA |
|---------|:------:|:---------:|:--------:|:----:|
| Stationary/Non-stationary flag | ✓ (`stationary=`) | ~ (Time column) | ~ (file type) | ✓ (`sim_ss`) |
| Time points in model | ✓ (`times=`) | ✗ | ✗ | ✓ (`int_tspan`) |
| Time in data | ✓ (`time=`) | ✓ (Time column) | ✓ (time column) | Array index |
| Pool size measurements | ✓ | ✓ (.mmet) | ~ (concentrations.tsv) | ✓ |
| Pool size bounds | ✓ | ✗ | ✗ | ✓ |

---

## 11. NMR Measurement Representation Comparison

NMR measurements track specific isotopomer patterns rather than mass distributions.

### 11.1 NMR Terminology

| Term | Description |
|------|-------------|
| **Singlet** | No 13C-13C coupling (adjacent carbons unlabeled) |
| **Doublet (D)** | Coupling to one adjacent 13C |
| **Doublet of Doublets (DD)** | Coupling to two adjacent 13Cs |
| **Triplet (T)** | Coupling to two adjacent 13Cs (special geometry) |
| **Cumomer** | Cumulative isotopomer (uses 0, 1, x notation) |

### 11.2 FluxML NMR Notation

Uses cumomer pattern notation:

```xml
<labelingmeasurement>
  <!-- Cumomer notation: 0=unlabeled, 1=labeled, x=either -->
  <group id="nmr_Ala_C2" scale="one">
    <textual>Ala:x1x</textual>  <!-- C2 labeling, any C1/C3 -->
  </group>

  <!-- Multiple patterns -->
  <group id="nmr_Glu" scale="one">
    <textual>Glu:xx1xx</textual>  <!-- C3 labeled -->
  </group>
</labelingmeasurement>

<data>
  <!-- Multiplet fractions -->
  <datum id="nmr_Ala_C2" type="S">0.35</datum>   <!-- Singlet -->
  <datum id="nmr_Ala_C2" type="DL">0.28</datum>  <!-- Doublet Left -->
  <datum id="nmr_Ala_C2" type="DR">0.22</datum>  <!-- Doublet Right -->
  <datum id="nmr_Ala_C2" type="DD">0.15</datum>  <!-- Doublet of Doublets -->
</data>
```

**Notation Elements:**
- Pattern: `Metabolite:cumomer_pattern`
- Cumomer: `x` = any, `0` = unlabeled, `1` = labeled
- Types: `S`, `DL`, `DR`, `DD`, `T`

### 11.3 influx_si NMR Notation

Uses arrow notation for coupling patterns:

```
Id    Comment    Specie    Fragment    Dataset    Isospecies    Value    SD    Time
                 Ala                   HSQC-C2    2->           0.35     0.02       # Singlet
                 Ala                   HSQC-C2    2->1          0.28     0.02       # Doublet coupled to C1
                 Ala                   HSQC-C2    2->3          0.22     0.02       # Doublet coupled to C3
                 Ala                   HSQC-C2    2->1,3        0.15     0.02       # DD coupled to C1 and C3
```

**Notation Elements:**
- `2->` : Observe C2, singlet (no coupling)
- `2->1` : Observe C2, doublet coupled to C1
- `2->3` : Observe C2, doublet coupled to C3
- `2->1,3` : Observe C2, doublet of doublets coupled to C1 and C3

### 11.4 NMR Notation Comparison

| Pattern | FluxML Cumomer | influx_si Arrow |
|---------|---------------|-----------------|
| C2 singlet (Ala) | `Ala:x1x` + type=S | `2->` |
| C2 doublet from C1 | `Ala:11x` + type=DL | `2->1` |
| C2 doublet from C3 | `Ala:x11` + type=DR | `2->3` |
| C2 doublet of doublets | `Ala:111` + type=DD | `2->1,3` |
| C3 triplet (5-carbon) | `Met:xx1xx` + type=T | `3->2,4` |

### 11.5 FreeFlux and INCA NMR Support

| Format | NMR Support |
|--------|-------------|
| **FreeFlux** | Limited - not documented in specification |
| **INCA** | Limited - primarily MS-focused |

### 11.6 NMR Conversion Strategy

```python
class NMRNotationConverter:
    """Convert between FluxML cumomer and influx_si arrow notations."""

    @staticmethod
    def fluxml_to_influx_si(metabolite: str, cumomer: str, multiplet_type: str) -> str:
        """
        Convert FluxML notation to influx_si notation.

        Example: Ala:x1x + DL -> 2->1
        """
        # Find observed position (the '1' in pattern)
        observed_pos = cumomer.index('1') + 1

        if multiplet_type == 'S':
            return f"{observed_pos}->"
        elif multiplet_type == 'DL':
            coupled_pos = observed_pos - 1
            return f"{observed_pos}->{coupled_pos}"
        elif multiplet_type == 'DR':
            coupled_pos = observed_pos + 1
            return f"{observed_pos}->{coupled_pos}"
        elif multiplet_type == 'DD':
            return f"{observed_pos}->{observed_pos-1},{observed_pos+1}"
        elif multiplet_type == 'T':
            return f"{observed_pos}->{observed_pos-1},{observed_pos+1}"

    @staticmethod
    def influx_si_to_fluxml(arrow_notation: str, n_atoms: int) -> Tuple[str, str]:
        """
        Convert influx_si notation to FluxML cumomer + type.

        Example: 2->1 -> x1x + DL
        """
        if '->' not in arrow_notation:
            raise ValueError(f"Invalid NMR notation: {arrow_notation}")

        observed, coupled = arrow_notation.split('->')
        observed_pos = int(observed)

        # Build cumomer pattern
        pattern = ['x'] * n_atoms
        pattern[observed_pos - 1] = '1'

        if not coupled:  # Singlet
            return ''.join(pattern), 'S'

        coupled_positions = [int(p) for p in coupled.split(',')]

        if len(coupled_positions) == 1:
            if coupled_positions[0] < observed_pos:
                return ''.join(pattern), 'DL'
            else:
                return ''.join(pattern), 'DR'
        else:
            return ''.join(pattern), 'DD'
```

### 11.7 NMR Conversion Limitations

| Conversion | Issue |
|------------|-------|
| FluxML → influx_si | Generally lossless |
| influx_si → FluxML | Need metabolite atom count |
| Either → FreeFlux | **BLOCKING** - limited NMR support |
| Either → INCA | **BLOCKING** - limited NMR support |

---

## 12. Symmetric Metabolite Architecture Recommendation

### 12.1 Current Problem

Symmetric metabolites are tightly coupled to reactions in the specifications:

| Format | Coupling |
|--------|----------|
| **FreeFlux** | Comma notation in reaction: `Suc(abcd,dcba)` |
| **influx_si** | Separate reactions: `fum_a`, `fum_b` + constraint |
| **FluxML** | `<variant>` child of `<reaction>` |
| **INCA** | Metabolite property: `m.mets{'Suc'}.sym` |

### 12.2 Recommended Architecture: Separate Classes

**Rationale:** Separating symmetric handling from reactions allows:
1. Format-independent symmetry representation
2. Cleaner conversion logic
3. Metabolite-centric symmetry (like INCA)
4. Proper handling when same metabolite appears in multiple reactions

**Proposed Class Structure:**

```python
@dataclass
class SymmetryDefinition:
    """
    Metabolite symmetry definition, INDEPENDENT of reactions.

    This should be defined once per metabolite, not per reaction.
    """
    metabolite_id: str
    name: str  # e.g., "C2_rotation", "plane_symmetry"

    # Position mapping: original_pos -> equivalent_pos
    # For succinate C2 symmetry: {1: 4, 2: 3, 3: 2, 4: 1}
    position_mapping: Dict[int, int]

    # Probability of each variant (default 0.5 for 2-fold symmetry)
    probability: float = 0.5

    def apply(self, atom_sequence: str) -> str:
        """Apply symmetry transformation to atom sequence."""
        result = list(atom_sequence)
        for orig, new in self.position_mapping.items():
            result[new - 1] = atom_sequence[orig - 1]
        return ''.join(result)


@dataclass
class AtomMapping:
    """
    Atom mapping SEPARATE from symmetry.

    Stores the canonical (non-scrambled) mapping.
    Symmetry variants are generated on-demand from SymmetryDefinition.
    """
    transitions: List[AtomTransition]

    # Reference to metabolites with symmetry (not embedded)
    symmetric_metabolites: List[str] = field(default_factory=list)

    def get_all_variants(self, symmetry_defs: Dict[str, SymmetryDefinition]) -> List['AtomMapping']:
        """Generate all symmetry variants using external symmetry definitions."""
        if not self.symmetric_metabolites:
            return [self]

        variants = [self]
        for met_id in self.symmetric_metabolites:
            if met_id in symmetry_defs:
                sym_def = symmetry_defs[met_id]
                new_variants = []
                for variant in variants:
                    # Apply symmetry transformation
                    transformed = variant.apply_symmetry(met_id, sym_def)
                    new_variants.extend([variant, transformed])
                variants = new_variants

        return variants


@dataclass
class Reaction:
    """
    Reaction WITHOUT embedded symmetry information.
    """
    id: str
    reactants: List[ReactionParticipant]
    products: List[ReactionParticipant]

    # Canonical atom mapping (no variants)
    atom_mapping: Optional[AtomMapping] = None

    # Just references to symmetric metabolites involved
    # Actual symmetry defined at metabolite level
    involves_symmetric_metabolites: List[str] = field(default_factory=list)


@dataclass
class Metabolite:
    """
    Metabolite WITH symmetry as a property.
    """
    id: str
    atoms: int

    # Symmetry definition (if symmetric)
    symmetry: Optional[SymmetryDefinition] = None

    @property
    def is_symmetric(self) -> bool:
        return self.symmetry is not None
```

### 12.3 Conversion with Separated Architecture

**FreeFlux comma notation → Universal:**
```python
def parse_freeflux_symmetric(reaction_str: str, model: FluxomicsModel):
    """
    Parse: Suc(abcd,dcba) -> Fum(abcd,dcba)

    Extracts symmetry to metabolite level, not reaction level.
    """
    # Parse reactants/products with variants
    for participant, variants in parse_participants_with_variants(reaction_str):
        if len(variants) > 1:
            # Create symmetry definition on METABOLITE
            sym_def = infer_symmetry_from_variants(variants)

            # Add to metabolite (not reaction)
            metabolite = model.get_metabolite(participant.id)
            if metabolite.symmetry is None:
                metabolite.symmetry = sym_def

            # Mark reaction as involving symmetric metabolite
            reaction.involves_symmetric_metabolites.append(participant.id)
```

**influx_si split reactions → Universal:**
```python
def parse_influx_si_symmetric(reactions: List[str], constraints: List[str], model: FluxomicsModel):
    """
    Parse: fum_a: Mal(ABCD) <-> Fum(ABCD)
           fum_b: Mal(ABCD) <-> Fum(DCBA)
           constraint: fum_a - fum_b == 0

    Merges into single reaction + metabolite symmetry.
    """
    # Detect symmetric reaction pairs from constraints
    symmetric_pairs = find_symmetric_pairs(reactions, constraints)

    for rxn_a, rxn_b in symmetric_pairs:
        # Infer symmetry from mapping difference
        sym_def = infer_symmetry_from_reactions(rxn_a, rxn_b)

        # Apply to metabolite
        symmetric_met = find_symmetric_metabolite(rxn_a, rxn_b)
        model.get_metabolite(symmetric_met).symmetry = sym_def

        # Create single canonical reaction
        canonical = merge_symmetric_reactions(rxn_a, rxn_b)
        canonical.involves_symmetric_metabolites.append(symmetric_met)
        model.reactions.append(canonical)
```

### 12.4 Benefits of Separation

| Benefit | Description |
|---------|-------------|
| **Single source of truth** | Symmetry defined once per metabolite |
| **Cleaner conversion** | Don't need to propagate symmetry through reactions |
| **INCA compatibility** | Matches INCA's metabolite-centric approach |
| **Validation** | Easy to check all reactions involving symmetric metabolite |
| **Flexibility** | Can add/remove symmetry without touching reactions |

### 12.5 Atom Mapping Separation Recommendation

Similarly, atom mapping should be a separate entity:

```python
@dataclass
class AtomMappingLibrary:
    """
    Central repository of atom mappings.

    Allows same mapping to be reused across:
    - Different representations (letter, InChI)
    - Different symmetric variants
    - Round-trip conversions
    """

    # Mapping ID -> canonical AtomMapping
    mappings: Dict[str, AtomMapping] = field(default_factory=dict)

    # Reaction ID -> Mapping ID (reference, not copy)
    reaction_mappings: Dict[str, str] = field(default_factory=dict)

    def get_mapping_for_reaction(self, reaction_id: str) -> Optional[AtomMapping]:
        if reaction_id in self.reaction_mappings:
            mapping_id = self.reaction_mappings[reaction_id]
            return self.mappings.get(mapping_id)
        return None
```

---

### 7.2 Priority Order

1. **Critical Path**: Atom mapping → Reactions → Metabolites
2. **Measurements**: MS fragments → Labeling input → Flux/pool measurements
3. **Advanced**: Symmetric metabolites → Constraints → Biomass
4. **Polish**: Multi-experiment handling → Error models → Metadata

### 7.3 Required Error Types

```python
class ConversionError(Exception): pass
class MultiTracerConversionError(ConversionError): pass
class ConstraintConversionError(ConversionError): pass
class SymmetryConversionError(ConversionError): pass
class NotationParseError(ConversionError): pass
class DataLossWarning(UserWarning): pass
class FeatureDroppedWarning(DataLossWarning): pass
class NotationSimplifiedWarning(DataLossWarning): pass
```
