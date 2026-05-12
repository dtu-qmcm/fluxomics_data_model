"""FluxML metabolite (pool) definitions.

Each metabolite in a FluxML model is represented as a pool element with an
optional number of labelling positions (``atoms``), an optional chemical
formula, compartment, and annotations.

The :data:`ELEMENTS_AND_MOLECULAR_WEIGHTS` dictionary covers all 109 elements
from hydrogen (H) to ununhexium (Uuh) and is used by the ``molecular_weight``
computed property to derive the molecular weight from the chemical formula.

Corresponds to ``fluxml/reactionnetwork/metabolitepools/pool``.
"""

from typing import Optional, List
from pydantic import BaseModel, Field, field_validator, computed_field
import re
import jax.numpy as jnp
from ..core.common import Annotation, JAXArray


# Atomic weights of elements (atomic mass units)
ELEMENTS_AND_MOLECULAR_WEIGHTS = {
    "H": 1.007940,
    "He": 4.002602,
    "Li": 6.941000,
    "Be": 9.012182,
    "B": 10.811000,
    "C": 12.010700,
    "N": 14.006700,
    "O": 15.999400,
    "F": 18.998403,
    "Ne": 20.179700,
    "Na": 22.989770,
    "Mg": 24.305000,
    "Al": 26.981538,
    "Si": 28.085500,
    "P": 30.973761,
    "S": 32.065000,
    "Cl": 35.453000,
    "Ar": 39.948000,
    "K": 39.098300,
    "Ca": 40.078000,
    "Sc": 44.955910,
    "Ti": 47.867000,
    "V": 50.941500,
    "Cr": 51.996100,
    "Mn": 54.938049,
    "Fe": 55.845000,
    "Co": 58.933200,
    "Ni": 58.693400,
    "Cu": 63.546000,
    "Zn": 65.409000,
    "Ga": 69.723000,
    "Ge": 72.640000,
    "As": 74.921600,
    "Se": 78.960000,
    "Br": 79.904000,
    "Kr": 83.798000,
    "Rb": 85.467800,
    "Sr": 87.620000,
    "Y": 88.905850,
    "Zr": 91.224000,
    "Nb": 92.906380,
    "Mo": 95.940000,
    "Tc": 98.000000,
    "Ru": 101.070000,
    "Rh": 102.905500,
    "Pd": 106.420000,
    "Ag": 107.868200,
    "Cd": 112.411000,
    "In": 114.818000,
    "Sn": 118.710000,
    "Sb": 121.760000,
    "Te": 127.600000,
    "I": 126.904470,
    "Xe": 131.293000,
    "Cs": 132.905450,
    "Ba": 137.327000,
    "La": 138.905500,
    "Ce": 140.116000,
    "Pr": 140.907650,
    "Nd": 144.240000,
    "Pm": 145.000000,
    "Sm": 150.360000,
    "Eu": 151.964000,
    "Gd": 157.250000,
    "Tb": 158.925340,
    "Dy": 162.500000,
    "Ho": 164.930320,
    "Er": 167.259000,
    "Tm": 168.934210,
    "Yb": 173.040000,
    "Lu": 174.967000,
    "Hf": 178.490000,
    "Ta": 180.947900,
    "W": 183.840000,
    "Re": 186.207000,
    "Os": 190.230000,
    "Ir": 192.217000,
    "Pt": 195.078000,
    "Au": 196.966550,
    "Hg": 200.590000,
    "Tl": 204.383300,
    "Pb": 207.200000,
    "Bi": 208.980380,
    "Po": 209.000000,
    "At": 210.000000,
    "Rn": 222.000000,
    "Fr": 223.000000,
    "Ra": 226.000000,
    "Ac": 227.000000,
    "Th": 232.038100,
    "Pa": 231.035880,
    "U": 238.028910,
    "Np": 237.000000,
    "Pu": 244.000000,
    "Am": 243.000000,
    "Cm": 247.000000,
    "Bk": 247.000000,
    "Cf": 251.000000,
    "Es": 252.000000,
    "Fm": 257.000000,
    "Md": 258.000000,
    "No": 259.000000,
    "Lr": 262.000000,
    "Rf": 261.000000,
    "Db": 262.000000,
    "Sg": 266.000000,
    "Bh": 264.000000,
    "Hs": 277.000000,
    "Mt": 268.000000,
    "Ds": 281.000000,
    "Rg": 272.000000,
    "Cn": 285.000000,
    "Uuq": 289.000000,
    "Uuh": 292.000000,
}


class Metabolite(BaseModel):
    """FluxML metabolite definition.

    Corresponds to fluxml/reactionnetwork/metabolites/metabolite
    """

    id: str = Field(description="Metabolite identifier")
    name: Optional[str] = Field(default=None, description="Metabolite name")
    atoms: Optional[int] = Field(
        default=0,
        ge=0,
        le=1024,
        description="Number of atoms for the labelling experiment",
    )
    weight: Optional[float] = Field(
        default=None, ge=0.0, le=50000.0, description="Molecular weight (g/mol)"
    )
    charge: Optional[int] = Field(
        default=None, ge=-100, le=100, description="Charge of the metabolite"
    )
    formula: Optional[str] = Field(
        default=None,
        description="Chemical elements invovled in labelling experiment",
    )
    compartment: Optional[str] = Field(
        default=None, description="Compartment where metabolite is located"
    )
    annotations: List[Annotation] = Field(
        default_factory=list, description="Annotations"
    )
    involved_in_variants: bool = Field(
        default=False,
        description="Flag indicating if this metabolite participates "
        "in variant reactions",
    )

    # JAX-compatible numerical representation
    jax_atoms: Optional[JAXArray] = Field(
        default=None, exclude=True, description="JAX atom array"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("formula")
    @classmethod
    def validate_formula(cls, v: Optional[str]) -> Optional[str]:
        """Validate chemical formula string."""
        if not v or v == "0" or v == "":
            return v
        pattern = r"^(([A-Z][a-z]?)([0-9.]+[0-9.]?|(?=[A-Z])?))+$"
        if not re.match(pattern, v):
            raise ValueError(f"Invalid chemical formula format: {v}")
        return v

    @staticmethod
    def _parse_formula(formula: str) -> dict[str, float]:
        """Parse chemical formula into element counts."""
        if not formula or formula in ["0", ""]:
            return {}

        elements = {}
        # Pattern to match element-count pairs
        pattern = r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)"
        matches = re.findall(pattern, formula)

        for element, count_str in matches:
            if count_str == "":
                count = 1.0
            else:
                count = float(count_str)
            elements[element] = elements.get(element, 0.0) + count

        return elements

    @computed_field
    @property
    def molecular_weight(self) -> Optional[float]:
        """Calculate molecular weight from formula if available."""
        if not self.formula or self.formula in ["0", ""]:
            return None

        try:
            elements = self._parse_formula(self.formula)
            total_weight = 0.0

            for element, count in elements.items():
                if element not in ELEMENTS_AND_MOLECULAR_WEIGHTS:
                    return None  # Unknown element, can't calculate
                total_weight += ELEMENTS_AND_MOLECULAR_WEIGHTS[element] * count

            return round(total_weight, 6)
        except Exception:
            return None

    @property
    def atom_vector(self) -> jnp.ndarray:
        """Get JAX array representation of atoms."""
        if self.jax_atoms is None:
            # Create atom vector based on configuration
            if self.atoms == 0:
                return jnp.array([])
            return jnp.ones(self.atoms, dtype=jnp.float32)
        return self.jax_atoms.to_jax_array()

    def with_jax_atoms(self, atoms: jnp.ndarray) -> "Metabolite":
        """Create new metabolite with JAX atom array."""
        return Metabolite(
            id=self.id,
            name=self.name,
            atoms=self.atoms,
            weight=self.weight,
            charge=self.charge,
            formula=self.formula,
            compartment=self.compartment,
            annotations=self.annotations,
            jax_atoms=JAXArray.from_jax_array(atoms),
        )
