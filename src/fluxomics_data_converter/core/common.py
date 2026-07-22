"""Shared building-block data structures used across all FluxML modules.

This module provides:

- :class:`Annotation`       — key/value metadata tag (``<annotation>``)
- :class:`TextualOrMath`    — a value that is either a plain string or MathML
- :class:`ErrorModel`       — measurement error formula wrapper
- :class:`JAXArray`         — Pydantic-serialisable wrapper for a JAX array
- :class:`TimeSeries`       — time-indexed arrays for non-stationary data
- :class:`AtomTransitionsNetwork` — ``dict`` subclass with a rich ``__repr__``
- :class:`DictList`         — combined list + O(1)-by-id dict collection

All Pydantic models in this module use ``frozen=True`` (immutable after
construction).
"""

from typing import Optional, Union, Any, List, TypeVar, Generic
from pydantic import BaseModel, Field, GetCoreSchemaHandler
from pydantic_core import core_schema
import jax.numpy as jnp
from itertools import islice


class Annotation(BaseModel):
    """FluxML annotation element for additional metadata.

    Corresponds to fluxml/reactionnetwork/metabolites/metabolite/annotation
    and fluxml/reactionnetwork/reaction/annotation.
    """

    name: str = Field(description="Annotation name")
    content: Optional[str] = Field(
        default=None, description="Annotation content"
    )

    class Config:
        frozen = True


class TextualOrMath(BaseModel):
    """FluxML textual or MathML content.

    Used for constraints and mathematical expressions.
    """

    textual: Optional[str] = Field(
        default=None, description="Textual representation"
    )
    mathml: Optional[str] = Field(
        default=None, description="MathML representation"
    )

    class Config:
        frozen = True

    def __init__(self, **data):
        """Initialise, requiring at least one of textual or mathml.

        Raises:
            ValueError: If neither textual nor mathml is provided.
        """
        super().__init__(**data)
        if not self.textual and not self.mathml:
            raise ValueError("Either textual or mathml must be provided")


class ErrorModel(BaseModel):
    """FluxML error model for measurement uncertainties.

    Corresponds to fluxml/experiments/measurement/model/*/errormodel
    """

    expression: TextualOrMath = Field(description="Error model expression")

    class Config:
        frozen = True


class JAXArray(BaseModel):
    """Pydantic-serialisable wrapper for a JAX array.

    JAX arrays cannot be directly stored inside a Pydantic model because
    they are not JSON-serialisable.  This class bridges the gap by
    storing the array's shape, dtype, and data as plain Python objects.
    Use :meth:`from_jax_array` to wrap a JAX array and :meth:`to_jax_array`
    to unwrap it.

    Example::

        arr = jnp.ones((3, 4), dtype=jnp.float32)
        wrapped = JAXArray.from_jax_array(arr)
        recovered = wrapped.to_jax_array()  # shape (3, 4), dtype float32
    """

    shape: tuple[int, ...] = Field(description="Array shape")
    dtype: str = Field(description="Array dtype")
    data: list = Field(description="Array data as nested lists")

    class Config:
        frozen = True

    @classmethod
    def from_jax_array(cls, arr: jnp.ndarray) -> "JAXArray":
        """Create from JAX array.

        Returns:
            A new ``JAXArray`` wrapping the given JAX array.
        """
        return cls(shape=arr.shape, dtype=str(arr.dtype), data=arr.tolist())

    def to_jax_array(self) -> jnp.ndarray:
        """Convert to JAX array.

        Returns:
            The reconstructed JAX array.
        """
        return jnp.array(self.data, dtype=self.dtype).reshape(self.shape)


class TimeSeries(BaseModel):
    """Time series data structure for non-stationary measurements."""

    times: "JAXArray" = Field(description="Time points")
    values: "JAXArray" = Field(description="Measurement values")
    errors: Optional["JAXArray"] = Field(
        default=None, description="Measurement errors"
    )

    class Config:
        frozen = True

    @classmethod
    def from_arrays(
        cls,
        times: jnp.ndarray,
        values: jnp.ndarray,
        errors: Optional[jnp.ndarray] = None,
    ) -> "TimeSeries":
        """Create from JAX arrays.

        Returns:
            A new ``TimeSeries`` built from the given JAX arrays.
        """
        return cls(
            times=JAXArray.from_jax_array(times),
            values=JAXArray.from_jax_array(values),
            errors=(
                JAXArray.from_jax_array(errors) if errors is not None else None
            ),
        )

    def to_arrays(
        self,
    ) -> tuple[jnp.ndarray, jnp.ndarray, Optional[jnp.ndarray]]:
        """Convert to JAX arrays.

        Returns:
            A tuple of ``(times, values, errors)`` JAX arrays, where
            ``errors`` may be ``None``.
        """
        return (
            self.times.to_jax_array(),
            self.values.to_jax_array(),
            self.errors.to_jax_array() if self.errors else None,
        )


class AtomTransitionsNetwork(dict):
    """Dictionary for atom transitions for a given metabolic network.

    Keys are reaction IDs, values are AtomTransition objects.
    """

    def __repr__(self) -> str:
        """Return string representation with summary information."""
        if not self:
            return "=== Atom Transitions ===\n  No atom transitions defined"

        lines = ["=== Atom Transitions ==="]

        # Total count
        lines.append(f"  Total: {len(self)}")

        # Count mappings with variants
        with_variants = 0
        for am in self.values():
            if hasattr(am, "maps") and len(am.maps) > 1:
                with_variants += 1

        if with_variants > 0:
            lines.append(f"  Reactions with multiple variants: {with_variants}")

        # Sample atom transitions
        sample_size = min(5, len(self))
        lines.append("")
        lines.append("  Sample atom transitions:")
        for rxn_id, am in list(self.items())[:sample_size]:
            if hasattr(am, "maps"):
                n_maps = len(am.maps)
                if n_maps > 1:
                    lines.append(f"    - {rxn_id}: {n_maps} variants")
                else:
                    # Show letter notation for single map
                    try:
                        notation = am.to_letter_notation()
                        if notation and len(notation) < 50:
                            lines.append(f"    - {rxn_id}: {notation}")
                        else:
                            lines.append(f"    - {rxn_id}: (mapping defined)")
                    except Exception:
                        lines.append(f"    - {rxn_id}: (mapping defined)")
            else:
                lines.append(f"    - {rxn_id}: (mapping defined)")

        if len(self) > sample_size:
            lines.append(f"    ... and {len(self) - sample_size} more")

        return "\n".join(lines)

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> core_schema.CoreSchema:
        """Get Pydantic core schema for serialization/validation.

        Returns:
            A Pydantic core schema that validates as a dict and returns
            an ``AtomTransitionsNetwork``.
        """
        # Get the schema for a dict
        dict_schema = handler.generate_schema(dict)

        # Return a schema that validates as a dict but returns an
        # AtomTransitionsNetwork
        return core_schema.no_info_after_validator_function(
            lambda v: cls(v),
            dict_schema,
        )


# Type variable for DictList items — must have an .id attribute.
T = TypeVar("T", bound=BaseModel)


class DictList(list, Generic[T]):
    """A combined list and dictionary providing O(1) lookup by ``id``.

    Behaves exactly like a ``list`` — supports indexing, slicing, iteration,
    ``len``, ``append``, ``extend``, etc. — but additionally maintains an
    internal ``{id: index}`` mapping so that individual items can be
    retrieved by their ``id`` attribute in O(1) time.

    Every item stored in a ``DictList`` must have an ``id`` attribute whose
    value is a unique string.  Attempting to insert a duplicate ID raises
    ``ValueError``.

    Example::

        metabolites = DictList([
            Metabolite(id="Glc", atoms=6),
            Metabolite(id="G6P", atoms=6),
        ])
        metabolites.get_by_id("G6P")   # O(1) lookup
        metabolites["G6P"]             # same, via __getitem__
        metabolites[0]                 # still works by integer index
        "Glc" in metabolites           # O(1) membership test
        metabolites.ids                # ["Glc", "G6P"] in insertion order

    The ``ids`` property preserves insertion order and is the correct source
    of metabolite/reaction orderings for numerical matrices (unlike iterating
    over a ``frozenset``).
    """

    def __init__(self, items: Optional[Union[List[T], "DictList[T]"]] = None):
        """Initialize a DictList.

        Parameters
        ----------
        items : list or DictList, optional
            Initial items to populate the DictList
        """
        super().__init__()
        self._dict: dict[str, int] = {}

        if items is not None:
            if isinstance(items, DictList):
                list.extend(self, items)
                self._dict = items._dict.copy()
            else:
                self.extend(items)

    def _generate_index(self) -> None:
        """Rebuild the _dict index from current list items."""
        self._dict = {item.id: idx for idx, item in enumerate(self)}

    def _check_id(self, item_id: str) -> None:
        """Check if an ID already exists in the DictList.

        Raises:
            ValueError: If the ID is already present in the list.
        """
        if item_id in self._dict:
            raise ValueError(f"ID '{item_id}' is already present in the list")

    def append(self, item: T) -> None:
        """Append an item to the end of the list."""
        item_id = item.id
        self._check_id(item_id)
        self._dict[item_id] = len(self)
        list.append(self, item)

    def extend(self, items: Union[List[T], "DictList[T]"]) -> None:
        """Extend list by appending elements from the iterable."""
        current_length = len(self)

        # First check all IDs are unique
        for item in items:
            self._check_id(item.id)

        # Then extend
        list.extend(self, items)

        # Update index
        for idx, item in enumerate(
            islice(self, current_length, None), current_length
        ):
            self._dict[item.id] = idx

    def insert(self, index: int, item: T) -> None:
        """Insert item before index."""
        self._check_id(item.id)
        list.insert(self, index, item)

        # Update indices for all items after insertion point
        for i, j in list(self._dict.items()):
            if j >= index:
                self._dict[i] = j + 1
        self._dict[item.id] = index

    def remove(self, item: Union[str, T]) -> None:
        """Remove first occurrence of item."""
        if isinstance(item, str):
            # Remove by ID
            index = self._dict[item]
        else:
            # Remove by object
            index = self.index(item)

        self.pop(index)

    def pop(self, index: int = -1) -> T:
        """Remove and return item at index (default last).

        Returns:
            The item removed from the list.
        """
        item = list.pop(self, index)

        # Remove from dict
        del self._dict[item.id]

        # Update indices if not popping from end
        if index != -1 and index < len(self):
            for key, idx in list(self._dict.items()):
                if idx > index:
                    self._dict[key] = idx - 1

        return item

    def clear(self) -> None:
        """Remove all items from the list."""
        list.clear(self)
        self._dict.clear()

    def index(self, item: Union[str, T], start: int = 0, stop: int = -1) -> int:
        """Return index of item in the list.

        Returns:
            The index of the item in the list.

        Raises:
            ValueError: If the item (looked up by ID) is not in the list.
        """
        if isinstance(item, str):
            # Look up by ID
            if item not in self._dict:
                raise ValueError(f"'{item}' is not in list")
            return self._dict[item]
        else:
            # Look up by object
            return list.index(
                self, item, start, stop if stop != -1 else len(self)
            )

    def get_by_id(self, item_id: str) -> T:
        """Get item by its ID attribute.

        Returns:
            The item with the given ID.

        Raises:
            KeyError: If no item with the given ID exists.
        """
        if item_id not in self._dict:
            raise KeyError(f"No item with ID '{item_id}'")
        return self[self._dict[item_id]]

    def get_by_any(self, key: Union[str, int, T]) -> T:
        """Get item by ID, index, or the item itself.

        Returns:
            The item matching the given key.

        Raises:
            ValueError: If the item cannot be found in the DictList.
        """
        if isinstance(key, int):
            return self[key]
        elif isinstance(key, str):
            return self.get_by_id(key)
        elif key in self:
            return key
        else:
            raise ValueError(f"Item {key} not found in DictList")

    def has_id(self, item_id: str) -> bool:
        """Check if an ID exists in the DictList.

        Returns:
            True if the ID exists, False otherwise.
        """
        return item_id in self._dict

    def list_attr(self, attr: str) -> List[Any]:
        """Return a list of the given attribute for every item."""
        return [getattr(item, attr) for item in self]

    def __contains__(self, item: Union[str, T]) -> bool:
        """Check if item or ID is in the DictList.

        Returns:
            True if the item or ID is present, False otherwise.
        """
        if isinstance(item, str):
            return item in self._dict
        else:
            return list.__contains__(self, item)

    def __getitem__(self, key: Union[int, slice, str]) -> Union[T, List[T]]:
        """Get item by index, slice, or ID.

        Returns:
            The item (for int/str keys) or a new ``DictList`` (for slices).
        """
        if isinstance(key, str):
            return self.get_by_id(key)
        elif isinstance(key, slice):
            # Return a new DictList for slices
            result = DictList[T]()
            result.extend(list.__getitem__(self, key))
            return result
        else:
            return list.__getitem__(self, key)

    def __setitem__(self, index: int, item: T) -> None:
        """Set item at index.

        Raises:
            TypeError: If the index is not an integer.
            ValueError: If the new item's ID conflicts with an existing
                ID at a different index.
        """
        if not isinstance(index, int):
            raise TypeError("DictList indices must be integers")

        # Remove old item's ID from dict
        old_item = self[index]
        if old_item.id in self._dict:
            del self._dict[old_item.id]

        # Check new item's ID doesn't conflict
        if item.id in self._dict and self._dict[item.id] != index:
            raise ValueError(
                f"ID '{item.id}' already exists at a different index"
            )

        # Set the item
        list.__setitem__(self, index, item)
        self._dict[item.id] = index

    def __delitem__(self, index: Union[int, slice]) -> None:
        """Delete item at index."""
        if isinstance(index, slice):
            # For slices, we need to regenerate the entire index
            list.__delitem__(self, index)
            self._generate_index()
        else:
            item = self[index]
            list.__delitem__(self, index)

            # Remove from dict and update subsequent indices
            del self._dict[item.id]
            for key, idx in list(self._dict.items()):
                if idx > index:
                    self._dict[key] = idx - 1

    def __repr__(self) -> str:
        """String representation of DictList with summary information.

        Returns:
            A multi-line string summarising the contents of the DictList.
        """
        if not self:
            return "DictList([])"

        # Detect item type from first item
        first_item = self[0]
        item_type = type(first_item).__name__

        lines = [f"=== {item_type}s ==="]

        # Type-specific summaries
        if item_type == "Metabolite":
            lines.extend(self._metabolite_summary())
        elif item_type == "Reaction":
            lines.extend(self._reaction_summary())
        else:
            # Generic summary
            sample_size = min(5, len(self))
            lines.append("  Sample items:")
            for item in list(self)[:sample_size]:
                lines.append(f"    - {item.id}")
            if len(self) > sample_size:
                lines.append(f"    ... and {len(self) - sample_size} more")

        return "\n".join(lines)

    def _metabolite_summary(self) -> List[str]:
        """Generate metabolite-specific summary lines.

        Returns:
            A list of summary lines describing the metabolites.
        """
        lines = []

        # Total count
        lines.append(f"  Total: {len(self)}")

        # Count metabolites by compartment and attributes
        compartments = {}
        with_atoms = 0
        with_formula = 0
        with_inchi = 0
        for m in self:
            comp = m.compartment or "unspecified"
            compartments[comp] = compartments.get(comp, 0) + 1
            if m.atoms and m.atoms > 0:
                with_atoms += 1
            if m.formula:
                with_formula += 1
            # Check for InChI in annotations
            if hasattr(m, "annotations") and m.annotations:
                for ann in m.annotations:
                    if ann.name and ann.name.lower() in ("inchi", "inchikey"):
                        with_inchi += 1
                        break

        # Show compartment distribution
        if len(compartments) > 1 or (
            len(compartments) == 1 and "unspecified" not in compartments
        ):
            lines.append("  By compartment:")
            for comp, count in sorted(compartments.items()):
                lines.append(f"    {comp}: {count}")

        lines.append(f"  With atom counts: {with_atoms}")
        lines.append(f"  With formula: {with_formula}")
        lines.append(f"  With InChI: {with_inchi}")

        # Sample metabolites
        sample_size = min(5, len(self))
        lines.append("")
        lines.append("  Sample metabolites:")
        for m in list(self)[:sample_size]:
            atoms_str = f", atoms={m.atoms}" if m.atoms else ""
            formula_str = f", formula={m.formula}" if m.formula else ""
            comp_str = f" [{m.compartment}]" if m.compartment else ""
            lines.append(f"    - {m.id}{comp_str}{atoms_str}{formula_str}")

        if len(self) > sample_size:
            lines.append(f"    ... and {len(self) - sample_size} more")

        return lines

    def _reaction_summary(self) -> List[str]:
        """Generate reaction-specific summary lines.

        Returns:
            A list of summary lines describing the reactions.
        """
        lines = []

        # Total count
        lines.append(f"  Total: {len(self)}")

        # Count reaction types
        reversible = sum(1 for r in self if r.reversibility)
        irreversible = len(self) - reversible
        variant_reactions = sum(1 for r in self if r.is_variant_reaction)
        total_variants = sum(
            r.n_variants for r in self if r.is_variant_reaction
        )

        lines.append(f"  Reversible: {reversible}")
        lines.append(f"  Irreversible: {irreversible}")
        if variant_reactions > 0:
            lines.append(
                f"  With atom transition variants: {variant_reactions} "
                f"({total_variants} total variants)"
            )

        # Sample reactions
        sample_size = min(5, len(self))
        lines.append("")
        lines.append("  Sample reactions:")
        for r in list(self)[:sample_size]:
            variant_str = (
                f" [{r.n_variants} variants]" if r.is_variant_reaction else ""
            )
            lines.append(f"    - {r.id}: {r.equation}{variant_str}")

        if len(self) > sample_size:
            lines.append(f"    ... and {len(self) - sample_size} more")

        return lines

    def __copy__(self) -> "DictList[T]":
        """Create a shallow copy of the DictList.

        Returns:
            A shallow copy of this DictList.
        """
        return DictList(self)

    def copy(self) -> "DictList[T]":
        """Create a shallow copy of the DictList.

        Returns:
            A shallow copy of this DictList.
        """
        return self.__copy__()

    @property
    def ids(self) -> List[str]:
        """Get list of all IDs in order."""
        return [item.id for item in self]

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> core_schema.CoreSchema:
        """Get Pydantic core schema for serialization/validation.

        Returns:
            A Pydantic core schema that validates as a list and returns
            a ``DictList``.
        """
        # Get the inner type from the generic
        if hasattr(source_type, "__args__") and source_type.__args__:
            inner_type = source_type.__args__[0]
        else:
            inner_type = Any

        # Get the schema for a list of the inner type
        list_schema = handler.generate_schema(List[inner_type])

        # Return a schema that validates as a list but returns a DictList
        return core_schema.no_info_after_validator_function(
            lambda v: cls(v),
            list_schema,
        )
