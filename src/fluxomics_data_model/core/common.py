"""
Common data structures used across FluxML models.
"""

from typing import Optional, Union, Any, List, TypeVar, Generic
from pydantic import BaseModel, Field, GetCoreSchemaHandler
from pydantic_core import core_schema
import jax.numpy as jnp
from itertools import islice


class Annotation(BaseModel):
    """
    FluxML annotation element for additional metadata.

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
    """
    FluxML textual or MathML content.

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
        super().__init__(**data)
        if not self.textual and not self.mathml:
            raise ValueError("Either textual or mathml must be provided")


class ErrorModel(BaseModel):
    """
    FluxML error model for measurement uncertainties.

    Corresponds to fluxml/experiments/measurement/model/*/errormodel
    """

    expression: TextualOrMath = Field(description="Error model expression")

    class Config:
        frozen = True


class JAXArray(BaseModel):
    """
    JAX array wrapper for Pydantic serialization.
    """

    shape: tuple[int, ...] = Field(description="Array shape")
    dtype: str = Field(description="Array dtype")
    data: list = Field(description="Array data as nested lists")

    class Config:
        frozen = True

    @classmethod
    def from_jax_array(cls, arr: jnp.ndarray) -> "JAXArray":
        """Create from JAX array."""
        return cls(shape=arr.shape, dtype=str(arr.dtype), data=arr.tolist())

    def to_jax_array(self) -> jnp.ndarray:
        """Convert to JAX array."""
        return jnp.array(self.data, dtype=self.dtype).reshape(self.shape)


class TimeSeries(BaseModel):
    """
    Time series data structure for non-stationary measurements.
    """

    times: JAXArray = Field(description="Time points")
    values: JAXArray = Field(description="Measurement values")
    errors: Optional[JAXArray] = Field(
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
        """Create from JAX arrays."""
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
        """Convert to JAX arrays."""
        return (
            self.times.to_jax_array(),
            self.values.to_jax_array(),
            self.errors.to_jax_array() if self.errors else None,
        )


# Type variable for DictList items
T = TypeVar("T", bound=BaseModel)


class DictList(list, Generic[T]):
    """
    A combined dict and list data structure.

    This object behaves like a list but has O(1) speed benefits
    of a dict when looking up elements by their id attribute.

    Items must have an 'id' attribute.
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
        """Check if an ID already exists in the DictList."""
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
        """Remove and return item at index (default last)."""
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
        """Return index of item in the list."""
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
        """Get item by its ID attribute."""
        if item_id not in self._dict:
            raise KeyError(f"No item with ID '{item_id}'")
        return self[self._dict[item_id]]

    def get_by_any(self, key: Union[str, int, T]) -> T:
        """Get item by ID, index, or the item itself."""
        if isinstance(key, int):
            return self[key]
        elif isinstance(key, str):
            return self.get_by_id(key)
        elif key in self:
            return key
        else:
            raise ValueError(f"Item {key} not found in DictList")

    def has_id(self, item_id: str) -> bool:
        """Check if an ID exists in the DictList."""
        return item_id in self._dict

    def list_attr(self, attr: str) -> List[Any]:
        """Return a list of the given attribute for every item."""
        return [getattr(item, attr) for item in self]

    def __contains__(self, item: Union[str, T]) -> bool:
        """Check if item or ID is in the DictList."""
        if isinstance(item, str):
            return item in self._dict
        else:
            return list.__contains__(self, item)

    def __getitem__(self, key: Union[int, slice, str]) -> Union[T, List[T]]:
        """Get item by index, slice, or ID."""
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
        """Set item at index."""
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
        """String representation of DictList."""
        return f"DictList({list.__repr__(self)})"

    def __copy__(self) -> "DictList[T]":
        """Create a shallow copy of the DictList."""
        return DictList(self)

    def copy(self) -> "DictList[T]":
        """Create a shallow copy of the DictList."""
        return self.__copy__()

    @property
    def ids(self) -> List[str]:
        """Get list of all IDs in order."""
        return [item.id for item in self]

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> core_schema.CoreSchema:
        """Get Pydantic core schema for serialization/validation."""
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
