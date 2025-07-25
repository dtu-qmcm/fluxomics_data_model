"""
FluxML measurement definitions.
"""

from typing import Optional, List, Dict, Literal
from pydantic import BaseModel, Field, field_validator
import jax.numpy as jnp
from .common import TextualOrMath, ErrorModel, JAXArray, TimeSeries


class Group(BaseModel):
    """
    FluxML measurement group for mass spectrometry data.

    Corresponds to fluxml/experiments/measurement/model/
    labelingmeasurement/group
    """

    id: str = Field(description="Group identifier")
    times: Optional[str] = Field(default=None, description="Time points")
    scale: Literal["auto", "one"] = Field(
        default="auto", description="Scaling method"
    )
    errormodel: Optional[ErrorModel] = Field(
        default=None, description="Error model"
    )
    expression: TextualOrMath = Field(description="Measurement expression")

    # JAX-compatible numerical representation
    time_points_array: Optional[JAXArray] = Field(
        default=None, exclude=True, description="Time points array"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def time_points(self) -> Optional[jnp.ndarray]:
        """Get time points as JAX array."""
        if self.time_points_array is None:
            if self.times is None:
                return None
            # Parse time string - simplified implementation
            try:
                times = [float(t.strip()) for t in self.times.split(",")]
                return jnp.array(times)
            except Exception:
                return None
        return self.time_points_array.to_jax_array()


class NetFlux(BaseModel):
    """
    FluxML net flux measurement.

    Corresponds to fluxml/experiments/measurement/model/
    fluxmeasurement/netflux
    """

    id: str = Field(description="Net flux identifier")
    errormodel: Optional[ErrorModel] = Field(
        default=None, description="Error model"
    )
    expression: TextualOrMath = Field(description="Net flux expression")

    class Config:
        frozen = True
        extra = "forbid"


class XchFlux(BaseModel):
    """
    FluxML exchange flux measurement.

    Corresponds to fluxml/experiments/measurement/model/
    fluxmeasurement/xchflux
    """

    id: str = Field(description="Exchange flux identifier")
    errormodel: Optional[ErrorModel] = Field(
        default=None, description="Error model"
    )
    expression: TextualOrMath = Field(description="Exchange flux expression")

    class Config:
        frozen = True
        extra = "forbid"


class MetaboliteSize(BaseModel):
    """
    FluxML metabolite size measurement.

    Corresponds to fluxml/experiments/measurement/model/
    metabolitesizemeasurement/metabolitesize
    """

    id: str = Field(description="Metabolite size identifier")
    errormodel: Optional[ErrorModel] = Field(
        default=None, description="Error model"
    )
    expression: TextualOrMath = Field(description="Metabolite size expression")

    class Config:
        frozen = True
        extra = "forbid"


class LabelingMeasurement(BaseModel):
    """
    FluxML labeling measurement collection.

    Corresponds to fluxml/experiments/measurement/model/labelingmeasurement
    """

    groups: List[Group] = Field(
        default_factory=list, description="Measurement groups"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @field_validator("groups")
    @classmethod
    def validate_unique_ids(cls, v: List[Group]) -> List[Group]:
        """Validate group IDs are unique."""
        ids = [group.id for group in v]
        if len(set(ids)) != len(ids):
            raise ValueError("Group IDs must be unique")
        return v


class FluxMeasurement(BaseModel):
    """
    FluxML flux measurement collection.

    Corresponds to fluxml/experiments/measurement/model/fluxmeasurement
    """

    net_fluxes: List[NetFlux] = Field(
        default_factory=list, description="Net flux measurements"
    )
    xch_fluxes: List[XchFlux] = Field(
        default_factory=list, description="Exchange flux measurements"
    )

    class Config:
        frozen = True
        extra = "forbid"


class MetaboliteSizeMeasurement(BaseModel):
    """
    FluxML metabolite size measurement collection.

    Corresponds to fluxml/experiments/measurement/model/
    metabolitesizemeasurement
    """

    metabolite_sizes: List[MetaboliteSize] = Field(
        default_factory=list, description="Metabolite size measurements"
    )

    class Config:
        frozen = True
        extra = "forbid"


class MeasurementModel(BaseModel):
    """
    FluxML measurement model.

    Corresponds to fluxml/experiments/measurement/model
    """

    labeling_measurement: Optional[LabelingMeasurement] = Field(
        default=None, description="Labeling measurements"
    )
    flux_measurement: Optional[FluxMeasurement] = Field(
        default=None, description="Flux measurements"
    )
    metabolitesize_measurement: Optional[MetaboliteSizeMeasurement] = Field(
        default=None, description="Metabolite size measurements"
    )

    class Config:
        frozen = True
        extra = "forbid"


class Datum(BaseModel):
    """
    FluxML measurement datum.

    Corresponds to fluxml/experiments/measurement/data/datum
    """

    id: str = Field(description="Datum identifier")
    value: float = Field(description="Measurement value")
    stddev: float = Field(description="Standard deviation")
    row: Optional[int] = Field(
        default=None, ge=1, le=256, description="Row index"
    )
    time: Optional[float] = Field(default=None, description="Time point")
    weight: Optional[str] = Field(
        default=None, description="Weight specification"
    )
    pos: Optional[int] = Field(
        default=None, ge=0, le=1024, description="Position"
    )
    type: Optional[Literal["S", "DL", "DR", "DD", "T"]] = Field(
        default=None, description="Data type"
    )

    class Config:
        frozen = True
        extra = "forbid"


class MeasurementData(BaseModel):
    """
    FluxML measurement data collection.

    Corresponds to fluxml/experiments/measurement/data
    """

    data: List[Datum] = Field(
        default_factory=list, description="Measurement data points"
    )

    # JAX-compatible numerical representation
    values_array: Optional[JAXArray] = Field(
        default=None, exclude=True, description="Values array"
    )
    errors_array: Optional[JAXArray] = Field(
        default=None, exclude=True, description="Errors array"
    )
    times_array: Optional[JAXArray] = Field(
        default=None, exclude=True, description="Times array"
    )

    class Config:
        frozen = True
        extra = "forbid"

    @property
    def values(self) -> jnp.ndarray:
        """Get measurement values as JAX array."""
        if self.values_array is None:
            return jnp.array([datum.value for datum in self.data])
        return self.values_array.to_jax_array()

    @property
    def errors(self) -> jnp.ndarray:
        """Get measurement errors as JAX array."""
        if self.errors_array is None:
            return jnp.array([datum.stddev for datum in self.data])
        return self.errors_array.to_jax_array()

    @property
    def times(self) -> Optional[jnp.ndarray]:
        """Get measurement times as JAX array."""
        if self.times_array is None:
            times = [
                datum.time for datum in self.data if datum.time is not None
            ]
            if not times:
                return None
            return jnp.array(times)
        return self.times_array.to_jax_array()

    def get_data_for_id(self, datum_id: str) -> List[Datum]:
        """Get all data points for a specific ID."""
        return [datum for datum in self.data if datum.id == datum_id]

    def to_time_series(self, datum_id: str) -> Optional[TimeSeries]:
        """Convert data for specific ID to time series."""
        data_points = self.get_data_for_id(datum_id)
        if not data_points:
            return None

        times = [dp.time for dp in data_points if dp.time is not None]
        values = [dp.value for dp in data_points]
        errors = [dp.stddev for dp in data_points]

        if not times:
            return None

        return TimeSeries.from_arrays(
            jnp.array(times), jnp.array(values), jnp.array(errors)
        )


class Measurement(BaseModel):
    """
    FluxML measurement specification.

    Corresponds to fluxml/experiments/measurement
    """

    model: MeasurementModel = Field(description="Measurement model")
    data: MeasurementData = Field(description="Measurement data")

    class Config:
        frozen = True
        extra = "forbid"

    def get_labeling_data(self) -> Optional[Dict[str, TimeSeries]]:
        """Get labeling measurement data as time series."""
        if not self.model.labeling_measurement:
            return None

        result = {}
        for group in self.model.labeling_measurement.groups:
            time_series = self.data.to_time_series(group.id)
            if time_series:
                result[group.id] = time_series

        return result if result else None
