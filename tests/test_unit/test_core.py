import pytest
from datetime import datetime

from fluxomics_data_converter.core.core import (
    FluxomicsData,
    Metadata,
    MetabolicNetworkModel,
    LabelingExperiments,
)
from fluxomics_data_converter.model.metabolite import Metabolite
from fluxomics_data_converter.model.reaction import Reaction


def _make_model():
    m1 = Metabolite(id="A", atoms=3)
    m2 = Metabolite(id="B", atoms=3)
    r1 = Reaction(id="r1", reactants=["A"], products=["B"])
    return MetabolicNetworkModel(metabolites=[m1, m2], reactions=[r1])


class TestMetadata:
    def test_creation(self):
        meta = Metadata(name="test", version="1.0")
        assert meta.name == "test"
        assert meta.version == "1.0"

    def test_date_parsing(self):
        meta = Metadata(date="2024-01-15 10:30:00")
        assert isinstance(meta.date, datetime)
        assert meta.date.year == 2024

    def test_all_optional(self):
        meta = Metadata()
        assert meta.name is None
        assert meta.date is None


class TestModel:
    def test_creation(self):
        model = _make_model()
        assert len(model.metabolites) == 2
        assert len(model.reactions) == 1

    def test_metabolite_ids(self):
        model = _make_model()
        assert model.metabolite_ids == frozenset(["A", "B"])

    def test_reaction_ids(self):
        model = _make_model()
        assert model.reaction_ids == frozenset(["r1"])

    def test_invalid_reactant_reference(self):
        with pytest.raises(ValueError, match="references unknown metabolite"):
            MetabolicNetworkModel(
                metabolites=[Metabolite(id="A")],
                reactions=[Reaction(id="r1", reactants=["A"], products=["C"])],
            )

    def test_invalid_product_reference(self):
        with pytest.raises(ValueError, match="references unknown metabolite"):
            MetabolicNetworkModel(
                metabolites=[Metabolite(id="A")],
                reactions=[Reaction(id="r1", reactants=["C"], products=["A"])],
            )

    def test_empty_model(self):
        model = MetabolicNetworkModel()
        assert len(model.metabolites) == 0
        assert len(model.reactions) == 0


class TestLabelingExperiments:
    def test_creation(self):
        exp = LabelingExperiments(name="exp1")
        assert exp.name == "exp1"
        assert exp.stationary is True

    def test_non_stationary(self):
        exp = LabelingExperiments(name="exp2", stationary=False)
        assert exp.stationary is False

    def test_traced_metabolites(self):
        from fluxomics_data_converter.experiment.tracer import (
            Tracers,
            LabelComposition,
        )

        exp = LabelingExperiments(
            name="exp1",
            tracers=[
                Tracers(
                    metabolite="glc",
                    labels=[LabelComposition(labeled_pattern="111111")],
                ),
            ],
        )
        assert exp.traced_metabolites == frozenset(["glc"])


class TestFluxomicsData:
    def test_creation(self):
        model = _make_model()
        dm = FluxomicsData(model=model)
        assert dm.model is model
        assert dm.metadata is None
        assert dm.experiments == []

    def test_with_metadata(self):
        model = _make_model()
        meta = Metadata(name="test")
        dm = FluxomicsData(model=model, metadata=meta)
        assert dm.metadata.name == "test"

    def test_metabolite_ids_property(self):
        model = _make_model()
        dm = FluxomicsData(model=model)
        assert dm.metabolite_ids == frozenset(["A", "B"])

    def test_reaction_ids_property(self):
        model = _make_model()
        dm = FluxomicsData(model=model)
        assert dm.reaction_ids == frozenset(["r1"])

    def test_repr(self):
        model = _make_model()
        meta = Metadata(name="test")
        dm = FluxomicsData(model=model, metadata=meta)
        repr_str = repr(dm)
        assert "Fluxomics Data Converter Summary" in repr_str
        assert "test" in repr_str

    def test_get_experiments(self):
        model = _make_model()
        exp = LabelingExperiments(name="exp1")
        dm = FluxomicsData(model=model, experiments=[exp])
        assert dm.get_experiments("exp1") is exp
        assert dm.get_experiments("missing") is None

    def test_experiments_names(self):
        model = _make_model()
        dm = FluxomicsData(
            model=model,
            experiments=[
                LabelingExperiments(name="a"),
                LabelingExperiments(name="b"),
            ],
        )
        assert dm.experiments_names == frozenset(["a", "b"])

    def test_duplicate_experiment_names_raises(self):
        model = _make_model()
        with pytest.raises(ValueError, match="unique"):
            FluxomicsData(
                model=model,
                experiments=[
                    LabelingExperiments(name="x"),
                    LabelingExperiments(name="x"),
                ],
            )
