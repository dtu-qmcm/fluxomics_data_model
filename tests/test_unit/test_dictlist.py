import pytest

from fluxomics_data_converter.core.common import DictList
from fluxomics_data_converter.model.metabolite import Metabolite


def _make_metabolite(mid, name=None, atoms=0):
    return Metabolite(id=mid, name=name, atoms=atoms)


class TestDictListCreation:
    def test_empty(self):
        dl = DictList()
        assert len(dl) == 0
        assert dl.ids == []

    def test_from_list(self):
        items = [_make_metabolite("a"), _make_metabolite("b")]
        dl = DictList(items)
        assert len(dl) == 2
        assert dl.ids == ["a", "b"]

    def test_from_dictlist(self):
        items = [_make_metabolite("a"), _make_metabolite("b")]
        dl1 = DictList(items)
        dl2 = DictList(dl1)
        assert len(dl2) == 2
        assert dl2.ids == ["a", "b"]


class TestDictListAppend:
    def test_append(self):
        dl = DictList()
        dl.append(_make_metabolite("glc"))
        assert len(dl) == 1
        assert dl.ids == ["glc"]

    def test_append_duplicate_raises(self):
        dl = DictList([_make_metabolite("glc")])
        with pytest.raises(ValueError, match="already present"):
            dl.append(_make_metabolite("glc"))


class TestDictListExtend:
    def test_extend(self):
        dl = DictList([_make_metabolite("a")])
        dl.extend([_make_metabolite("b"), _make_metabolite("c")])
        assert len(dl) == 3
        assert dl.ids == ["a", "b", "c"]

    def test_extend_duplicate_raises(self):
        dl = DictList([_make_metabolite("a")])
        with pytest.raises(ValueError):
            dl.extend([_make_metabolite("a")])


class TestDictListGetById:
    def test_get_by_id(self):
        dl = DictList([_make_metabolite("glc", name="glucose")])
        item = dl.get_by_id("glc")
        assert item.name == "glucose"

    def test_get_by_id_missing(self):
        dl = DictList([_make_metabolite("glc")])
        with pytest.raises(KeyError):
            dl.get_by_id("missing")


class TestDictListHasId:
    def test_has_id_true(self):
        dl = DictList([_make_metabolite("glc")])
        assert dl.has_id("glc") is True

    def test_has_id_false(self):
        dl = DictList([_make_metabolite("glc")])
        assert dl.has_id("missing") is False


class TestDictListContains:
    def test_contains_by_string(self):
        dl = DictList([_make_metabolite("glc")])
        assert "glc" in dl

    def test_contains_by_object(self):
        m = _make_metabolite("glc")
        dl = DictList([m])
        assert m in dl


class TestDictListIndex:
    def test_index_by_string(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        assert dl.index("b") == 1

    def test_index_by_object(self):
        m = _make_metabolite("a")
        dl = DictList([m])
        assert dl.index(m) == 0

    def test_index_missing_raises(self):
        dl = DictList([_make_metabolite("a")])
        with pytest.raises(ValueError):
            dl.index("missing")


class TestDictListGetItem:
    def test_getitem_by_index(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        assert dl[1].id == "b"

    def test_getitem_by_string(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        assert dl["b"].id == "b"

    def test_getitem_slice(self):
        dl = DictList(
            [
                _make_metabolite("a"),
                _make_metabolite("b"),
                _make_metabolite("c"),
            ]
        )
        result = dl[1:]
        assert len(result) == 2
        assert result.ids == ["b", "c"]


class TestDictListRemove:
    def test_remove_by_string(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        dl.remove("a")
        assert len(dl) == 1
        assert dl.ids == ["b"]

    def test_remove_by_object(self):
        m = _make_metabolite("a")
        dl = DictList([m, _make_metabolite("b")])
        dl.remove(m)
        assert len(dl) == 1


class TestDictListPop:
    def test_pop_last(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        item = dl.pop()
        assert item.id == "b"
        assert len(dl) == 1

    def test_pop_index(self):
        dl = DictList(
            [
                _make_metabolite("a"),
                _make_metabolite("b"),
                _make_metabolite("c"),
            ]
        )
        item = dl.pop(1)
        assert item.id == "b"
        assert dl.ids == ["a", "c"]


class TestDictListIds:
    def test_ids_property(self):
        dl = DictList([_make_metabolite("x"), _make_metabolite("y")])
        assert dl.ids == ["x", "y"]


class TestDictListCopy:
    def test_copy(self):
        dl = DictList([_make_metabolite("a")])
        cp = dl.copy()
        assert len(cp) == 1
        assert cp.ids == ["a"]
        assert cp is not dl


class TestDictListSetItem:
    def test_setitem(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        dl[0] = _make_metabolite("c")
        assert dl.ids == ["c", "b"]

    def test_setitem_duplicate_at_different_index_raises(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        with pytest.raises(ValueError):
            dl[0] = _make_metabolite("b")


class TestDictListListAttr:
    def test_list_attr(self):
        dl = DictList(
            [_make_metabolite("a", atoms=3), _make_metabolite("b", atoms=6)]
        )
        assert dl.list_attr("atoms") == [3, 6]


class TestDictListGetByAny:
    def test_get_by_any_string(self):
        dl = DictList([_make_metabolite("glc", name="glucose")])
        assert dl.get_by_any("glc").name == "glucose"

    def test_get_by_any_int(self):
        dl = DictList([_make_metabolite("a"), _make_metabolite("b")])
        assert dl.get_by_any(1).id == "b"

    def test_get_by_any_object(self):
        m = _make_metabolite("a")
        dl = DictList([m])
        assert dl.get_by_any(m).id == "a"
