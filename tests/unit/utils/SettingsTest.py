from traits.api import Int, List, Str, Float, TraitError, ListStr
import openpnm as op
from openpnm.utils import SettingsAttr, TypedList, TypedSet
import pytest


class SettingsTest:

    def setup_class(self): ...

    def test_standard_initialization(self):
        class S1:
            r"""
            This is a docstring
            """
            a = 1
            b = 2
            d = TypedList(types=[str])

        sets1 = SettingsAttr(S1)
        assert "This is a docstring" in sets1.__doc__

    def test_inheritance_and_immutability(self):
        class S2:
            r"""
            This is a docstring
            """
            a = 2
            b = 3

        class S3(S2):
            r"""
            Different docstring
            """
            b = 3
            c = 4

        sets2 = SettingsAttr(S2)
        sets3 = SettingsAttr(S3)
        sets3.b = 44
        assert sets3._attrs == ['a', 'b', 'c']
        assert sets2.b != sets3.b
        assert "Different docstring" in sets3.__doc__

        with pytest.raises(Exception):
            sets3.b = 'nope'

    def test_adding_new_attr_type_is_enforced(self):
        class S4:
            r"""
            This is a docstring
            """
            a = 2
            b = 3.5

        sets4 = SettingsAttr(S4)
        sets4.c = 1.5
        with pytest.raises(Exception):
            sets4.c = "string"

    def test_update_from_dataclass(self):
        class S7:
            r"""
            This is a docstring
            """
            a = 1
            b = 2

        class Data:
            a = 22
            c = 5.5

        sets7 = SettingsAttr(S7)
        assert sets7.a == 1
        sets7._update(Data())
        assert sets7.a == 22

    def test_initial_None(self):
        class S8:
            r"""
            This is a docstring
            """
            a = 1
            b = None

        sets8 = SettingsAttr(S8)
        sets8.b = 2.2
        assert sets8.b == 2.2
        with pytest.raises(Exception):
            sets8.b = 'str'

    def test_typed_set_inferred_type_after_init(self):
        s = TypedSet()
        s.add(0)
        s2 = set((0, 2))
        assert s2.difference(s) == set([2])
        s.add(2)
        assert s2.difference(s) == set()
        s.add(2)
        assert len(s) == 2  # Ensure length stays the same
        with pytest.raises(Exception):
            s.add('1')

    def test_typed_set_given_multiple_types_during_init(self):
        s = TypedSet(types=[int, float])
        s.add(1)
        s.add(2.0)
        assert len(s) == 2
        with pytest.raises(Exception):
            s.add([2])


if __name__ == '__main__':

    t = SettingsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
