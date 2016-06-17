from OpenPNM.Base import Tools


class ToolsTest:
    class PrintableListTest:
        def setup_class(self):
            self.list = Tools.PrintableList(['first', 'second'])

        def test_str(self):
            actual_string = self.list.__str__()
            expected_string = \
                '------------------------------------------------------------\n' + \
                '1\t: first\n' + \
                '2\t: second\n' + \
                '------------------------------------------------------------'
            assert actual_string == expected_string

        def test_hidden_items(self):
            list_ = Tools.PrintableList(['pore.first', 'pore._second'])
            s = list_.__str__()
            assert 'pore.first' in s
            assert 'pore._second' not in s

    class PrintableDictTest:
        def setup_class(self):
            self.dict = Tools.PrintableDict()
            self.dict['first_key'] = 'first_value'
            self.dict['second_key'] = 'second_value'

        def test_str(self):
            actual_string = self.dict.__str__()
            expected_string = \
                '------------------------------------------------------------\n' + \
                'key                       value\n' + \
                '------------------------------------------------------------\n' + \
                'first_key                 first_value\n' + \
                'second_key                second_value\n' + \
                '------------------------------------------------------------'
            assert actual_string == expected_string

        def test_representation(self):
            a = self.dict.__repr__()
            assert type(a) is str
