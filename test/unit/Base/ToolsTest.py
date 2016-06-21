from OpenPNM.Base import Tools


class PrintableListTest:
    def setup_class(self):
        pass

    def test_str(self):
        list_ = Tools.PrintableList(['first', 'second'])
        actual_string = list_.__str__()
        expected_string = \
            '------------------------------------------------------------\n' + \
            '1\t: first\n' + \
            '2\t: second\n' + \
            '------------------------------------------------------------'
        assert actual_string == expected_string

    def test_hidden_items(self):
        list_ = Tools.PrintableList(['pore.first', 'pore._second'])
        actual_string = list_.__str__()
        assert 'pore.first' in actual_string
        assert 'pore._second' not in actual_string


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
