from OpenPNM.Base import Tools


class ToolsTest:
    class PrintableListTest:
        def setup_class(self):
            self.list = Tools.PrintableList(['first', 'second'])

        def test_str(self):
            actual_string = self.list.__str__()
            expected_string = '\n' + \
                '------------------------------------------------------------\n' + \
                '1\t: first' + \
                '2\t: second' + \
                '------------------------------------------------------------\n'
            assert actual_string == expected_string
