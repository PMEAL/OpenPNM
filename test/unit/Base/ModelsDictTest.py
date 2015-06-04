from OpenPNM.Base import ModelsDict
from OpenPNM.Base.__ModelsDict__ import ModelWrapper


class ModelsDictTest:

    class ModelWrapperTest:
        def test_str_with_none_model(self):
            actual_string = ModelWrapper(model=None).__str__()
            expected_string = 'No model specified.'
            assert expected_string == actual_string

        def test_find_master_with_more_than_one_master(self):
            pass

    class ModelsDictTest:
        def setup_class(self):
            self.models_dict = ModelsDict()

        def test_str(self):
            actual_string = self.models_dict.__str__()
            expected_string = \
                '------------------------------------------------------------\n' + \
                '#     Property Name                  Regeneration Mode\n' + \
                '------------------------------------------------------------\n' + \
                '------------------------------------------------------------'
            assert expected_string == actual_string

        def test_str_with_models(self):
            pass
