import pytest
import numpy as np
from openpnm.materials import BundleOfTubes


class BundleOfTubesTest:

    def setup_class(self):
        pass

    def test_instantiate_with_defaults(self):
        net, geo, phase = BundleOfTubes(shape=[30, 30])
        assert net.Np == 1800
        assert net.Nt == (net.Np / 2)

    def test_instantiate_custom_shape_and_size(self):
        _ = BundleOfTubes(
            shape=30,
            spacing=0.001,
            length=0.01,
            psd_params={
                "distribution": "weibull",
                "scale": 0.0002,
                "shape": 2,
                "loc": 1e-12,
            },
        )

    def test_instantiate_custom_distribution(self):
        _ = BundleOfTubes(
            shape=30,
            spacing=0.001,
            length=0.01,
            psd_params={
                "distribution": "lognorm",
                "scale": 0.0002,
                "s": 2,
                "loc": 1e-12,
            },
        )

    def test_instantiate_with_exceptions(self):
        with pytest.raises(Exception):
            _ = BundleOfTubes(shape=[30, 30, 30], spacing=0.001, length=0.01)
        with pytest.raises(Exception):
            _ = BundleOfTubes(shape=30, spacing=[0.001, 0.001], length=0.01)

    def test_instantiate_with_settings(self):
        net, geo, phase = BundleOfTubes(
            shape=30,
            spacing=0.001,
            length=0.01,
            psd_params={"distribution": "normal", "loc": 0.02, "scale": 2},
            settings={"adjust_psd": "clip"},
        )
        assert geo["throat.size_distribution"].max() > geo["throat.diameter"].max()

        net, geo, phase = BundleOfTubes(
            shape=30,
            spacing=0.001,
            length=0.01,
            psd_params={"distribution": "norm", "loc": 0.02, "scale": 0.001},
            settings={"adjust_psd": "clip"},
        )
        assert geo["throat.size_distribution"].max() > geo["throat.diameter"].max()

        net, geo, phase = BundleOfTubes(
            shape=30,
            spacing=1.0,
            length=0.01,
            psd_params={"distribution": "normal", "loc": 0.02, "scale": 0.001},
            settings={"adjust_psd": ''},
        )
        assert np.all(geo["throat.size_distribution"] == geo["throat.diameter"])


if __name__ == "__main__":

    t = BundleOfTubesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
    self = t
