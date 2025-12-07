import logging
import os
import numpy as np
from pathlib import Path
from openpnm.models.phase import _phasedocs


logger = logging.getLogger(__name__)


__all__ = ["gaseous_species_in_water"]


@_phasedocs
def gaseous_species_in_water(
    phase,
    T="throat.temperature",
):
    r"""
    Calculate Henry's law constant for gaseous species dissolved in water

    Parameters
    ----------
    %(phase)s
    %(T)s

    Returns
    -------
    H : ndarray
        A numpy ndarray containing Henry's law constant (Kpx) [atm/mol-frac]

    Notes
    -----
    The constant for the correlation a lookup using the chemical formula
    stored in ``phase.params['formula']``.

    References
    ----------
    Yaws, Carl L., et al. "Solubility & Henry's Law constants for sulfur
    compounds in water: unlike traditional methods, the new correlation and
    data presented here are appropriate for very low concentrations."
    Chemical Engineering 110.8 (2003): 60-65.

    """
    import pandas as pd
    fname = "gas_water_henry.csv"
    path = Path(os.path.realpath(__file__), "../")
    path = Path(path.resolve(), fname)
    df = pd.read_csv(path)
    row = df[df.Formula == phase.params['formula']]
    A, B, C, D = row.iloc[0, 3:7].astype(float)
    Tmin, Tmax = row.iloc[0, 7:9].astype(float)
    T = phase[T]
    if (T.min() < Tmin) or (T.max() > Tmax):
        logger.critical("The correlation is only accurate for temperatures in the "
                        + f"range of {Tmin:.1f}K and {Tmax:.1f}K!")
    Hpx_log10 = A + B/T + C*np.log10(T) + D*T
    return 10**Hpx_log10
