"""
Exact solutions of the pressure equation, used to check the numerical solvers.

Both assume: uniform rock, an aquifer that is infinitely large, a single well
switched on at t = 0 with a constant rate, and Δp = 0 at t = 0.
"""

import numpy as np
from scipy.special import erfc, exp1


def overpressure_1d(x, t, volume_rate, flow_area, params):
    """Overpressure [Pa] at distance x [m] from a well in an infinite 1D aquifer.

    The well injects `volume_rate` [m³/s] through a cross-section `flow_area`
    [m²]; half of it flows to each side. Solution (heat equation with a
    constant point source):

        Δp = (Q μ / (A k)) · [ sqrt(D t / π) · exp(−x² / 4Dt) − |x|/2 · erfc(|x| / (2 sqrt(D t))) ]
    """
    x = np.abs(np.asarray(x, dtype=float))
    diffusivity = params.diffusivity
    spread = np.sqrt(diffusivity * t)  # how far the pressure has travelled [m]

    scale = volume_rate * params.viscosity / (flow_area * params.permeability)  # [Pa/m]
    shape = spread / np.sqrt(np.pi) * np.exp(-((x / (2 * spread)) ** 2)) - x / 2 * erfc(x / (2 * spread))
    return scale * shape


def overpressure_theis(r, t, volume_rate, params):
    """Overpressure [Pa] at distance r [m] from a well in an infinite 2D aquifer.

    Theis (1935), the classic solution of well hydraulics:

        Δp = Q μ / (4 π k H) · E1( r² / (4 D t) )

    E1 is the exponential integral (scipy.special.exp1).
    """
    r = np.asarray(r, dtype=float)
    u = r**2 / (4 * params.diffusivity * t)
    scale = volume_rate * params.viscosity / (4 * np.pi * params.permeability * params.thickness)
    return scale * exp1(u)
