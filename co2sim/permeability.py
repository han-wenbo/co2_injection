"""
Permeability fields k(x) and k(x, y): how easily fluid flows through each cell.

Real rock is not uniform. These helpers build a few typical situations:
sealing faults (thin walls of almost impermeable rock) and a randomly
varying ("heterogeneous") sandstone. Start from a uniform field, e.g.
np.full(n_cells, params.permeability), and add features to it.
"""

import numpy as np
from scipy.ndimage import gaussian_filter

# A sealing fault is filled with crushed, clay-rich rock: about 1e-18 m²,
# i.e. a hundred thousand times less permeable than the sandstone around it.
SEALING_FAULT_PERMEABILITY = 1e-18  # [m²]


def add_fault_1d(k, x, fault_x, fault_width, fault_permeability=SEALING_FAULT_PERMEABILITY):
    """Copy of the 1D field k with a fault of width `fault_width` [m] centred at x = fault_x.

    x : cell centres [m] of the 1D grid
    """
    k = k.copy()
    k[np.abs(x - fault_x) <= fault_width / 2] = fault_permeability
    return k


def add_fault_2d(k, x, y, fault_x, fault_y_start, fault_y_end, fault_permeability=SEALING_FAULT_PERMEABILITY):
    """Copy of the 2D field k with a straight north–south fault at x = fault_x.

    The fault is one cell wide and runs from y = fault_y_start to fault_y_end.
    If it does not reach the edges, fluid can flow around its ends.
    x, y : cell centres [m] of the 2D grid; k has shape (len(y), len(x))
    """
    k = k.copy()
    fault_column = np.argmin(np.abs(x - fault_x))
    fault_rows = (y >= fault_y_start) & (y <= fault_y_end)
    k[fault_rows, fault_column] = fault_permeability
    return k


def random_heterogeneous_2d(shape, mean, spread_decades=0.5, correlation_cells=8.0, seed=1):
    """Random but smooth permeability field (log-normal), like a real sandstone.

    log10(k) = log10(mean) + spread_decades · (smooth random noise with std 1)

    spread_decades : 0.5 means k typically varies by a factor 10^0.5 ≈ 3 up or down
    correlation_cells : size of the "patches" of similar rock, in cells
    seed : fixed seed, so every run uses the same field
    """
    rng = np.random.default_rng(seed)
    noise = gaussian_filter(rng.standard_normal(shape), sigma=correlation_cells)
    noise /= noise.std()  # smoothing shrinks the noise; scale it back to std 1
    return mean * 10 ** (spread_decades * noise)
