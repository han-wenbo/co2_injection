"""
Building blocks shared by the 1D and the 2D model (finite-volume method).

The idea in one paragraph
-------------------------
Cut the aquifer into small cells. Each cell stores one number: its
*overpressure* Δp = p − p0, the pressure above the natural pressure p0.
Two neighbouring cells exchange fluid through the face between them:

    flow from cell a to cell b = T_ab · (Δp_a − Δp_b)          [m³/s]

This is Darcy's law: fluid moves from high to low pressure. T_ab is the
*transmissibility* of the face [m³/(Pa·s)]. It is large when the rock is
permeable (k), the fluid is runny (1/μ), the face is big (its area) and the
two cell centres are close (1/distance).

A cell's pressure then rises by

    (fluid flowing in − fluid flowing out + fluid from a well) / storage

where storage = φ · c_t · cell volume [m³/Pa] is how much extra fluid fits
into the cell per pascal of pressure increase.
"""

import numpy as np
from scipy import sparse


def cell_centres(length, n_cells):
    """Positions of the centres of `n_cells` equal cells covering [0, length] [m]."""
    return (np.arange(n_cells) + 0.5) * (length / n_cells)


def face_transmissibility(k_face, viscosity, face_area, distance):
    """Transmissibility T = (k / μ) · area / distance  [m³/(Pa·s)] (Darcy's law)."""
    return k_face / viscosity * face_area / distance


def harmonic_mean(k_a, k_b):
    """Permeability of the face between two cells with permeabilities k_a and k_b.

    We use the harmonic mean, not the ordinary average: fluid has to pass
    through *both* half-cells in a row, so an almost impermeable cell blocks
    the face. Example: harmonic_mean(100, 0.001) ≈ 0.002, not 50.
    """
    return 2 * k_a * k_b / (k_a + k_b)


def assemble_flow_matrix(n_cells, cell_a, cell_b, transmissibility, boundary_transmissibility):
    """Build the sparse "flow matrix" L of a grid of cells.

    L is defined by  (L @ Δp)[i] = net flow OUT of cell i  [m³/s],
    i.e. L @ Δp == -net_inflow(Δp). Implicit time stepping needs it.

    Parameters
    ----------
    n_cells : number of cells (unknowns)
    cell_a, cell_b : for every face, the indices of the two cells it connects
    transmissibility : for every face, its transmissibility T
    boundary_transmissibility : for every cell, the transmissibility towards
        the outside of the model (0 for cells not on an open boundary)

    Every face contributes
        +T to L[a, a] and L[b, b]   (cell a and b lose fluid when their own pressure is high)
        -T to L[a, b] and L[b, a]   (and gain fluid when the neighbour's pressure is high)
    """
    rows = np.concatenate([cell_a, cell_b, cell_a, cell_b])
    cols = np.concatenate([cell_a, cell_b, cell_b, cell_a])
    values = np.concatenate([transmissibility, transmissibility, -transmissibility, -transmissibility])

    # coo_matrix adds up entries that land on the same (row, col): exactly what we need
    face_part = sparse.coo_matrix((values, (rows, cols)), shape=(n_cells, n_cells))
    boundary_part = sparse.diags(boundary_transmissibility)
    return (face_part + boundary_part).tocsc()
