"""
1D model: the aquifer as a row of cells along x.

        well
         ↓
    | 0 | 1 | 2 | ... | i | ... | n-1 |
    x = 0                       x = length

Each cell is dx long, `width` wide and `thickness` high. Pressure can only
change along x, which is what "1D" means. Think of a long, narrow aquifer
channel: the injected fluid can only escape to the left or to the right.
"""

import numpy as np

from .finite_volume import assemble_flow_matrix, cell_centres, face_transmissibility, harmonic_mean


class Aquifer1D:
    def __init__(self, params, length, width, n_cells, wells, permeability=None, boundary="open"):
        """
        Parameters
        ----------
        params : AquiferParameters
        length : length of the aquifer along x [m]
        width : width of the aquifer strip [m] (flow happens through width × thickness)
        n_cells : number of cells
        wells : list of Well
        permeability : k of every cell [m²], array of length n_cells.
            None means "params.permeability everywhere".
        boundary : "open"   – the aquifer continues beyond both ends, where the
                              pressure stays at its natural value p0
                   "closed" – both ends are sealed: no fluid can leave
        """
        if boundary not in ("open", "closed"):
            raise ValueError(f"boundary must be 'open' or 'closed', not {boundary!r}")

        self.params = params
        self.wells = list(wells)
        self.boundary = boundary
        self.length = length
        self.shape = (n_cells,)

        # ---- grid ----
        self.dx = length / n_cells
        self.x = cell_centres(length, n_cells)  # [m]

        if permeability is None:
            permeability = np.full(n_cells, params.permeability)
        self.permeability = np.asarray(permeability, dtype=float)
        if self.permeability.shape != self.shape:
            raise ValueError(f"permeability must have shape {self.shape}")

        flow_area = width * params.thickness  # area of every face [m²]
        cell_volume = flow_area * self.dx  # [m³]

        # ---- storage of every cell: φ · c_t · volume [m³/Pa] ----
        self.storage = np.full(n_cells, params.porosity * params.total_compressibility * cell_volume)

        # ---- transmissibility of the n-1 faces between cell i and cell i+1 ----
        k_face = harmonic_mean(self.permeability[:-1], self.permeability[1:])
        self.face_transmissibility = face_transmissibility(k_face, params.viscosity, flow_area, self.dx)

        # ---- open ends: a face to the outside, half a cell away from the centre ----
        self.boundary_transmissibility = np.zeros(n_cells)
        if boundary == "open":
            for end_cell in (0, n_cells - 1):
                self.boundary_transmissibility[end_cell] = face_transmissibility(
                    self.permeability[end_cell], params.viscosity, flow_area, self.dx / 2
                )

        # ---- which cell contains each well ----
        self.well_cells = [self.cell_index(well.x) for well in self.wells]

    def cell_index(self, x):
        """Index of the cell that contains position x [m]."""
        if not 0 <= x <= self.length:
            raise ValueError(f"x = {x} m is outside the aquifer [0, {self.length}] m")
        return min(int(x // self.dx), self.shape[0] - 1)

    def net_inflow(self, overpressure):
        """Net volume of fluid flowing INTO each cell from its neighbours [m³/s].

        This is the term ∂/∂x( k/μ ∂p/∂x ) of the equation, multiplied by the
        cell volume.
        """
        dp = overpressure

        # flow through every interior face, positive = to the right (cell i → i+1)
        flow_right = self.face_transmissibility * (dp[:-1] - dp[1:])

        inflow = np.zeros_like(dp)
        inflow[:-1] -= flow_right  # the cell left of a face loses that fluid ...
        inflow[1:] += flow_right  # ... and the cell right of it receives it

        # through an open end, fluid leaves towards the outside, where Δp = 0
        inflow -= self.boundary_transmissibility * dp
        return inflow

    def source(self, time):
        """Volume injected into each cell by the wells [m³/s] at `time` [s].

        This is the term q of the equation, multiplied by the cell volume.
        """
        q = np.zeros(self.shape)
        for well, cell in zip(self.wells, self.well_cells):
            q[cell] += well.volume_rate(time, self.params.co2_density)
        return q

    def flow_matrix(self):
        """Sparse matrix L with  L @ Δp == -net_inflow(Δp)  (for implicit solvers)."""
        n = self.shape[0]
        left_cell = np.arange(n - 1)  # face i connects cell i ...
        right_cell = left_cell + 1  # ... with cell i+1
        return assemble_flow_matrix(
            n, left_cell, right_cell, self.face_transmissibility, self.boundary_transmissibility
        )
