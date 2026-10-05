"""
2D model: the aquifer seen from above (map view), cut into nx × ny cells.

    y ↑
      | (0,ny-1)  ...  (nx-1,ny-1) |
      |   ...      ★ well    ...   |
      | (0,0)     ...    (nx-1,0)  |
      +----------------------------→ x

Arrays have shape (ny, nx): `array[j, i]` is the cell at (x[i], y[j]).
Row j is the y-direction, column i the x-direction (the layout that
matplotlib's imshow uses).

Compared with 1D, every cell now has up to four neighbours (left, right,
below, above), so fluid can spread out in all directions.
"""

import numpy as np

from .finite_volume import assemble_flow_matrix, cell_centres, face_transmissibility, harmonic_mean


class Aquifer2D:
    def __init__(self, params, size_x, size_y, nx, ny, wells, permeability=None, boundary="open"):
        """
        Parameters
        ----------
        params : AquiferParameters
        size_x, size_y : size of the modelled area [m]
        nx, ny : number of cells in x and y
        wells : list of Well
        permeability : k of every cell [m²], array of shape (ny, nx).
            None means "params.permeability everywhere".
        boundary : "open"   – the aquifer continues beyond the edges, where the
                              pressure stays at its natural value p0
                   "closed" – all four edges are sealed: no fluid can leave
        """
        if boundary not in ("open", "closed"):
            raise ValueError(f"boundary must be 'open' or 'closed', not {boundary!r}")

        self.params = params
        self.wells = list(wells)
        self.boundary = boundary
        self.size_x, self.size_y = size_x, size_y
        self.shape = (ny, nx)

        # ---- grid ----
        self.dx = size_x / nx
        self.dy = size_y / ny
        self.x = cell_centres(size_x, nx)  # [m]
        self.y = cell_centres(size_y, ny)

        if permeability is None:
            permeability = np.full(self.shape, params.permeability)
        k = self.permeability = np.asarray(permeability, dtype=float)
        if k.shape != self.shape:
            raise ValueError(f"permeability must have shape {self.shape}")

        H, mu = params.thickness, params.viscosity
        cell_volume = self.dx * self.dy * H  # [m³]

        # ---- storage of every cell: φ · c_t · volume [m³/Pa] ----
        self.storage = np.full(self.shape, params.porosity * params.total_compressibility * cell_volume)

        # ---- faces between left/right neighbours: shape (ny, nx-1) ----
        # fluid crosses an area dy × H and travels a distance dx
        k_x = harmonic_mean(k[:, :-1], k[:, 1:])
        self.transmissibility_x = face_transmissibility(k_x, mu, self.dy * H, self.dx)

        # ---- faces between lower/upper neighbours: shape (ny-1, nx) ----
        # fluid crosses an area dx × H and travels a distance dy
        k_y = harmonic_mean(k[:-1, :], k[1:, :])
        self.transmissibility_y = face_transmissibility(k_y, mu, self.dx * H, self.dy)

        # ---- open edges: a face to the outside, half a cell away ----
        self.boundary_transmissibility = np.zeros(self.shape)
        if boundary == "open":
            to_left_right = face_transmissibility(k, mu, self.dy * H, self.dx / 2)
            to_bottom_top = face_transmissibility(k, mu, self.dx * H, self.dy / 2)
            self.boundary_transmissibility[:, 0] += to_left_right[:, 0]  # left edge
            self.boundary_transmissibility[:, -1] += to_left_right[:, -1]  # right edge
            self.boundary_transmissibility[0, :] += to_bottom_top[0, :]  # bottom edge
            self.boundary_transmissibility[-1, :] += to_bottom_top[-1, :]  # top edge

        # ---- which cell contains each well ----
        self.well_cells = [self.cell_index(well.x, well.y) for well in self.wells]

    def cell_index(self, x, y):
        """(row j, column i) of the cell that contains the point (x, y) [m]."""
        if not (0 <= x <= self.size_x and 0 <= y <= self.size_y):
            raise ValueError(f"point ({x}, {y}) m is outside the modelled area")
        ny, nx = self.shape
        i = min(int(x // self.dx), nx - 1)
        j = min(int(y // self.dy), ny - 1)
        return j, i

    def net_inflow(self, overpressure):
        """Net volume of fluid flowing INTO each cell from its neighbours [m³/s].

        2D version of the term ∂/∂x( k/μ ∂p/∂x ) + ∂/∂y( k/μ ∂p/∂y ),
        multiplied by the cell volume.
        """
        dp = overpressure

        # flow through the faces, positive = towards +x (right) or +y (up)
        flow_right = self.transmissibility_x * (dp[:, :-1] - dp[:, 1:])
        flow_up = self.transmissibility_y * (dp[:-1, :] - dp[1:, :])

        inflow = np.zeros_like(dp)
        inflow[:, :-1] -= flow_right  # cell left of a face loses ...
        inflow[:, 1:] += flow_right  # ... cell right of it gains
        inflow[:-1, :] -= flow_up  # cell below a face loses ...
        inflow[1:, :] += flow_up  # ... cell above it gains

        # through open edges, fluid leaves towards the outside, where Δp = 0
        inflow -= self.boundary_transmissibility * dp
        return inflow

    def source(self, time):
        """Volume injected into each cell by the wells [m³/s] at `time` [s] (the q term)."""
        q = np.zeros(self.shape)
        for well, cell in zip(self.wells, self.well_cells):
            q[cell] += well.volume_rate(time, self.params.co2_density)
        return q

    def flow_matrix(self):
        """Sparse matrix L with  L @ Δp.ravel() == -net_inflow(Δp).ravel()."""
        ny, nx = self.shape
        index = np.arange(nx * ny).reshape(self.shape)  # cell (j, i) has number index[j, i]

        # every face connects cell_a with cell_b
        cell_a = np.concatenate([index[:, :-1].ravel(), index[:-1, :].ravel()])
        cell_b = np.concatenate([index[:, 1:].ravel(), index[1:, :].ravel()])
        transmissibility = np.concatenate([self.transmissibility_x.ravel(), self.transmissibility_y.ravel()])

        return assemble_flow_matrix(
            nx * ny, cell_a, cell_b, transmissibility, self.boundary_transmissibility.ravel()
        )
