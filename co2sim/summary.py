"""Print the key numbers of a simulation before it starts."""

import math

from .parameters import MPA, MILLIDARCY, SECONDS_PER_YEAR
from .solvers import max_stable_time_step


def print_summary(model, solver, end_time):
    p = model.params
    if len(model.shape) == 1:
        grid = f"{model.shape[0]} cells, dx = {model.dx:.0f} m"
    else:
        ny, nx = model.shape
        grid = f"{nx} × {ny} cells, dx = {model.dx:.0f} m, dy = {model.dy:.0f} m"
    total_rate = sum(w.rate_mt_per_year for w in model.wells)

    print("-" * 66)
    print(f"grid                 : {grid}, {model.boundary} boundary")
    print(f"rock                 : φ = {p.porosity}, k = {p.permeability / MILLIDARCY:.0f} mD, H = {p.thickness:.0f} m")
    print(f"fluid                : μ = {p.viscosity * 1e3:.2f} mPa·s, c_t = {p.total_compressibility:.1e} 1/Pa")
    print(f"wells                : {len(model.wells)}, total {total_rate:.2f} Mt CO2 per year")
    print(f"natural pressure p0  : {p.initial_pressure / MPA:.2f} MPa (at {p.depth:.0f} m depth)")
    print(f"safety limit         : {p.max_safe_pressure / MPA:.2f} MPa "
          f"({p.safety_factor:.0%} of fracture pressure {p.fracture_pressure / MPA:.1f} MPa)")
    print(f"diffusivity D        : {p.diffusivity:.2f} m²/s -> pressure spreads "
          f"~{math.sqrt(p.diffusivity * SECONDS_PER_YEAR) / 1e3:.1f} km in 1 year")
    print(f"time stepping        : {solver.name}, dt = {solver.dt / 86400:.2f} days "
          f"(explicit stability limit: {max_stable_time_step(model) / 86400:.3f} days)")
    print(f"simulated time       : {end_time / SECONDS_PER_YEAR:.0f} years "
          f"= {math.ceil(end_time / solver.dt)} time steps")
    print("-" * 66)
