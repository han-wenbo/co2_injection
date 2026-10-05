"""
Check the solvers against exact (analytical) solutions and measure how fast
the error shrinks when the grid and the time step are refined.

    python validate.py

Prints the results and saves three figures in results/:

1. 1D: one well vs. the exact 1D solution            -> validation_1d.png
2. 2D: one well vs. the Theis solution               -> validation_2d.png
3. Convergence in space and in time                  -> convergence.png
4. Explicit vs. backward Euler vs. Crank–Nicolson: accuracy and run time (printed)
"""

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # only save figures, never open a window
import matplotlib.pyplot as plt
import numpy as np

from co2sim.analytical import overpressure_1d, overpressure_theis
from co2sim.model_1d import Aquifer1D
from co2sim.model_2d import Aquifer2D
from co2sim.parameters import MPA, SECONDS_PER_YEAR, AquiferParameters
from co2sim.solvers import make_solver, max_stable_time_step, number_of_steps, run_to_end
from co2sim.wells import Well

PARAMS = AquiferParameters()
RESULTS = Path("results")
YEAR = SECONDS_PER_YEAR

# The exact solutions assume an infinitely large aquifer, so the test aquifers
# are made large enough that their edges hardly matter during the test.
LENGTH_1D = 200e3  # [m]
WIDTH_1D = 5e3  # [m]
SIZE_2D = 60e3  # [m]


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def relative_error(numerical, exact):
    """Largest difference between the two, relative to the largest exact value."""
    return np.max(np.abs(numerical - exact)) / np.max(np.abs(exact))


def observed_orders(sizes, errors):
    """Convergence order between successive refinements: log(e1/e2) / log(h1/h2)."""
    return [np.log(errors[i] / errors[i + 1]) / np.log(sizes[i] / sizes[i + 1]) for i in range(len(errors) - 1)]


def one_well_1d(n_cells):
    well = Well(x=LENGTH_1D / 2)
    return Aquifer1D(PARAMS, LENGTH_1D, WIDTH_1D, n_cells, [well]), well


def exact_1d(model, well, time):
    rate = well.volume_rate(0.0, PARAMS.co2_density)
    return overpressure_1d(model.x - well.x, time, rate, WIDTH_1D * PARAMS.thickness, PARAMS)


def numerical_overpressure(model, method, dt, end_time, initial_overpressure=None):
    solver = make_solver(method, model, dt)
    return run_to_end(model, solver, end_time, initial_overpressure) - PARAMS.initial_pressure


# --------------------------------------------------------------------------
# 1. 1D well vs. exact solution
# --------------------------------------------------------------------------
def check_1d_against_exact():
    print("\n1) 1D well vs. exact solution (1001 cells, backward Euler, dt = 0.01 year)")
    model, well = one_well_1d(n_cells=1001)
    fig, ax = plt.subplots(figsize=(8, 4.5))

    for years, colour in [(1, "tab:blue"), (5, "tab:orange"), (10, "tab:green")]:
        numerical = numerical_overpressure(model, "implicit", 0.01 * YEAR, years * YEAR)
        exact = exact_1d(model, well, years * YEAR)
        print(f"   after {years:2d} years: max overpressure {exact.max() / MPA:5.2f} MPa, "
              f"relative error {relative_error(numerical, exact):.1e}")

        x_km = (model.x - well.x) / 1e3
        ax.plot(x_km, exact / MPA, color=colour, label=f"exact, {years} yr")
        ax.plot(x_km[::20], numerical[::20] / MPA, "o", color=colour, ms=4, label=f"numerical, {years} yr")

    ax.set_xlabel("distance from the well [km]")
    ax.set_ylabel("overpressure Δp [MPa]")
    ax.set_title("1D: numerical solution (dots) vs. exact solution (lines)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "validation_1d.png", dpi=120)


# --------------------------------------------------------------------------
# 2. 2D well vs. Theis solution
# --------------------------------------------------------------------------
def check_2d_against_theis():
    print("\n2) 2D well vs. Theis solution (241 × 241 cells, backward Euler, dt = 0.01 year)")
    centre = SIZE_2D / 2
    well = Well(x=centre, y=centre)
    model = Aquifer2D(PARAMS, SIZE_2D, SIZE_2D, 241, 241, [well])
    rate = well.volume_rate(0.0, PARAMS.co2_density)

    # compare along the line from the well towards +x
    row, column = model.well_cells[0]
    distance = model.x[column + 1:] - model.x[column]
    away_from_well = (distance >= 1e3) & (distance <= 15e3)  # the well cell itself is a 250 m average

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for years, colour in [(1, "tab:blue"), (3, "tab:orange")]:
        overpressure = numerical_overpressure(model, "implicit", 0.01 * YEAR, years * YEAR)
        numerical = overpressure[row, column + 1:]
        exact = overpressure_theis(distance, years * YEAR, rate, PARAMS)
        error = relative_error(numerical[away_from_well], exact[away_from_well])
        print(f"   after {years} years, 1–15 km from the well: relative error {error:.1e}")

        ax.plot(distance / 1e3, exact / MPA, color=colour, label=f"Theis (exact), {years} yr")
        ax.plot(distance[::6] / 1e3, numerical[::6] / MPA, "o", color=colour, ms=4,
                label=f"numerical, {years} yr")

    ax.set_xscale("log")
    ax.set_xlabel("distance from the well [km]")
    ax.set_ylabel("overpressure Δp [MPa]")
    ax.set_title("2D: numerical solution (dots) vs. Theis solution (lines)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "validation_2d.png", dpi=120)


# --------------------------------------------------------------------------
# 3. Convergence studies
# --------------------------------------------------------------------------
def convergence_in_space():
    """Refine the grid; the time step is tiny so that only the space error is left.
    Expected: error ∝ dx², i.e. order 2."""
    print("\n3a) convergence in space (1D well after 1 year, Crank–Nicolson, dt = 0.001 year)")
    cell_counts = [251, 501, 1001, 2001, 4001]
    spacings, errors = [], []
    for n_cells in cell_counts:
        model, well = one_well_1d(n_cells)
        numerical = numerical_overpressure(model, "crank-nicolson", 0.001 * YEAR, 1 * YEAR)
        spacings.append(model.dx)
        errors.append(relative_error(numerical, exact_1d(model, well, 1 * YEAR)))

    orders = observed_orders(spacings, errors)
    for i, (dx, error) in enumerate(zip(spacings, errors)):
        order = f"   order {orders[i - 1]:.2f}" if i > 0 else ""
        print(f"   dx = {dx:6.1f} m: relative error {error:.2e}{order}")
    return spacings, errors


def gaussian_bump(x, time, centre, height, width):
    """Exact solution for an initial Gaussian pressure bump with no wells: it
    spreads out and gets lower, but its total volume stays the same."""
    spread_sq = width**2 + 2 * PARAMS.diffusivity * time
    return height * width / np.sqrt(spread_sq) * np.exp(-((x - centre) ** 2) / (2 * spread_sq))


def convergence_in_time():
    """Refine the time step on a fine grid. We use a smooth initial bump
    instead of a well: a well switched on suddenly is not smooth, which
    spoils the accuracy of Crank–Nicolson (see 4).
    Expected: backward Euler error ∝ dt (order 1), Crank–Nicolson ∝ dt² (order 2)."""
    print("\n3b) convergence in time (smooth pressure bump after 1 year, 8001 cells)")
    model = Aquifer1D(PARAMS, LENGTH_1D, WIDTH_1D, 8001, wells=[])
    bump = dict(centre=LENGTH_1D / 2, height=1 * MPA, width=5e3)
    exact = gaussian_bump(model.x, 1 * YEAR, **bump)

    time_steps_years = [0.2, 0.1, 0.05, 0.025]
    errors = {}
    for method in ["implicit", "crank-nicolson"]:
        errors[method] = []
        for dt_years in time_steps_years:
            start = gaussian_bump(model.x, 0.0, **bump)
            numerical = numerical_overpressure(model, method, dt_years * YEAR, 1 * YEAR, start)
            errors[method].append(relative_error(numerical, exact))
        orders = observed_orders(time_steps_years, errors[method])
        print(f"   {method:15s} errors " + "  ".join(f"{e:.1e}" for e in errors[method])
              + "   orders " + "  ".join(f"{o:.2f}" for o in orders))
    return time_steps_years, errors


def plot_convergence(spacings, space_errors, time_steps, time_errors):
    fig, (ax_space, ax_time) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax_space.loglog(spacings, space_errors, "o-", label="numerical")
    reference = space_errors[0] * (np.array(spacings) / spacings[0]) ** 2
    ax_space.loglog(spacings, reference, "k--", label="slope 2 (error ∝ dx²)")
    ax_space.set_xticks(spacings, [f"{dx:.0f}" for dx in spacings])
    ax_space.minorticks_off()
    ax_space.set_xlabel("cell size dx [m]")
    ax_space.set_ylabel("relative error")
    ax_space.set_title("convergence in space")
    ax_space.legend()

    for method, marker in [("implicit", "o-"), ("crank-nicolson", "s-")]:
        ax_time.loglog(time_steps, time_errors[method], marker, label=method)
    steps = np.array(time_steps) / time_steps[0]
    ax_time.loglog(time_steps, time_errors["implicit"][0] * steps, "k--", label="slope 1")
    ax_time.loglog(time_steps, time_errors["crank-nicolson"][0] * steps**2, "k:", label="slope 2")
    ax_time.set_xticks(time_steps, [f"{dt:g}" for dt in time_steps])
    ax_time.minorticks_off()
    ax_time.set_xlabel("time step dt [years]")
    ax_time.set_title("convergence in time")
    ax_time.legend()

    fig.tight_layout()
    fig.savefig(RESULTS / "convergence.png", dpi=120)


# --------------------------------------------------------------------------
# 4. Explicit vs. implicit vs. Crank–Nicolson
# --------------------------------------------------------------------------
def compare_methods():
    print("\n4) the three methods on the 1D well problem (1001 cells, 5 years)")
    model, well = one_well_1d(n_cells=1001)
    exact = exact_1d(model, well, 5 * YEAR)
    dt_explicit = 0.9 * max_stable_time_step(model)

    runs = [
        ("explicit", dt_explicit),
        ("implicit", 0.05 * YEAR),
        ("implicit", 0.005 * YEAR),
        ("crank-nicolson", 0.05 * YEAR),
        ("crank-nicolson", 0.005 * YEAR),
    ]
    print(f"   {'method':15s} {'dt [days]':>10s} {'steps':>7s} {'run time [s]':>13s} {'rel. error':>11s}")
    for method, dt in runs:
        start = time.perf_counter()
        numerical = numerical_overpressure(model, method, dt, 5 * YEAR)
        seconds = time.perf_counter() - start
        steps = number_of_steps(5 * YEAR, dt)
        print(f"   {method:15s} {dt / 86400:10.2f} {steps:7d} {seconds:13.3f} "
              f"{relative_error(numerical, exact):11.1e}")
    print("   explicit: accurate, but needs thousands of tiny steps to stay stable.\n"
          "   implicit methods: 100 large steps are enough. For the same dt, Crank–Nicolson\n"
          "   is more accurate than backward Euler, but it can oscillate near the well when\n"
          "   dt is much larger than dx²/D (the well is switched on suddenly).")


def main():
    RESULTS.mkdir(exist_ok=True)
    check_1d_against_exact()
    check_2d_against_theis()
    spacings, space_errors = convergence_in_space()
    time_steps, time_errors = convergence_in_time()
    plot_convergence(spacings, space_errors, time_steps, time_errors)
    compare_methods()
    print(f"\nfigures saved in {RESULTS}/")


if __name__ == "__main__":
    main()
