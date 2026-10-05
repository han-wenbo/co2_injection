"""
Automatic checks of the models and solvers.

Run from the project folder:   python -m pytest
All tests use small grids and finish in a few seconds.
"""

import numpy as np
import pytest

from co2sim.analytical import overpressure_1d, overpressure_theis
from co2sim.model_1d import Aquifer1D
from co2sim.model_2d import Aquifer2D
from co2sim.parameters import MPA, SECONDS_PER_YEAR, AquiferParameters
from co2sim.solvers import (
    SOLVERS,
    ExplicitEuler,
    make_solver,
    max_stable_time_step,
    number_of_steps,
    pick_time_step,
    run_to_end,
)
from co2sim.wells import Well

PARAMS = AquiferParameters()
YEAR = SECONDS_PER_YEAR


def random_rock(shape, seed=0):
    """Permeability that varies randomly between 10 and 1000 mD (to test non-uniform rock)."""
    rng = np.random.default_rng(seed)
    return PARAMS.permeability * 10 ** rng.uniform(-1, 1, size=shape)


def overpressure_after(model, method, dt, end_time):
    return run_to_end(model, make_solver(method, model, dt), end_time) - PARAMS.initial_pressure


# --------------------------------------------------------------------------
# The two ways of computing the flow must agree
# --------------------------------------------------------------------------
@pytest.mark.parametrize("boundary", ["open", "closed"])
def test_flow_matrix_matches_net_inflow_1d(boundary):
    model = Aquifer1D(PARAMS, 10e3, 1e3, 41, wells=[], permeability=random_rock(41), boundary=boundary)
    overpressure = np.random.default_rng(1).random(41) * MPA
    assert np.allclose(model.flow_matrix() @ overpressure, -model.net_inflow(overpressure))


@pytest.mark.parametrize("boundary", ["open", "closed"])
def test_flow_matrix_matches_net_inflow_2d(boundary):
    model = Aquifer2D(PARAMS, 10e3, 8e3, 13, 9, wells=[], permeability=random_rock((9, 13)), boundary=boundary)
    overpressure = np.random.default_rng(1).random((9, 13)) * MPA
    assert np.allclose(model.flow_matrix() @ overpressure.ravel(), -model.net_inflow(overpressure).ravel())


# --------------------------------------------------------------------------
# Physics: fluid is neither created nor destroyed
# --------------------------------------------------------------------------
@pytest.mark.parametrize("method", list(SOLVERS))
def test_closed_aquifer_keeps_all_injected_fluid(method):
    """With sealed edges nothing can leave: the extra fluid stored in all
    cells (storage · Δp) must equal the volume the well injected."""
    well = Well(x=3e3, y=4e3)
    model = Aquifer2D(PARAMS, 10e3, 10e3, 15, 15, [well], permeability=random_rock((15, 15)), boundary="closed")
    dt = pick_time_step(method, model, implicit_dt=0.01 * YEAR)
    end_time = 0.2 * YEAR

    overpressure = overpressure_after(model, method, dt, end_time)

    stored = np.sum(model.storage * overpressure)
    injected = well.volume_rate(0.0, PARAMS.co2_density) * number_of_steps(end_time, dt) * dt
    assert stored == pytest.approx(injected, rel=1e-8)


def test_pressure_returns_to_natural_after_wells_stop():
    """Open aquifer: once the well is switched off, the extra pressure leaks away."""
    well = Well(x=5e3, rate_mt_per_year=1.0, stop_year=1.0)
    model = Aquifer1D(PARAMS, 10e3, 1e3, 51, [well])
    overpressure = overpressure_after(model, "implicit", 0.01 * YEAR, 5 * YEAR)
    assert np.abs(overpressure).max() < 1e-3 * MPA


# --------------------------------------------------------------------------
# Stability of the explicit method
# --------------------------------------------------------------------------
def test_explicit_stability_limit_matches_textbook_formula():
    """Uniform rock, closed ends: dt_max = dx² / (2 D) in 1D and dx² / (4 D) in 2D with dx = dy."""
    model_1d = Aquifer1D(PARAMS, 10e3, 1e3, 50, wells=[], boundary="closed")
    model_2d = Aquifer2D(PARAMS, 10e3, 10e3, 20, 20, wells=[], boundary="closed")
    D = PARAMS.diffusivity
    assert max_stable_time_step(model_1d) == pytest.approx(model_1d.dx**2 / (2 * D))
    assert max_stable_time_step(model_2d) == pytest.approx(model_2d.dx**2 / (4 * D))


def test_explicit_euler_refuses_an_unstable_time_step():
    model = Aquifer1D(PARAMS, 10e3, 1e3, 50, wells=[])
    with pytest.raises(ValueError):
        ExplicitEuler(model, 1.01 * max_stable_time_step(model))


# --------------------------------------------------------------------------
# Comparison with exact solutions
# --------------------------------------------------------------------------
@pytest.mark.parametrize("method", list(SOLVERS))
def test_1d_well_matches_exact_solution(method):
    length, width = 100e3, 5e3
    well = Well(x=length / 2)
    model = Aquifer1D(PARAMS, length, width, 501, [well])
    dt = pick_time_step(method, model, implicit_dt=0.005 * YEAR)

    numerical = overpressure_after(model, method, dt, 1 * YEAR)
    rate = well.volume_rate(0.0, PARAMS.co2_density)
    exact = overpressure_1d(model.x - well.x, 1 * YEAR, rate, width * PARAMS.thickness, PARAMS)

    assert np.max(np.abs(numerical - exact)) < 5e-3 * exact.max()


def test_2d_well_matches_theis_solution():
    size = 30e3
    well = Well(x=size / 2, y=size / 2)
    model = Aquifer2D(PARAMS, size, size, 121, 121, [well])
    numerical = overpressure_after(model, "implicit", 0.01 * YEAR, 0.5 * YEAR)

    row, column = model.well_cells[0]
    distance = model.x[column:] - model.x[column]
    compare = (distance >= 1e3) & (distance <= 5e3)
    rate = well.volume_rate(0.0, PARAMS.co2_density)
    exact = overpressure_theis(distance[compare], 0.5 * YEAR, rate, PARAMS)

    assert numerical[row, column:][compare] == pytest.approx(exact, rel=0.02)


# --------------------------------------------------------------------------
# Linearity and symmetry
# --------------------------------------------------------------------------
def test_two_wells_add_up():
    """The equation is linear: the pressure from two wells together is the
    sum of the pressures from each well alone (superposition)."""
    well_a = Well(x=3e3, y=5e3, rate_mt_per_year=1.0)
    well_b = Well(x=7e3, y=4e3, rate_mt_per_year=0.5)
    rock = random_rock((21, 21))

    def overpressure_with(wells):
        model = Aquifer2D(PARAMS, 10e3, 10e3, 21, 21, wells, permeability=rock)
        return overpressure_after(model, "implicit", 0.02 * YEAR, 1 * YEAR)

    together = overpressure_with([well_a, well_b])
    assert together == pytest.approx(overpressure_with([well_a]) + overpressure_with([well_b]), rel=1e-9, abs=1e-6)


def test_centred_well_gives_symmetric_pressure():
    model = Aquifer2D(PARAMS, 10e3, 10e3, 21, 21, [Well(x=5e3, y=5e3)])
    overpressure = overpressure_after(model, "implicit", 0.02 * YEAR, 1 * YEAR)
    assert np.allclose(overpressure, overpressure[::-1, :])  # mirror top–bottom
    assert np.allclose(overpressure, overpressure[:, ::-1])  # mirror left–right
    assert np.allclose(overpressure, overpressure.T)  # mirror along the diagonal


def test_well_switches_on_and_off():
    well = Well(x=0.0, rate_mt_per_year=1.0, start_year=2.0, stop_year=5.0)
    assert well.volume_rate(1.9 * YEAR, 700.0) == 0.0
    assert well.volume_rate(3.0 * YEAR, 700.0) == pytest.approx(1e9 / YEAR / 700.0)
    assert well.volume_rate(5.0 * YEAR, 700.0) == 0.0
