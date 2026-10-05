"""
Time stepping: how to get from the pressure at time t to the pressure at t + dt.

Both models turn the PDE into one equation per cell:

    storage · dΔp/dt  =  net_inflow(Δp)  +  source(t)
     (φ c_t · V)          (flow term · V)     (q · V)

i.e. "how fast the pressure rises × how much fits in = what flows in".

Three classic ways to march this forward in time:

* ExplicitEuler  – compute the flow with the OLD pressure.
                   One line of code, but blows up if dt is too large.
* BackwardEuler  – compute the flow with the NEW pressure (implicit).
                   Stable for any dt; each step solves a sparse linear system.
* CrankNicolson  – average of old and new. Stable, and more accurate in time.

All solvers work on the overpressure Δp = p − p0 and have the same interface:
    new_overpressure = solver.step(overpressure, time)
"""

import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


def max_stable_time_step(model):
    """Largest dt for which explicit Euler stays stable [s].

    In one explicit step a cell sends  dt · T · (Δp_cell − Δp_neighbour)  to
    every neighbour. If dt is too large, a high-pressure cell gives away so
    much that it ends up LOWER than its neighbours; next step it overshoots
    back up, and the error grows every step. Every cell is safe when

        dt ≤ storage / (sum of the transmissibilities of all its faces)

    The diagonal of the flow matrix is exactly that sum.
    """
    face_sum = model.flow_matrix().diagonal()
    return float(np.min(model.storage.ravel() / face_sum))


class ExplicitEuler:
    name = "explicit Euler"

    def __init__(self, model, dt):
        dt_max = max_stable_time_step(model)
        if dt > dt_max:
            raise ValueError(
                f"dt = {dt:.4g} s is too large for explicit Euler; the stability limit is {dt_max:.4g} s"
            )
        self.model = model
        self.dt = dt

    def step(self, overpressure, time):
        # rate of change of every cell's pressure [Pa/s], using the OLD pressure
        rate = (self.model.net_inflow(overpressure) + self.model.source(time)) / self.model.storage
        return overpressure + self.dt * rate


class ThetaMethod:
    """Implicit time stepping, shared by backward Euler and Crank–Nicolson.

    With matrices,  S dΔp/dt = −L Δp + q,  where S = diag(storage) and
    L = flow matrix. Evaluate the right-hand side as a weighted average of
    the old and the new time (weight θ for the new one):

        (S/dt + θ·L) Δp_new = (S/dt − (1−θ)·L) Δp_old + θ·q_new + (1−θ)·q_old

    θ = 1   → backward Euler
    θ = 1/2 → Crank–Nicolson

    The matrix on the left never changes, so we factorise it (LU) once and
    every step only does a cheap forward/backward substitution.
    """

    def __init__(self, model, dt, theta):
        self.model = model
        self.dt = dt
        self.theta = theta

        storage_over_dt = sparse.diags(model.storage.ravel() / dt)
        flow = model.flow_matrix()
        self._lu = splu((storage_over_dt + theta * flow).tocsc())
        self._rhs_matrix = (storage_over_dt - (1 - theta) * flow).tocsr()

    def step(self, overpressure, time):
        q_old = self.model.source(time).ravel()
        q_new = self.model.source(time + self.dt).ravel()

        rhs = self._rhs_matrix @ overpressure.ravel() + self.theta * q_new + (1 - self.theta) * q_old
        return self._lu.solve(rhs).reshape(overpressure.shape)


class BackwardEuler(ThetaMethod):
    name = "backward Euler"

    def __init__(self, model, dt):
        super().__init__(model, dt, theta=1.0)


class CrankNicolson(ThetaMethod):
    name = "Crank–Nicolson"

    def __init__(self, model, dt):
        super().__init__(model, dt, theta=0.5)


SOLVERS = {
    "explicit": ExplicitEuler,
    "implicit": BackwardEuler,
    "crank-nicolson": CrankNicolson,
}


def make_solver(method, model, dt):
    """Create a solver by name: 'explicit', 'implicit' or 'crank-nicolson'."""
    if method not in SOLVERS:
        raise ValueError(f"unknown method {method!r}; choose one of {list(SOLVERS)}")
    return SOLVERS[method](model, dt)


def pick_time_step(method, model, implicit_dt):
    """Explicit Euler must stay under its stability limit; the implicit methods use `implicit_dt`."""
    if method == "explicit":
        return 0.9 * max_stable_time_step(model)
    return implicit_dt


def number_of_steps(end_time, dt):
    """Time steps needed to reach end_time (rounded to a whole number, at least 1)."""
    return max(1, round(end_time / dt))


def steps_per_snapshot(dt, snapshot_interval):
    """How many time steps lie between two saved snapshots (at least 1)."""
    return max(1, round(snapshot_interval / dt))


def count_snapshots(end_time, dt, snapshot_interval):
    """Number of (time, pressure) pairs that simulate() yields, including t = 0."""
    total = number_of_steps(end_time, dt)
    return 1 + math.ceil(total / steps_per_snapshot(dt, snapshot_interval))


def simulate(model, solver, end_time, snapshot_interval, initial_overpressure=None):
    """Run the simulation; yield (time [s], pressure [Pa]) about every `snapshot_interval` seconds.

    This is a *generator*: whoever loops over it (e.g. the animation) gets
    each snapshot as soon as it has been computed. That is what makes the
    animation run in real time.

    end_time is rounded to a whole number of time steps.
    initial_overpressure : Δp at t = 0 [Pa]. None means p = p0 everywhere
        (the aquifer before any injection).
    """
    p0 = model.params.initial_pressure
    if initial_overpressure is None:
        overpressure = np.zeros(model.shape)
    else:
        overpressure = np.array(initial_overpressure, dtype=float)

    total_steps = number_of_steps(end_time, solver.dt)
    snapshot_every = steps_per_snapshot(solver.dt, snapshot_interval)

    yield 0.0, p0 + overpressure
    for step in range(1, total_steps + 1):
        time_before_step = (step - 1) * solver.dt
        overpressure = solver.step(overpressure, time_before_step)
        if step % snapshot_every == 0 or step == total_steps:
            yield step * solver.dt, p0 + overpressure


def run_to_end(model, solver, end_time, initial_overpressure=None):
    """Run without animation; return the final pressure [Pa]."""
    *_, (_, final_pressure) = simulate(model, solver, end_time, end_time, initial_overpressure)
    return final_pressure
