"""
2D (map view) simulation of CO2 injection into a saline aquifer.

Run from the project folder:

    python run_2d.py                            interactive window: play, pause, replay,
                                                change the speed and the parameters
    python run_2d.py --scenario channel         start the window with another scenario
    python run_2d.py --save results/2d.gif      save an animation file instead (no window)

Scenarios
    single          one well in the centre
    fault           one well next to a sealing fault
    channel         one well between two parallel sealing faults
    heterogeneous   one well in rock whose permeability varies
    multi-well      four wells close together

The command-line options only set the starting values; in the window
everything can be changed.
"""

import argparse
from dataclasses import dataclass, replace

import numpy as np

from co2sim.animation import save_animation_2d
from co2sim.finite_volume import cell_centres
from co2sim.interactive import ChoiceControl, NumberControl, open_viewer_2d
from co2sim.model_2d import Aquifer2D
from co2sim.parameters import MILLIDARCY, SECONDS_PER_YEAR, AquiferParameters
from co2sim.permeability import add_fault_2d, random_heterogeneous_2d
from co2sim.solvers import SOLVERS, count_snapshots, make_solver, pick_time_step, simulate
from co2sim.summary import print_summary
from co2sim.wells import Well

# --------------------------------------------------------------------------
# Fixed settings (edit here)
# --------------------------------------------------------------------------
SIZE = 50e3  # the modelled area is SIZE × SIZE [m]
N_CELLS = 201  # cells per side (odd, so a centred well sits in the middle of a cell)
SIMULATED_YEARS = 40

FAULT_DISTANCE = 2e3  # "fault": distance from the well to the fault [m]
FAULT_LENGTH = 30e3  # [m]; fluid can flow around the ends of the fault
CHANNEL_WIDTH = 4e3  # "channel": distance between the two parallel faults [m]
WELL_SPACING = 4e3  # "multi-well": distance between neighbouring wells [m]

IMPLICIT_DT_YEARS = 0.05  # time step of the implicit methods [years] (~18 days)
YEARS_PER_FRAME = 0.25  # saved animations: one frame every 3 months

SCENARIOS = ("single", "fault", "channel", "heterogeneous", "multi-well")


# --------------------------------------------------------------------------
# Settings that can be changed in the window
# --------------------------------------------------------------------------
@dataclass
class Settings:
    scenario: str = "single"
    boundary: str = "open"
    method: str = "implicit"
    rate: float = 1.0  # injection rate per well [Mt CO2 per year]
    injection_years: float = 20  # the wells are switched off after this many years
    permeability_md: float = 100.0  # permeability of the rock [mD]
    porosity: float = 0.20


NUMBER_CONTROLS = [
    NumberControl("rate", "injection rate [Mt/yr]", 0.1, 5.0, step=0.1),
    NumberControl("injection_years", "injection time [years]", 1, 35, step=1),
    NumberControl("permeability_md", "permeability k [mD]", 10, 1000, log_scale=True),
    NumberControl("porosity", "porosity φ", 0.05, 0.35, step=0.01),
]
CHOICE_CONTROLS = [
    ChoiceControl("scenario", "scenario", SCENARIOS),
    ChoiceControl("boundary", "boundary", ("open", "closed")),
    ChoiceControl("method", "method", tuple(SOLVERS)),
]


def build_scenario(name, params, rate, injection_years):
    """Return (wells, permeability) for the chosen scenario."""
    centre = SIZE / 2
    x = y = cell_centres(SIZE, N_CELLS)
    uniform_rock = np.full((N_CELLS, N_CELLS), params.permeability)

    def injector(x_well, y_well):
        return Well(x=x_well, y=y_well, rate_mt_per_year=rate, stop_year=injection_years)

    if name == "single":
        return [injector(centre, centre)], uniform_rock

    if name == "fault":
        permeability = add_fault_2d(uniform_rock, x, y, fault_x=centre + FAULT_DISTANCE,
                                    fault_y_start=centre - FAULT_LENGTH / 2,
                                    fault_y_end=centre + FAULT_LENGTH / 2)
        return [injector(centre, centre)], permeability

    if name == "channel":
        # two faults over the full height: the fluid can only escape north or south
        permeability = add_fault_2d(uniform_rock, x, y, centre - CHANNEL_WIDTH / 2, 0, SIZE)
        permeability = add_fault_2d(permeability, x, y, centre + CHANNEL_WIDTH / 2, 0, SIZE)
        return [injector(centre, centre)], permeability

    if name == "heterogeneous":
        permeability = random_heterogeneous_2d((N_CELLS, N_CELLS), mean=params.permeability)
        return [injector(centre, centre)], permeability

    if name == "multi-well":
        half = WELL_SPACING / 2
        corners = [(-half, -half), (half, -half), (-half, half), (half, half)]
        return [injector(centre + dx, centre + dy) for dx, dy in corners], uniform_rock

    raise ValueError(f"unknown scenario {name!r}")


def build_simulation(settings):
    """Create the model and the solver for the given settings: (model, solver, end_time)."""
    params = replace(AquiferParameters(), permeability=settings.permeability_md * MILLIDARCY,
                     porosity=settings.porosity)
    wells, permeability = build_scenario(settings.scenario, params, settings.rate, settings.injection_years)
    aquifer = Aquifer2D(params, SIZE, SIZE, N_CELLS, N_CELLS, wells, permeability, boundary=settings.boundary)

    dt = pick_time_step(settings.method, aquifer, IMPLICIT_DT_YEARS * SECONDS_PER_YEAR)
    solver = make_solver(settings.method, aquifer, dt)
    end_time = SIMULATED_YEARS * SECONDS_PER_YEAR

    print(f"\nscenario: {settings.scenario}")
    print_summary(aquifer, solver, end_time)
    return aquifer, solver, end_time


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenario", choices=SCENARIOS, default="single", help="which situation to simulate")
    parser.add_argument("--boundary", choices=["open", "closed"], default="open",
                        help="open: aquifer continues beyond the edges; closed: edges are sealed")
    parser.add_argument("--method", choices=list(SOLVERS), default="implicit", help="time-stepping method")
    parser.add_argument("--rate", type=float, default=1.0, help="injection rate per well [Mt CO2 per year]")
    parser.add_argument("--save", metavar="FILE", help="save an animation (.gif or .mp4) instead of opening the window")
    return parser.parse_args()


def main():
    args = parse_arguments()
    settings = Settings(scenario=args.scenario, boundary=args.boundary, method=args.method, rate=args.rate)

    if args.save:
        aquifer, solver, end_time = build_simulation(settings)
        frame_interval = YEARS_PER_FRAME * SECONDS_PER_YEAR
        save_animation_2d(
            aquifer,
            make_snapshots=lambda: simulate(aquifer, solver, end_time, frame_interval),
            n_frames=count_snapshots(end_time, solver.dt, frame_interval),
            path=args.save,
            end_time=end_time,
        )
    else:
        open_viewer_2d(settings, build_simulation, NUMBER_CONTROLS, CHOICE_CONTROLS)


if __name__ == "__main__":
    main()
