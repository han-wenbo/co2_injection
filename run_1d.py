"""
1D simulation of CO2 injection into a saline aquifer.

Run from the project folder:

    python run_1d.py                            interactive window: play, pause, replay,
                                                change the speed and the parameters
    python run_1d.py --boundary closed --fault  start the window with other settings
    python run_1d.py --save results/1d.gif      save an animation file instead (no window)

The command-line options only set the starting values; in the window
everything can be changed.
"""

import argparse
from dataclasses import dataclass, replace

import numpy as np

from co2sim.animation import save_animation_1d
from co2sim.finite_volume import cell_centres
from co2sim.interactive import ChoiceControl, NumberControl, open_viewer_1d
from co2sim.model_1d import Aquifer1D
from co2sim.parameters import MILLIDARCY, SECONDS_PER_YEAR, AquiferParameters
from co2sim.permeability import add_fault_1d
from co2sim.solvers import SOLVERS, count_snapshots, make_solver, pick_time_step, simulate
from co2sim.summary import print_summary
from co2sim.wells import Well

# --------------------------------------------------------------------------
# Fixed settings (edit here)
# --------------------------------------------------------------------------
LENGTH = 50e3  # aquifer length [m]
WIDTH = 5e3  # aquifer width [m]; in 1D the fluid can only flow along x
N_CELLS = 251  # odd, so the well sits exactly in the middle cell
SIMULATED_YEARS = 40

FAULT_DISTANCE = 5e3  # distance from the well to the fault [m]
FAULT_WIDTH = 400.0  # [m]

IMPLICIT_DT_YEARS = 0.05  # time step of the implicit methods [years] (~18 days)
YEARS_PER_FRAME = 0.25  # saved animations: one frame every 3 months


# --------------------------------------------------------------------------
# Settings that can be changed in the window
# --------------------------------------------------------------------------
@dataclass
class Settings:
    fault: str = "none"  # "none" or "sealing fault" (5 km right of the well)
    boundary: str = "open"
    method: str = "implicit"
    rate: float = 1.0  # injection rate [Mt CO2 per year] (Sleipner, Norway: about 1)
    injection_years: float = 20  # the well is switched off after this many years
    permeability_md: float = 100.0  # permeability of the rock [mD]
    porosity: float = 0.20


NUMBER_CONTROLS = [
    NumberControl("rate", "injection rate [Mt/yr]", 0.1, 5.0, step=0.1),
    NumberControl("injection_years", "injection time [years]", 1, 35, step=1),
    NumberControl("permeability_md", "permeability k [mD]", 10, 1000, log_scale=True),
    NumberControl("porosity", "porosity φ", 0.05, 0.35, step=0.01),
]
CHOICE_CONTROLS = [
    ChoiceControl("fault", "fault", ("none", "sealing fault")),
    ChoiceControl("boundary", "boundary", ("open", "closed")),
    ChoiceControl("method", "method", tuple(SOLVERS)),
]


def build_simulation(settings):
    """Create the model and the solver for the given settings: (model, solver, end_time)."""
    params = replace(AquiferParameters(), permeability=settings.permeability_md * MILLIDARCY,
                     porosity=settings.porosity)
    well = Well(x=LENGTH / 2, rate_mt_per_year=settings.rate, stop_year=settings.injection_years)

    permeability = np.full(N_CELLS, params.permeability)  # the same rock everywhere ...
    if settings.fault == "sealing fault":  # ... except for an (optional) sealing fault
        permeability = add_fault_1d(permeability, cell_centres(LENGTH, N_CELLS),
                                    fault_x=well.x + FAULT_DISTANCE, fault_width=FAULT_WIDTH)

    aquifer = Aquifer1D(params, LENGTH, WIDTH, N_CELLS, [well], permeability, boundary=settings.boundary)
    dt = pick_time_step(settings.method, aquifer, IMPLICIT_DT_YEARS * SECONDS_PER_YEAR)
    solver = make_solver(settings.method, aquifer, dt)
    end_time = SIMULATED_YEARS * SECONDS_PER_YEAR

    print_summary(aquifer, solver, end_time)
    return aquifer, solver, end_time


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--boundary", choices=["open", "closed"], default="open",
                        help="open: aquifer continues beyond the ends; closed: ends are sealed")
    parser.add_argument("--fault", action="store_true", help="add a sealing fault 5 km right of the well")
    parser.add_argument("--method", choices=list(SOLVERS), default="implicit", help="time-stepping method")
    parser.add_argument("--rate", type=float, default=1.0, help="injection rate [Mt CO2 per year]")
    parser.add_argument("--save", metavar="FILE", help="save an animation (.gif or .mp4) instead of opening the window")
    return parser.parse_args()


def main():
    args = parse_arguments()
    settings = Settings(fault="sealing fault" if args.fault else "none", boundary=args.boundary,
                        method=args.method, rate=args.rate)

    if args.save:
        aquifer, solver, end_time = build_simulation(settings)
        frame_interval = YEARS_PER_FRAME * SECONDS_PER_YEAR
        save_animation_1d(
            aquifer,
            make_snapshots=lambda: simulate(aquifer, solver, end_time, frame_interval),
            n_frames=count_snapshots(end_time, solver.dt, frame_interval),
            path=args.save,
            end_time=end_time,
        )
    else:
        open_viewer_1d(settings, build_simulation, NUMBER_CONTROLS, CHOICE_CONTROLS)


if __name__ == "__main__":
    main()
