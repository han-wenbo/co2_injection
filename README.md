# Project 1 – Underground storage of CO₂

When CO₂ is pumped into a deep saline aquifer, the pores are already full of
brine, so the pressure rises around the well and spreads outwards. If it gets
too high anywhere, the caprock can fracture. This code simulates that pressure,
in 1D and in 2D (map view), and shows it as a live animation.

It solves the equation from the project slides,

```
φ c_t ∂p/∂t = ∇·( k/μ ∇p ) + q
```

with the finite-volume method and three time-stepping schemes (explicit Euler,
backward Euler, Crank–Nicolson). This is the **sequential (single-core)**
version; the parallel version will be compared against it.

## How to run

### 1. Requirements

Python 3.10 or newer (developed with 3.12). Check your version with

```bash
python3 --version
```

### 2. Download the code

```bash
git clone git@github.com:han-wenbo/co2_injection.git
cd co2_injection
```

### 3. Install the packages (only the first time)

This creates a private Python environment in the folder `.venv` and installs
NumPy, SciPy, Matplotlib, Pillow and pytest into it.

macOS / Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Windows (PowerShell):

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Every time you open a **new terminal**, activate the environment again
(`source .venv/bin/activate` on macOS / Linux, `.venv\Scripts\Activate.ps1` on
Windows). Your prompt then starts with `(.venv)`.

### 4. Run

| command | what happens |
|---|---|
| `python run_1d.py` | opens the 1D window |
| `python run_2d.py` | opens the 2D window |
| `python run_2d.py --scenario channel` | opens the 2D window with another starting scenario |
| `python run_2d.py --save results/2d.gif` | no window: writes an animation file (`.gif` or `.mp4`) |
| `python run_1d.py --help` | lists all options |
| `python validate.py` | compares with exact solutions, saves figures to `results/` |
| `python -m pytest` | runs the automatic tests (takes about a second) |

The window opens **paused at year 0: press Play**. The terminal prints the key
numbers of the simulation (grid, time step, pressure limit) every time it
starts or restarts.

### The interactive window

| control | what it does |
|---|---|
| Play / Pause | start or stop the animation; at the end, Play replays it |
| Restart | back to year 0 |
| year slider | drag to jump to any moment (computed years are kept in memory, so going back is instant) |
| speed slider | 0.1 to 10 simulated years per second |
| injection rate, injection time, permeability, porosity | sliders; the simulation restarts with the new value |
| scenario / fault, boundary, method | radio buttons; the simulation restarts |
| map shows (2D) | pressure rise on a log scale, pressure on a linear scale, or the rock permeability |

Command-line options only set the starting values, e.g.
`run_2d.py --scenario channel` or `run_1d.py --boundary closed --fault`.
Add `--save results/name.gif` (or `.mp4`) to write an animation file instead of
opening the window. Fixed settings (grid size, area, simulated years) are at the
top of each run script.

| scenario | what you see |
|---|---|
| 1D, no fault | one well in a 50 km long aquifer strip |
| 1D, closed boundary | both ends sealed: the pressure has nowhere to go |
| 1D, sealing fault | a fault 5 km from the well traps the pressure on one side |
| 2D single | one well, 50 km × 50 km area |
| 2D fault | a sealing fault 2 km from the well |
| 2D channel | the well sits between two parallel faults (4 km apart) |
| 2D heterogeneous | rock whose permeability varies randomly |
| 2D multi-well | four wells 4 km apart: their pressures add up |

Peak pressures with the defaults (1 Mt/yr per well for 20 years, safety limit 22.95 MPa):

| case | peak pressure | safe? |
|---|---|---|
| 1D open / closed / fault | 20.8 / 23.0 / 25.0 MPa | yes / **no** (just above) / **no** |
| 2D single / fault / channel | 17.8 / 18.4 / 23.5 MPa | yes / yes / **no** |
| 2D heterogeneous / multi-well | 16.7 / 19.7 MPa | yes / yes |

Same injection rate, very different outcomes: whether the pressure can escape
(open or sealed boundaries, faults) decides whether the injection is safe.

## Project structure

```
co2sim/                 the simulation library
  parameters.py         rock, fluid and safety-limit values (all physical numbers live here)
  wells.py              wells: the source term q
  finite_volume.py      shared building blocks: transmissibility, harmonic mean, flow matrix
  model_1d.py           the aquifer as a row of cells
  model_2d.py           the aquifer as a grid of cells (map view)
  permeability.py       faults and heterogeneous rock
  solvers.py            time stepping: explicit Euler, backward Euler, Crank–Nicolson
  analytical.py         exact solutions (1D line source, Theis) for validation
  plots.py              draws the 1D and 2D pictures
  interactive.py        the interactive window (play, pause, replay, speed, parameters)
  animation.py          saves an animation file (.gif / .mp4)
  summary.py            prints the key numbers before a run
run_1d.py, run_2d.py    start the window or save an animation (fixed settings at the top)
validate.py             accuracy and convergence checks, saves figures to results/
tests/                  automatic tests (pytest)
```

## How the code works

**Finite volumes.** The aquifer is cut into cells. Each cell stores its
*overpressure* Δp = p − p0 (pressure above the natural pressure). Between two
neighbouring cells fluid flows from high to low pressure (Darcy's law):

```
flow from cell a to cell b = T_ab · (Δp_a − Δp_b)        T = (k/μ) · face area / distance
```

The face permeability is the harmonic mean of the two cells, so a thin
almost-impermeable fault blocks the flow. Each cell then obeys

```
storage · dΔp/dt = net_inflow(Δp) + source(t)        storage = φ · c_t · cell volume
```

| equation term | code |
|---|---|
| φ c_t ∂p/∂t | `model.storage` × rate of change of Δp |
| ∂/∂x( k/μ ∂p/∂x ) | `model.net_inflow(Δp)` (or `-model.flow_matrix() @ Δp`) |
| q | `model.source(t)` from the wells |

**Time stepping** (`solvers.py`):

* *explicit Euler*: `Δp_new = Δp + dt · (net_inflow + source) / storage`. One line, but
  only stable if `dt ≤ storage / (sum of face transmissibilities)` for every cell
  (`dx²/2D` in 1D, `dx²/4D` in 2D).
* *backward Euler* and *Crank–Nicolson*: solve `(S/dt + θL) Δp_new = (S/dt − (1−θ)L) Δp_old + q`
  (θ = 1 or ½) with a sparse LU factorisation that is computed once and reused.
  Stable for any `dt`.

**Boundaries.** `open`: the aquifer continues beyond the model edge, where the
pressure stays at p0. `closed`: no flow through the edge (a sealed compartment).

## Default parameters

Typical values for a deep sandstone aquifer, not one specific site (see `parameters.py`).

| symbol | meaning | value |
|---|---|---|
| φ | porosity | 0.20 |
| k | permeability | 100 mD ≈ 1e-13 m² |
| H | aquifer thickness | 100 m |
| μ | brine viscosity | 0.5 mPa·s |
| c_t | total compressibility (rock + brine) | 1e-9 1/Pa |
| p0 | natural pressure at 1500 m depth | 15.45 MPa |
| limit | 90 % of fracture pressure (17 kPa/m × 1500 m) | 22.95 MPa |
| D = k/(φμc_t) | pressure diffusivity | ≈ 1 m²/s (pressure spreads ~5.6 km in a year) |
| injection | per well, CO₂ density 700 kg/m³ | 1 Mt/yr ≈ 0.045 m³/s |

## Validation (`validate.py`)

| check | result |
|---|---|
| 1D well vs. exact solution (1, 5, 10 years) | relative error 1e-3 … 1e-4 |
| 2D well vs. Theis solution, 1–15 km from the well | relative error ≈ 1e-3 |
| refine the grid (dx 800 → 50 m) | error ∝ dx², observed order 2.00 |
| refine the time step, backward Euler | error ∝ dt, observed order 1.00 |
| refine the time step, Crank–Nicolson | error ∝ dt², observed order 2.00 |

Crank–Nicolson is the most accurate per step, but it can oscillate near the
well when `dt ≫ dx²/D`, because the well is switched on suddenly. The run
scripts therefore use backward Euler by default.

## Sequential performance (baseline for the parallel comparison)

2D default grid (201 × 201 cells), 30 simulated years, one core:

| method | time steps | run time |
|---|---|---|
| explicit Euler | 100 663 | ≈ 16 s |
| backward Euler | 600 | ≈ 1.1 s |
| Crank–Nicolson | 600 | ≈ 1.2 s |

(measured on the development laptop; run `validate.py` for the 1D numbers)

## Assumptions and limitations

* One fluid (single phase): we model the pressure, not where the CO₂ plume is.
* Constant rock and fluid properties in time, horizontal aquifer of constant thickness, no gravity.
* The value in a well's cell is the average over the cell, not the pressure
  inside the borehole (which is higher; see the Peaceman well model).
* The 1D model is Cartesian (a strip). A radial 1D model would match a single well more closely.
