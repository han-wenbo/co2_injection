"""
Drawing the pressure: one class for the 1D picture, one for the 2D picture.

Both are used by the saved animations (animation.py) and by the interactive
viewer (interactive.py), in three steps:

    plot = Plot2D(fig, ax_map, ax_history)                # create the empty picture
    plot.show_model(model, end_time)                      # draw wells, faults, axis ranges
    plot.show_snapshot(time, pressure, years, highest)    # update for one moment in time

Colours: dark = natural pressure, yellow = safety limit, red = above the limit.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap, LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from .parameters import MILLIDARCY, MPA, SECONDS_PER_YEAR

KM = 1e3  # we plot distances in kilometres


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------
def pressure_colormap():
    """viridis from low (dark blue) to the safety limit (yellow); red above the limit."""
    colormap = plt.get_cmap("viridis").copy()
    colormap.set_over("red")
    colormap.set_under(colormap(0.0))
    return colormap


def safety_status(pressure, params):
    """Text and colour that say whether the pressure is safe everywhere."""
    highest = pressure.max() / MPA
    if pressure.max() > params.max_safe_pressure:
        return f"max pressure {highest:.2f} MPa: ABOVE the safety limit!", "red"
    return f"max pressure {highest:.2f} MPa: below the safety limit", "green"


def injection_status(wells, time):
    return "wells injecting" if any(well.is_active(time) for well in wells) else "wells switched off"


def is_barrier(permeability):
    """Cells that are almost impermeable compared with the typical rock (faults)."""
    return permeability < 1e-3 * np.median(permeability)


def log_permeability_md(model):
    """log10 of the permeability in millidarcy (2 means 100 mD)."""
    return np.log10(model.permeability / MILLIDARCY)


def limits_mpa(model):
    """(natural pressure p0, safety limit) in MPa."""
    return model.params.initial_pressure / MPA, model.params.max_safe_pressure / MPA


# --------------------------------------------------------------------------
# 1D: pressure curve (top) and the aquifer as a coloured strip (bottom)
# --------------------------------------------------------------------------
class Plot1D:
    def __init__(self, fig, ax_curve, ax_strip, colorbar_ax=None):
        self.fig, self.ax_curve, self.ax_strip = fig, ax_curve, ax_strip

        (self.curve,) = ax_curve.plot([], [], color="tab:blue", lw=2, label="pressure p(x, t)")
        self.natural_line = ax_curve.axhline(0, color="gray", ls=":", label="natural pressure p0")
        self.limit_line = ax_curve.axhline(0, color="red", ls="--", label="safety limit (90 % of fracture pressure)")
        ax_curve.set_ylabel("pressure [MPa]")
        ax_curve.tick_params(labelbottom=False)

        self.strip = ax_strip.imshow(np.zeros((1, 2)), aspect="auto", cmap=pressure_colormap())
        colorbar_place = {"cax": colorbar_ax} if colorbar_ax else {"ax": [ax_curve, ax_strip], "shrink": 0.6}
        fig.colorbar(self.strip, label="pressure [MPa]", extend="max", **colorbar_place)
        ax_strip.set_yticks([])
        ax_strip.set_ylabel("aquifer")
        ax_strip.set_xlabel("position x [km]")

        self.model_artists = []  # wells and faults: redrawn whenever the model changes

    def show_model(self, model, end_time=None):
        """Draw everything that depends on the model but not on time."""
        self.model = model
        p0, limit = limits_mpa(model)
        x_km = model.x / KM
        length_km = model.length / KM

        for artist in self.model_artists:
            artist.remove()
        self.model_artists = []

        self.curve.set_data(x_km, np.full_like(x_km, p0))
        self.natural_line.set_ydata([p0, p0])
        self.limit_line.set_ydata([limit, limit])
        for well in model.wells:
            self.model_artists.append(self.ax_curve.axvline(well.x / KM, color="black", lw=0.8))
            self.model_artists.append(self.ax_curve.annotate(
                " well", (well.x / KM, 0.02), xycoords=("data", "axes fraction")))
        barrier = is_barrier(model.permeability)
        for x_barrier in model.x[barrier]:
            self.model_artists.append(self.ax_curve.axvspan(
                (x_barrier - model.dx / 2) / KM, (x_barrier + model.dx / 2) / KM, color="gray", alpha=0.4, lw=0))

        handles = [self.curve, self.natural_line, self.limit_line]
        if barrier.any():
            handles.append(Patch(color="gray", alpha=0.4, label="sealing fault"))
        self.ax_curve.legend(handles=handles, loc="upper left", fontsize=8)

        self.ax_curve.set_xlim(0, length_km)
        self.ax_curve.set_ylim(p0 - 0.5, limit + 2.5)  # room for the legend above the limit line
        self.strip.set_extent([0, length_km, 0, 1])
        self.strip.set_clim(p0, limit)
        self.ax_strip.set_xlim(0, length_km)

    def show_snapshot(self, time, pressure, history_years=None, history_highest=None):
        """Update the picture for one moment in time (the history is not used in 1D)."""
        pressure_mpa = pressure / MPA
        self.curve.set_ydata(pressure_mpa)
        self.strip.set_data(pressure_mpa[np.newaxis, :])

        bottom, top = self.ax_curve.get_ylim()
        if pressure_mpa.max() > top:  # make room if the curve leaves the plot
            self.ax_curve.set_ylim(bottom, pressure_mpa.max() + 1.0)

        text, colour = safety_status(pressure, self.model.params)
        self.ax_curve.set_title(text, color=colour, fontsize=10)
        self.fig.suptitle(f"1D aquifer  –  year {time / SECONDS_PER_YEAR:5.1f}  "
                          f"({injection_status(self.model.wells, time)})")


# --------------------------------------------------------------------------
# 2D: map of the aquifer (left) and the highest pressure over time (right)
# --------------------------------------------------------------------------
# what the map can show: the pressure rise on a log scale, the pressure on a linear scale, or the rock
MAP_RISE, MAP_PRESSURE, MAP_ROCK = "rise Δp (log)", "pressure", "permeability"
MAP_MODES = (MAP_RISE, MAP_PRESSURE, MAP_ROCK)


class Plot2D:
    def __init__(self, fig, ax_map, ax_history, colorbar_ax=None):
        self.fig, self.ax_map, self.ax_history = fig, ax_map, ax_history
        self.map_mode = MAP_RISE

        # ---- map: an image whose colours we update, plus faults and wells on top ----
        self.image = ax_map.imshow(np.ones((2, 2)), origin="lower", cmap=pressure_colormap())
        colorbar_place = {"cax": colorbar_ax} if colorbar_ax else {"ax": ax_map, "shrink": 0.85}
        self.colorbar = fig.colorbar(self.image, extend="max", **colorbar_place)
        self.barrier_overlay = ax_map.imshow(np.ma.masked_all((2, 2)), origin="lower",
                                             cmap=ListedColormap(["white"]))
        (self.well_markers,) = ax_map.plot([], [], ls="none", marker="v", markersize=10, color="white",
                                           markeredgecolor="black", label="injection well")
        ax_map.set_xlabel("x [km]")
        ax_map.set_ylabel("y [km]")

        # ---- history: highest pressure anywhere in the aquifer ----
        (self.history_line,) = ax_history.plot([], [], color="tab:blue", lw=2,
                                               label="highest pressure in the aquifer")
        self.natural_line = ax_history.axhline(0, color="gray", ls=":", label="natural pressure p0")
        self.limit_line = ax_history.axhline(0, color="red", ls="--", label="safety limit")
        self.injection_shading = None
        ax_history.set_xlabel("time [years]")
        ax_history.set_ylabel("pressure [MPa]")

    def show_model(self, model, end_time):
        """Draw everything that depends on the model but not on time."""
        self.model = model
        p0, limit = limits_mpa(model)
        extent = [0, model.size_x / KM, 0, model.size_y / KM]

        # ---- map ----
        self.image.set_extent(extent)
        barrier = is_barrier(model.permeability)
        self.barrier_overlay.set_data(np.ma.masked_where(~barrier, np.ones(model.shape)))
        self.barrier_overlay.set_extent(extent)
        self.well_markers.set_data([w.x / KM for w in model.wells], [w.y / KM for w in model.wells])
        self.ax_map.set_xlim(extent[0], extent[1])
        self.ax_map.set_ylim(extent[2], extent[3])
        handles = [self.well_markers]
        if barrier.any():
            handles.append(Line2D([], [], color="white", lw=4, label="sealing fault"))
        self.ax_map.legend(handles=handles, loc="upper right", fontsize=8, facecolor="lightgray")

        # ---- history ----
        self.history_line.set_data([], [])
        self.natural_line.set_ydata([p0, p0])
        self.limit_line.set_ydata([limit, limit])
        if self.injection_shading is not None:
            self.injection_shading.remove()
        start = min(well.start_year for well in model.wells)
        stop = min(max(well.stop_year for well in model.wells), end_time / SECONDS_PER_YEAR)
        self.injection_shading = self.ax_history.axvspan(start, stop, color="tab:orange", alpha=0.12,
                                                         label="injection period")
        self.ax_history.set_xlim(0, end_time / SECONDS_PER_YEAR)
        self.ax_history.set_ylim(p0 - 0.5, limit + 2.5)  # room for the legend above the limit line
        self.ax_history.legend(handles=[self.history_line, self.natural_line, self.limit_line,
                                        self.injection_shading], loc="upper right", fontsize=8)

        self.set_map_mode(self.map_mode)

    def set_map_mode(self, mode):
        """Choose what the map shows: one of MAP_MODES."""
        self.map_mode = mode
        p0, limit = limits_mpa(self.model)

        if mode == MAP_RISE:
            # logarithmic colours: both the big rise at the well and the small rise far away are visible
            allowed_rise = limit - p0
            self.image.set_cmap(pressure_colormap())
            self.image.set_norm(LogNorm(vmin=0.01, vmax=allowed_rise))
            self.colorbar.set_ticks([0.01, 0.1, 1, allowed_rise],
                                    labels=["0.01", "0.1", "1", f"{allowed_rise:.1f}\n(limit)"])
            self.colorbar.set_label("pressure rise Δp = p − p0 [MPa]")
            self.ax_map.set_title("pressure rise (map view)")

        elif mode == MAP_PRESSURE:
            self.image.set_cmap(pressure_colormap())
            self.image.set_norm(Normalize(vmin=p0, vmax=limit))
            ticks = list(range(int(np.ceil(p0)), int(limit) + 1, 2))
            self.colorbar.set_ticks(ticks + [limit], labels=[str(t) for t in ticks] + [f"{limit:.1f}\n(limit)"])
            self.colorbar.set_label("pressure [MPa]")
            self.ax_map.set_title("pressure (map view)")

        elif mode == MAP_ROCK:
            log_k = log_permeability_md(self.model)
            low, high = np.floor(log_k.min()), np.ceil(log_k.max())
            if low == high:  # uniform rock: give the colour scale some width
                low, high = low - 1, high + 1
            self.image.set_cmap("cividis")
            self.image.set_norm(Normalize(vmin=low, vmax=high))
            ticks = np.arange(low, high + 1)
            self.colorbar.set_ticks(ticks, labels=[f"{10 ** t:g}" for t in ticks])
            self.colorbar.set_label("permeability k [mD]")
            self.ax_map.set_title("rock permeability (map view)")
            self.image.set_data(log_k)

        else:
            raise ValueError(f"unknown map mode {mode!r}; choose one of {MAP_MODES}")

    def show_snapshot(self, time, pressure, history_years, history_highest):
        """Update the picture for one moment in time.

        history_years, history_highest : the highest pressure [MPa] at every
            earlier snapshot, for the line on the right
        """
        p0, _ = limits_mpa(self.model)
        if self.map_mode == MAP_RISE:
            rise = pressure / MPA - p0
            self.image.set_data(np.maximum(rise, 1e-6))  # a log scale cannot show 0; tiny values stay dark
        elif self.map_mode == MAP_PRESSURE:
            self.image.set_data(pressure / MPA)
        # the permeability (MAP_ROCK) does not change with time

        self.history_line.set_data(history_years, history_highest)
        bottom, top = self.ax_history.get_ylim()
        if max(history_highest) > top:
            self.ax_history.set_ylim(bottom, max(history_highest) + 1.0)

        text, colour = safety_status(pressure, self.model.params)
        self.ax_history.set_title(text, color=colour, fontsize=10)
        self.fig.suptitle(f"2D aquifer  –  year {time / SECONDS_PER_YEAR:5.1f}  "
                          f"({injection_status(self.model.wells, time)})")
