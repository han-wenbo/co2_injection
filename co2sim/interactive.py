"""
Interactive viewer: watch the simulation, pause, replay, jump to any year,
change the playback speed, and change the parameters while it runs.

Controls at the bottom of the window
    Play / Pause        start or stop the animation (at the end, Play replays it)
    Restart             go back to year 0
    year slider         drag to jump to any moment
    speed slider        simulated years per second of real time
    parameter sliders   the simulation restarts with the new values
    and buttons

How it works
    * SnapshotCache runs the simulation only as far as needed and keeps every
      snapshot, so replaying or jumping back costs nothing.
    * A timer fires 25 times per second. Each time, it moves the displayed year
      forward according to the speed and draws the matching snapshot.
    * When a parameter changes, build_simulation(settings) from the run script
      creates a new model and solver, and a new cache starts at year 0.
"""

import time as clock
from dataclasses import dataclass, replace

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button, RadioButtons, Slider

from .parameters import MPA, SECONDS_PER_YEAR
from .plots import MAP_MODES, Plot1D, Plot2D
from .solvers import count_snapshots, number_of_steps, simulate, steps_per_snapshot


# --------------------------------------------------------------------------
# Descriptions of the controls (the run scripts list which ones they want)
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class NumberControl:
    """A slider that sets one number in the settings."""

    field: str  # name of the field in the settings dataclass
    label: str
    minimum: float
    maximum: float
    step: float | None = None  # None = any value
    log_scale: bool = False  # the slider moves in powers of ten (good for permeability)


@dataclass(frozen=True)
class ChoiceControl:
    """A group of radio buttons that sets one option in the settings."""

    field: str
    label: str
    options: tuple


def format_number(value):
    """Short text for a slider value: 1000 -> '1000', 2.5 -> '2.5', 0.1 -> '0.1'."""
    return f"{value:.0f}" if value >= 10 else f"{value:.2g}"


# --------------------------------------------------------------------------
# Memory of everything computed so far
# --------------------------------------------------------------------------
class SnapshotCache:
    """Runs a simulation only as far as needed and remembers every snapshot."""

    def __init__(self, model, solver, end_time, snapshot_interval):
        self._snapshots = simulate(model, solver, end_time, snapshot_interval)
        self.count = count_snapshots(end_time, solver.dt, snapshot_interval)
        self.years_between = steps_per_snapshot(solver.dt, snapshot_interval) * solver.dt / SECONDS_PER_YEAR
        self.last_year = number_of_steps(end_time, solver.dt) * solver.dt / SECONDS_PER_YEAR

        self.times = []  # [s]
        self.pressures = []  # [Pa]
        self.years = []  # same as times, in years (for the history line)
        self.highest = []  # highest pressure in the aquifer at each snapshot [MPa]

    def index_at(self, year):
        """Number of the last snapshot at or before `year`."""
        return min(int(year / self.years_between + 1e-9), self.count - 1)

    def get(self, index):
        """Snapshot number `index` as (time, pressure); computes the missing ones first."""
        while len(self.times) <= index:
            time, pressure = next(self._snapshots)
            self.times.append(time)
            self.pressures.append(pressure.astype(np.float32))  # half the memory, precise enough to draw
            self.years.append(time / SECONDS_PER_YEAR)
            self.highest.append(pressure.max() / MPA)
        return self.times[index], self.pressures[index]


# --------------------------------------------------------------------------
# The viewer
# --------------------------------------------------------------------------
class InteractiveViewer:
    TICK_SECONDS = 0.04  # redraw 25 times per second
    SNAPSHOT_YEARS = 0.1  # keep one snapshot every 0.1 simulated years
    REBUILD_DELAY = 0.3  # [s] wait until a slider has stopped moving before restarting

    def __init__(self, fig, plot, settings, build_simulation, number_controls, choice_controls, map_modes=None):
        """
        fig, plot : the figure and the Plot1D / Plot2D drawn in its upper part
        settings : dataclass with the starting values of all controls
        build_simulation : function settings -> (model, solver, end_time)
        number_controls, choice_controls : which sliders and radio buttons to show
        map_modes : 2D only, the options of the "map shows" buttons
        """
        self.fig = fig
        self.plot = plot
        self.settings = settings
        self.build_simulation = build_simulation

        self.playing = False
        self.speed = 2.0  # simulated years per second of real time
        self.display_year = 0.0
        self.rebuild_requested_at = None  # time of the last parameter change, None = nothing to do
        self.moving_year_slider = False  # True while the animation itself moves the year slider

        self.add_playback_controls()
        self.add_number_controls(number_controls)
        self.add_choice_controls(choice_controls, map_modes)
        self.rebuild()

        self.last_tick = clock.perf_counter()
        self.timer = fig.canvas.new_timer(interval=int(self.TICK_SECONDS * 1000))
        self.timer.add_callback(self.tick)
        self.timer.start()
        fig.canvas.mpl_connect("close_event", lambda event: self.timer.stop())

    # ---- building the controls ------------------------------------------
    def add_playback_controls(self):
        self.play_button = Button(self.fig.add_axes([0.03, 0.335, 0.08, 0.04]), "Play")
        self.play_button.on_clicked(lambda event: self.toggle_play())
        self.restart_button = Button(self.fig.add_axes([0.12, 0.335, 0.08, 0.04]), "Restart")
        self.restart_button.on_clicked(lambda event: self.restart())

        self.year_slider = Slider(self.fig.add_axes([0.32, 0.345, 0.50, 0.025]), "year", 0, 1, valinit=0)
        self.year_slider.on_changed(self.on_year_dragged)

        # the speed slider moves in powers of ten: from 0.1 to 10 years per second
        self.speed_slider = Slider(self.fig.add_axes([0.32, 0.30, 0.50, 0.025]), "speed [years per second]",
                                   -1, 1, valinit=np.log10(self.speed))
        self.speed_slider.valtext.set_text(format_number(self.speed))
        self.speed_slider.on_changed(self.on_speed_changed)

    def add_number_controls(self, controls):
        self.number_sliders = []
        for row, control in enumerate(controls):
            ax = self.fig.add_axes([0.20, 0.225 - 0.05 * row, 0.27, 0.03])
            value = getattr(self.settings, control.field)
            if control.log_scale:
                slider = Slider(ax, control.label, np.log10(control.minimum), np.log10(control.maximum),
                                valinit=np.log10(value))
                slider.valtext.set_text(format_number(value))
            else:
                slider = Slider(ax, control.label, control.minimum, control.maximum,
                                valinit=value, valstep=control.step)
            self.connect_number_slider(slider, control)
            self.number_sliders.append(slider)

    def connect_number_slider(self, slider, control):
        def changed(slider_value):
            value = 10**slider_value if control.log_scale else slider_value
            if control.log_scale:
                slider.valtext.set_text(format_number(value))
            self.change_setting(control.field, value)

        slider.on_changed(changed)

    def add_choice_controls(self, controls, map_modes):
        groups = [(c.label, c.options, getattr(self.settings, c.field), self.choice_callback(c.field))
                  for c in controls]
        if map_modes:
            groups.append(("map shows", map_modes, self.plot.map_mode, self.on_map_mode_changed))

        self.radio_groups = []
        for column, (label, options, current, callback) in enumerate(groups):
            ax = self.fig.add_axes([0.56 + 0.11 * column, 0.02, 0.105, 0.22])
            ax.set_title(label, fontsize=9)
            radio = RadioButtons(ax, options, active=list(options).index(current))
            for text in radio.labels:
                text.set_fontsize(8)
            radio.on_clicked(callback)
            self.radio_groups.append(radio)

    def choice_callback(self, field):
        def chosen(option):
            self.change_setting(field, option)

        return chosen

    # ---- reacting to the controls ---------------------------------------
    def change_setting(self, field, value):
        """Store the new value; the simulation restarts a moment later (see tick)."""
        self.settings = replace(self.settings, **{field: value})
        self.rebuild_requested_at = clock.perf_counter()

    def on_year_dragged(self, year):
        if self.moving_year_slider:  # the animation moved the slider, not the mouse
            return
        self.set_playing(False)
        self.display_year = year
        self.show_current_snapshot()

    def on_speed_changed(self, slider_value):
        self.speed = 10**slider_value
        self.speed_slider.valtext.set_text(format_number(self.speed))

    def on_map_mode_changed(self, mode):
        self.plot.set_map_mode(mode)
        self.show_current_snapshot()

    def toggle_play(self):
        if not self.playing and self.display_year >= self.cache.last_year:
            self.display_year = 0.0  # at the end, Play replays from the start
        self.set_playing(not self.playing)

    def set_playing(self, playing):
        self.playing = playing
        self.play_button.label.set_text("Pause" if playing else "Play")
        self.last_tick = clock.perf_counter()

    def restart(self):
        self.display_year = 0.0
        self.show_current_snapshot()

    # ---- the heart of the viewer ------------------------------------------
    def rebuild(self):
        """Start a new simulation with the current settings."""
        self.rebuild_requested_at = None
        model, solver, end_time = self.build_simulation(self.settings)
        self.cache = SnapshotCache(model, solver, end_time, self.SNAPSHOT_YEARS * SECONDS_PER_YEAR)
        self.plot.show_model(model, end_time)
        self.year_slider.valmax = self.cache.last_year
        self.year_slider.ax.set_xlim(0, self.cache.last_year)
        self.display_year = 0.0
        self.show_current_snapshot()

    def tick(self):
        """Called by the timer 25 times per second."""
        now = clock.perf_counter()
        elapsed = min(now - self.last_tick, 2 * self.TICK_SECONDS)  # if drawing is slow, play slower
        self.last_tick = now

        if self.rebuild_requested_at is not None and now - self.rebuild_requested_at > self.REBUILD_DELAY:
            self.rebuild()

        if self.playing:
            self.display_year = min(self.display_year + self.speed * elapsed, self.cache.last_year)
            self.show_current_snapshot()
            if self.display_year >= self.cache.last_year:
                self.set_playing(False)

    def show_current_snapshot(self):
        index = self.cache.index_at(self.display_year)
        time, pressure = self.cache.get(index)
        self.plot.show_snapshot(time, pressure, self.cache.years[: index + 1], self.cache.highest[: index + 1])

        self.moving_year_slider = True
        self.year_slider.set_val(time / SECONDS_PER_YEAR)
        self.moving_year_slider = False
        self.fig.canvas.draw_idle()


# --------------------------------------------------------------------------
# Windows for 1D and 2D: plots in the upper part, controls in the lower part
# --------------------------------------------------------------------------
FIGURE_SIZE = (12, 7.6)  # inches; small enough for a laptop screen


def open_viewer_1d(settings, build_simulation, number_controls, choice_controls, show=True):
    fig = plt.figure(figsize=FIGURE_SIZE)
    plot = Plot1D(fig,
                  ax_curve=fig.add_axes([0.07, 0.62, 0.80, 0.29]),
                  ax_strip=fig.add_axes([0.07, 0.47, 0.80, 0.10]),
                  colorbar_ax=fig.add_axes([0.89, 0.47, 0.012, 0.44]))
    viewer = InteractiveViewer(fig, plot, settings, build_simulation, number_controls, choice_controls)
    if show:
        plt.show()
    return viewer


def open_viewer_2d(settings, build_simulation, number_controls, choice_controls, show=True):
    fig = plt.figure(figsize=FIGURE_SIZE)
    plot = Plot2D(fig,
                  ax_map=fig.add_axes([0.04, 0.45, 0.33, 0.46]),
                  ax_history=fig.add_axes([0.53, 0.47, 0.44, 0.42]),
                  colorbar_ax=fig.add_axes([0.38, 0.47, 0.012, 0.42]))
    viewer = InteractiveViewer(fig, plot, settings, build_simulation, number_controls, choice_controls,
                               map_modes=MAP_MODES)
    if show:
        plt.show()
    return viewer
