"""
Save a simulation as an animation file (.gif or .mp4), without opening a window.

To watch, pause, replay and change parameters while it runs, use the
interactive viewer instead (interactive.py; run_1d.py / run_2d.py without --save).
"""

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from .parameters import MPA, SECONDS_PER_YEAR
from .plots import MAP_RISE, Plot1D, Plot2D


def save_animation_1d(model, make_snapshots, n_frames, path, end_time, fps=15):
    """make_snapshots() must return a generator of (time, pressure), see solvers.simulate."""
    fig, (ax_curve, ax_strip) = plt.subplots(2, 1, figsize=(10, 6.5), gridspec_kw={"height_ratios": [3, 1]})
    plot = Plot1D(fig, ax_curve, ax_strip)
    plot.show_model(model, end_time)

    def draw_frame(snapshot):
        time, pressure = snapshot
        plot.show_snapshot(time, pressure)

    write_animation(fig, draw_frame, make_snapshots, n_frames, path, fps)


def save_animation_2d(model, make_snapshots, n_frames, path, end_time, map_mode=MAP_RISE, fps=15):
    fig, (ax_map, ax_history) = plt.subplots(1, 2, figsize=(12, 5))
    plot = Plot2D(fig, ax_map, ax_history)
    plot.show_model(model, end_time)
    plot.set_map_mode(map_mode)
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    years, highest = [], []  # highest pressure so far, for the line on the right

    def draw_frame(snapshot):
        time, pressure = snapshot
        years.append(time / SECONDS_PER_YEAR)
        highest.append(pressure.max() / MPA)
        plot.show_snapshot(time, pressure, years, highest)

    write_animation(fig, draw_frame, make_snapshots, n_frames, path, fps)


def write_animation(fig, draw_frame, make_snapshots, n_frames, path, fps):
    animation = FuncAnimation(fig, draw_frame, frames=make_snapshots, save_count=n_frames,
                              repeat=False, cache_frame_data=False)
    writer = "ffmpeg" if str(path).endswith(".mp4") else "pillow"
    animation.save(path, writer=writer, fps=fps)
    plt.close(fig)
    print(f"animation saved to {path}")
