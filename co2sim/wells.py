"""
Wells: the source term q of the pressure equation.

A well pushes fluid into the aquifer (injection, positive rate) or pulls it
out (extraction, negative rate). Everywhere else q = 0.
"""

import math
from dataclasses import dataclass

from .parameters import SECONDS_PER_YEAR


@dataclass(frozen=True)
class Well:
    # position [m]; the 1D model ignores y
    x: float
    y: float = 0.0

    # CO2 mass injected per year, in million tonnes (Mt).
    # For scale: the Sleipner project in Norway injects about 1 Mt per year.
    rate_mt_per_year: float = 1.0

    # the pump is switched on at start_year and off at stop_year
    start_year: float = 0.0
    stop_year: float = math.inf

    def is_active(self, time: float) -> bool:
        """Is the pump running at `time` [s]?"""
        return self.start_year <= time / SECONDS_PER_YEAR < self.stop_year

    def volume_rate(self, time: float, co2_density: float) -> float:
        """Volume pushed into the aquifer per second [m³/s] at `time` [s]."""
        if not self.is_active(time):
            return 0.0
        kilograms_per_second = self.rate_mt_per_year * 1e9 / SECONDS_PER_YEAR
        return kilograms_per_second / co2_density
