"""
Physical parameters of the aquifer, the fluid and the safety limit.

Everything is in SI units (metre, second, pascal, kilogram) unless the name
says otherwise, e.g. `rate_mt_per_year`.

The default values describe a *typical* deep saline aquifer used for CO2
storage: a sandstone layer about 1.5 km deep and 100 m thick whose pores are
already full of salty water (brine). They are representative textbook values,
not the data of one specific site.
"""

from dataclasses import dataclass

# --------------------------------------------------------------------------
# Unit conversions
# --------------------------------------------------------------------------
SECONDS_PER_YEAR = 365.25 * 24 * 3600  # [s]
MPA = 1e6  # 1 megapascal in pascal
MILLIDARCY = 9.869233e-16  # 1 millidarcy (oil-field unit of permeability) in m²
GRAVITY = 9.81  # [m/s²]


@dataclass(frozen=True)
class AquiferParameters:
    """Rock, fluid and pressure-limit properties.

    The pressure equation from the project slides is

        φ c_t ∂p/∂t = ∂/∂x ( k/μ ∂p/∂x ) + q

    φ, c_t, k and μ are fields of this class. q comes from the wells
    (see wells.py).
    """

    # ---- rock -------------------------------------------------------------
    # φ [-]: fraction of the rock volume that is pore space (sandstone: 0.1–0.3)
    porosity: float = 0.20

    # k [m²]: how easily fluid flows *through* the rock (good sandstone: 10–1000 mD).
    # This is the default value; the models also accept a k that varies in space.
    permeability: float = 100 * MILLIDARCY

    # H [m]: vertical thickness of the aquifer layer
    thickness: float = 100.0

    # [m]: depth of the aquifer below the surface
    depth: float = 1500.0

    # ---- fluid ------------------------------------------------------------
    # μ [Pa·s]: viscosity of the brine at reservoir temperature (~50 °C)
    viscosity: float = 5e-4

    # c_t [1/Pa]: how much the rock and the brine "give" when squeezed.
    # Rock ~4.5e-10 plus brine ~4e-10, so roughly 1e-9 per pascal.
    total_compressibility: float = 1e-9

    # [kg/m³]: brine density, sets the natural (hydrostatic) pressure at depth
    brine_density: float = 1050.0

    # [kg/m³]: CO2 density at reservoir conditions (~15 MPa, ~50 °C),
    # converts an injected mass (tonnes) into an injected volume (m³)
    co2_density: float = 700.0

    # ---- safety limit -----------------------------------------------------
    # [Pa/m]: the pressure that fractures the rock grows by about 17 kPa
    # per metre of depth (typical range 14–20 kPa/m)
    fracture_gradient: float = 17e3

    # Stay below 90 % of the fracture pressure. This is the rule the US EPA
    # uses for CO2 injection wells (40 CFR 146.88).
    safety_factor: float = 0.9

    # ---- derived quantities ----------------------------------------------
    @property
    def initial_pressure(self) -> float:
        """p0 [Pa]: natural pressure before injection (weight of the brine above)."""
        return self.brine_density * GRAVITY * self.depth

    @property
    def fracture_pressure(self) -> float:
        """[Pa]: pressure at which the rock would crack."""
        return self.fracture_gradient * self.depth

    @property
    def max_safe_pressure(self) -> float:
        """[Pa]: the pressure must stay below this everywhere."""
        return self.safety_factor * self.fracture_pressure

    @property
    def diffusivity(self) -> float:
        """D = k / (φ μ c_t) [m²/s]: how fast pressure spreads.

        After a time t the pressure signal has travelled roughly sqrt(D·t) metres.
        """
        return self.permeability / (self.porosity * self.viscosity * self.total_compressibility)
