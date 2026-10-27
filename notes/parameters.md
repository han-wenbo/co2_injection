# Model parameters: what the equation needs, the values we use, and where they come from

The pressure equation of the project,

```
φ c_t ∂p/∂t = ∇·( k/μ ∇p ) + q
```

contains only four physical coefficients (φ, c_t, k, μ) and one source term (q).
Everything else in the simulator is geometry, initial/boundary conditions, or the
safety criterion the pressure is compared with. Section 1 lists exactly what the
model needs, section 2 the values in `co2sim/parameters.py` with their literature
ranges, section 3 what is *not* needed, section 4 the corrections to the draft
literature table, and section 5 the verified references.

## 1. What the model needs

| group | parameter | symbol | used for |
|---|---|---|---|
| rock | porosity | φ | storage term φ c_t |
| rock | permeability (field k(x) or k(x, y)) | k | flow term k/μ |
| rock | permeability of a sealing fault | k_fault | scenarios with faults |
| rock + fluid | total compressibility = pore + brine | c_t = c_r + c_f | storage term φ c_t |
| fluid | viscosity of the resident brine | μ | flow term k/μ |
| fluid | CO₂ density at reservoir conditions | ρ_CO₂ | converts the injected mass rate [Mt/yr] into a volume rate [m³/s]; q = volume rate / cell volume |
| well | injection rate per well | ṁ | source term q |
| well | injection duration, well positions | — | source term q |
| geometry | aquifer thickness (and strip width in 1D) | H (W) | cell volume, face areas |
| geometry | domain size | L | must be large enough for the pressure footprint (see 2.3) |
| initial state | depth and brine density | z, ρ_brine | natural (hydrostatic) pressure p0 = ρ_brine g z |
| safety | fracture-pressure gradient and safety factor | G_frac, s | p_max = s · G_frac · z; the run is "safe" while p < p_max everywhere |
| boundary | open or closed edges | — | scenario choice, not a measured quantity |

Units of q: dividing the equation by φ c_t shows that q has units of 1/s, i.e.
m³ of fluid per m³ of rock per second. The literature reports the *well* rate
(Mt/yr or kg/s), and the code computes q = (ṁ / ρ_CO₂) / V_cell itself. A
value of q "per unit volume" is meaningless without the cell size.

## 2. Values used in the simulator and their literature support

Literature values were read from the papers themselves where we have the full
text (Burke 2011 [B11], Zhou et al. 2008 [Z08], Zhou & Birkholzer 2011 [ZB11]),
otherwise from the verified citation of the paper (marked "draft table").

### 2.1 Coefficients of the equation

| parameter | our value (`parameters.py`) | literature values and ranges | source |
|---|---|---|---|
| porosity φ | 0.20 | 0.12 (base case), 0.10 (analytical estimate); 0.05–0.35 for candidate sandstones | [Z08] Table 1; [ZB11]; general range per [BK05] (draft table) |
| permeability k | 100 mD = 9.87·10⁻¹⁴ m² | 1·10⁻¹³ m² base case, sensitivity 5·10⁻¹⁴ – 1·10⁻¹² m²; USGS "Class II" favourable targets 1 mD – 1 D | [Z08] Table 2; [B11] permeability classes |
| sealing fault k_fault | 1·10⁻¹⁸ m² (≈ 10⁻³ mD) | measured sealing unit 1.0·10⁻¹⁸ m²; sensitivity 10⁻²⁰ – 10⁻¹⁷ m² | [Z08] |
| brine viscosity μ | 5.0·10⁻⁴ Pa·s (0.5 cP) | 0.5·10⁻³ Pa·s; 0.5–0.8 cP depending on temperature and salinity | [ZB11]; [B15] (draft table) |
| total compressibility c_t | 1.0·10⁻⁹ Pa⁻¹ | pore 4.5·10⁻¹⁰ + brine 3.5·10⁻¹⁰ = 8.0·10⁻¹⁰ Pa⁻¹; pore 3.7·10⁻¹⁰ + brine 3.4·10⁻¹⁰ = 7.1·10⁻¹⁰ Pa⁻¹; pore compressibility varied 4.5·10⁻¹¹ – 4.5·10⁻⁹ Pa⁻¹ | [Z08] Table 1 and 2; [ZB11] |

Our c_t is at the upper end of the Birkholzer-group values (they would give
D ≈ 1.2–1.4 m²/s instead of our 1.0 m²/s); the difference is well inside the
spread of the sensitivity studies.

### 2.2 Source term

| parameter | our value | literature values | source |
|---|---|---|---|
| injection rate per well ṁ | 1.0 Mt/yr (adjustable 0.1–5 in the viewer) | Sleipner ≈ 1 Mt/yr; 1 and 2 Mt/yr per well; 3.8 Mt/yr (120 kg/s) single site; 5 Mt/yr per site × 20 sites (Illinois Basin scenario) | [BK05]; [U22]; [Z08]; [ZB11] |
| CO₂ density ρ_CO₂ | 700 kg/m³ (at ≈ 15 MPa, ≈ 50 °C) | 628 kg/m³ at 25.5 MPa and 93 °C (2 440 m); 600–750 kg/m³ at 1–2 km depth | [B11] Table 1 |
| volume rate per well | 1 Mt/yr → 0.045 m³/s ≈ 3 900 m³/day | — | computed |
| injection duration | 20 yr (adjustable 1–35) | 30 yr; 50 yr; "time scale of interest, usually 30–50 years" | [Z08]; [ZB11] |
| simulated time | 40 yr (pressure recovery after shut-in) | up to 1 000 yr post-injection in [U22] | — |

### 2.3 Geometry and boundary

| parameter | our value | literature values | source |
|---|---|---|---|
| aquifer thickness H | 100 m | 250 m (base case); Mt. Simon Sandstone 300–700 m | [Z08]; [ZB11] |
| 1D strip width W | 5 km | — (modelling choice) | — |
| depth z | 1 500 m | 1 200 m; USGS storage window 914–3 962 m (3 000–13 000 ft), base case 2 438 m (8 000 ft); CO₂ stays dense only below ≈ 800–900 m | [Z08]; [B11] |
| temperature (not in the equation) | ≈ 50 °C implied by ρ_CO₂ and μ | 45 °C at 1 200 m; 30 °C/km gradient → 93 °C at 2 438 m | [Z08]; [B11] |
| domain size L | 1D: 50 km; 2D: 50 × 50 km (validation: 100–200 km) | radial extents 10, 20 (base), 30, 100 km; pressure footprint R ≈ √(8 D t) ≈ 190 km after 50 yr for D = 2.8 m²/s | [Z08] Table 2; [ZB11] |
| boundary | open (p = p0 at the edge) or closed (no flow) | closed compartment vs. open basin is the central distinction of [Z08]; the 10 km closed case exceeds the limit after 8 yr | [Z08] |

Domain-size rule: with our D ≈ 1 m²/s the footprint √(8 D t) is ≈ 16 km after
1 yr and ≈ 70 km after 20 yr. The 50 km default domain therefore behaves as an
aquifer in contact with a constant-pressure region 25 km from the well, not as
an infinite aquifer. Measured effect on the single-well peak pressure at 20 yr:
17.75 MPa (50 km domain) vs. 17.87 MPa (100 km and 150 km), i.e. 0.12 MPa or
about 5 % of the pressure rise. `validate.py` uses 100–200 km domains for that reason.

### 2.4 Initial pressure and safety limit

| parameter | our value | literature values | source |
|---|---|---|---|
| brine density ρ_brine | 1 050 kg/m³ → hydrostatic gradient 10.3 kPa/m | 0.465 psi/ft = 10.5 kPa/m; ≈ 10.0 kPa/m (12.0 MPa at 1 200 m) | [B11]; [Z08] |
| natural pressure p0 | 15.45 MPa at 1 500 m | 12.0 MPa at 1 200 m; 25.5 MPa at 2 438 m; 8–25 MPa over typical depths | [Z08]; [B11]; draft table |
| fracture-pressure gradient G_frac | 17 kPa/m | 12.9 kPa/m (Mt. Simon, 130 % of hydrostatic) and 18.1 kPa/m (Dundee Limestone, 181 %); 18.8–22.3 kPa/m (16–19 ppg, Eaton method, at 2 438 m), 20.9 kPa/m used | [Z08]; [B11] Fig. 10–11 |
| safety factor s | 0.9 | injection pressure ≤ 90 % of fracture pressure (US EPA Class VI rule, 40 CFR 146.88(a)); [Z08] instead keep Δp ≤ 50 % of p0 | EPA; [Z08] |
| maximum pressure p_max | 0.9 × 17 kPa/m × 1 500 m = 22.95 MPa = 1.49 p0 | sustainable pressures 15.5 and 21.7 MPa at 1 200 m (1.3 and 1.8 × p0) | [Z08] |
| allowed buildup Δp_max = p_max − p0 | 7.5 MPa (0.49 p0) | 3.5 and 9.7 MPa; chosen value 6.0 MPa (50 % of p0 at 1 200 m); 25.5 MPa at 2 438 m with no safety factor | [Z08]; [B11] |

### 2.5 Derived quantity

| quantity | formula | our value |
|---|---|---|
| pressure diffusivity D | k / (φ μ c_t) | 0.99 m²/s (Birkholzer-group values give 2.8 m²/s with φ = 0.10, c_t = 7.1·10⁻¹⁰) |

## 3. Not needed in this (single-phase) model

| quantity | why it is not needed | when it would be |
|---|---|---|
| CO₂ viscosity (5.0·10⁻⁵ Pa·s at 25.5 MPa, 93 °C [B11]; 0.03–0.08 cP) | the equation tracks the pressure of the resident brine; CO₂ only fills the plume near the well | two-phase model, or a near-well mobility correction |
| CO₂ compressibility (1.66·10⁻⁸ Pa⁻¹ [B11]) | same reason; the aquifer is ≫ 99 % brine | two-phase model |
| relative permeability, capillary pressure (van Genuchten α, m in [Z08]) | two-phase quantities | plume simulation |
| seal (caprock) permeability | top and bottom are assumed perfectly sealed (conservative: no vertical pressure bleed-off) | 3D or layered model; [Z08] shows seal k ≥ 10⁻¹⁸ m² relieves pressure noticeably |
| temperature, salinity | enter only through ρ_CO₂, ρ_brine and μ, which we set directly | property correlations |

## 4. Corrections to the draft literature table

1. **Citation of Birkholzer, Zhou & Tsang (2009).** Correct: *Int. J. Greenhouse Gas Control* 3(2), 181–194, doi:10.1016/j.ijggc.2008.08.002. "3(6), 745–756" is a different paper: Birkholzer & Zhou (2009), *Basin-scale hydrogeologic impacts of CO₂ storage: Capacity and regulatory implications*, IJGGC 3(6), 745–756.
2. **USGS report has an author.** Burke, L.A., 2011, *Carbon dioxide fluid-flow modeling and injectivity calculations*, USGS Scientific Investigations Report 2011–5083, 16 p.
3. **Units of q.** In this equation q is in 1/s (volume per rock volume per second), not kg/(m³·s); and a per-volume value depends on the grid cell size. The range "10⁻⁴–10⁻² kg/(m³·s)" could not be traced to Bentham & Kirby (2005) or any other listed source — drop it and specify the well rate instead.
4. **Rate conversion.** 1 Mt/yr is ≈ 3 900–4 600 m³/day (ρ_CO₂ = 700–600 kg/m³) and 2 Mt/yr ≈ 7 800–9 100 m³/day; "5 000–20 000 m³/day" corresponds to about 1–5 Mt/yr.
5. **"Max pressure buildup Δp_max = 1.3–1.8 × p0, baseline 1.5 × p0 (15 MPa)" mixes two quantities.** 1.3–1.8 × p0 is the *fracture (closure) pressure* itself ([Z08]: 130 % and 181 % of the hydrostatic gradient). The *buildup* is the difference: for p0 = 10 MPa and p_max = 1.5 p0, Δp_max = 5 MPa, not 15 MPa. [Z08] chose Δp_max = 0.5 p0 (6 MPa at 1 200 m).
6. **Domain length "1–10 km baseline".** Too small for an open aquifer: the pressure footprint √(8 D t) passes 10 km within months. [Z08] used 10–100 km and showed that a 10 km closed compartment exceeds the limit after 8 yr. The USGS report normalises its time-scale calculation to a 1 km distance; it is not a domain recommendation.
7. **Missing rows that the model needs:** aquifer thickness, depth, brine density (→ p0), CO₂ density (→ volume rate), fracture gradient + safety factor, injection duration, fault/seal permeability.
8. **CO₂ viscosity** is correct as quoted ([B11]: 5.0·10⁻⁵ Pa·s) but is not used by a single-phase pressure model (section 3).
9. Values that check out: porosity, permeability, brine viscosity, total compressibility, initial pressure ranges; the citations Ansari et al. (2022), Wang et al. (2024), Urych et al. (2022), Jin et al. (2012), Bentham & Kirby (2005), Szulczewski et al. (2012), Zhou et al. (2008) all exist as cited. The specific ranges attributed to Birkholzer et al. (2009), Birkholzer et al. (2015), Jin et al. (2012) and Szulczewski et al. (2012) were not checked against the full texts.

## 5. References

Full texts read for this document:

- [B11] Burke, L.A. (2011). *Carbon dioxide fluid-flow modeling and injectivity calculations.* U.S. Geological Survey Scientific Investigations Report 2011–5083, 16 p. https://pubs.usgs.gov/sir/2011/5083/ — reference given in the lecture.
- [Z08] Zhou, Q., Birkholzer, J.T., Tsang, C.-F., Rutqvist, J. (2008). *A method for quick assessment of CO₂ storage capacity in closed and semi-closed saline formations.* Int. J. Greenhouse Gas Control 2(4), 626–639. Open-access report LBNL-63820: https://www.osti.gov/servlets/purl/932805
- [ZB11] Zhou, Q., Birkholzer, J.T. (2011). *On scale and magnitude of pressure build-up induced by large-scale geologic storage of CO₂.* Greenhouse Gases: Science and Technology 1(1), 11–20. Preprint: https://digital.library.unt.edu/ark:/67531/metadc837589/

Citations verified (existence, journal, volume, pages):

- [B09] Birkholzer, J.T., Zhou, Q., Tsang, C.-F. (2009). *Large-scale impact of CO₂ storage in deep saline aquifers: A sensitivity study on pressure response in stratified systems.* Int. J. Greenhouse Gas Control 3(2), 181–194. doi:10.1016/j.ijggc.2008.08.002
- [B15] Birkholzer, J.T., Oldenburg, C.M., Zhou, Q. (2015). *CO₂ migration and pressure evolution in deep saline aquifers.* Int. J. Greenhouse Gas Control 40, 203–220.
- [S12] Szulczewski, M.L., MacMinn, C.W., Herzog, H.J., Juanes, R. (2012). *Lifetime of carbon capture and storage as a climate-change mitigation technology.* PNAS 109(14), 5185–5189.
- [BK05] Bentham, M., Kirby, G.A. (2005). *CO₂ storage in saline aquifers.* Oil & Gas Science and Technology 60(3), 559–567.
- [J12] Jin, M., Pickup, G., Mackay, E., Todd, A., Monaghan, A., Naylor, M. (2012). *Static and dynamic estimates of CO₂-storage capacity in two saline formations in the UK.* SPE Journal 17(4), 1108–1118. doi:10.2118/131609-PA
- [U22] Urych, T., Chećko, J., Magdziarczyk, M., Smoliński, A. (2022). *Numerical simulations of carbon dioxide storage in selected geological structures in north-western Poland.* Frontiers in Energy Research 10, 827794. doi:10.3389/fenrg.2022.827794
- [A22] Ansari, E., Holubnyak, E., Hasiuk, F., Raney, J. (2022). *Comparing CO₂ and variable-TDS waste-fluid injection into deep saline geologic formations.* Frontiers in Energy Research 10, 864013. doi:10.3389/fenrg.2022.864013
- [W24] Wang, Q., Zhang, D., Li, Y., Li, C., Tang, H. (2024). *Numerical simulation study of CO₂ storage capacity in deep saline aquifers.* Science and Technology for Energy Transition 79, 12. doi:10.2516/stet/2024005
- US EPA (2010). *Federal Requirements Under the Underground Injection Control (UIC) Program for Carbon Dioxide (CO₂) Geologic Sequestration (GS) Wells; Final Rule.* Federal Register 75(237), 77230 (10 Dec 2010), codified as 40 CFR 146.88(a): injection pressure must not exceed 90 % of the fracture pressure of the injection zone (Class VI wells). https://www.law.cornell.edu/cfr/text/40/146.88
