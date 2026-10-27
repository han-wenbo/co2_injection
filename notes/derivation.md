# From the pressure equation to the code: a step-by-step derivation

We start from the equation on the project slides and derive, step by step,
exactly what the simulator computes. The symbols of the equation are kept
throughout; every new symbol is defined where it first appears.

$$
\phi\, c_t\, \frac{\partial p(x,t)}{\partial t} = \frac{\partial}{\partial x}\left( \frac{k(x)}{\mu}\, \frac{\partial p(x,t)}{\partial x} \right) + q(x,t)
$$

| symbol | meaning | unit | default value |
|---|---|---|---|
| $p(x,t)$ | pressure at position $x$ and time $t$ | Pa | starts at $p_0 = 15.45$ MPa |
| $\phi$ | porosity: pore fraction of the rock | – | 0.20 |
| $c_t$ | total compressibility of rock and brine | 1/Pa | $10^{-9}$ |
| $k(x)$ | permeability: how easily fluid flows through the rock | m² | $9.87\cdot10^{-14}$ (100 mD) |
| $\mu$ | brine viscosity | Pa·s | $5\cdot10^{-4}$ |
| $q(x,t)$ | injected fluid volume per rock volume per second | 1/s | non-zero only at the well |

---

## Step 0 – What each term means, and its unit

| term | meaning | unit |
|---|---|---|
| $\phi c_t \,\partial p/\partial t$ | **storage**: extra fluid packed into the pores as the pressure rises | (1/Pa)(Pa/s) = 1/s |
| $\partial_x\!\left(\frac{k}{\mu}\partial_x p\right)$ | **flow**: net fluid arriving from the neighbours (Darcy's law) | (m²/(Pa·s))(Pa/m)/m = 1/s |
| $q$ | **source**: fluid pushed in by the well | 1/s |

All three terms are "volume of fluid per volume of rock per second", so the
equation is a volume balance: *stored = flowed in + injected*.

## Step 1 – Subtract the natural pressure

Before injection the pressure is the hydrostatic pressure of the brine column,

$$
p_0 = \rho_w\, g\, z ,
$$

with brine density $\rho_w$, gravity $g$ and depth $z$. The aquifer is
horizontal, so $p_0$ is the same everywhere and does not change in time. Define
the **overpressure** (pressure rise)
$$
\Delta p(x,t) = p(x,t) - p_0 .
$$

Because $\partial p_0/\partial t = 0$ and $\partial p_0/\partial x = 0$, inserting
$p = p_0 + \Delta p$ gives exactly the same equation for $\Delta p$:
$$
\phi\, c_t\, \frac{\partial \Delta p}{\partial t} = \frac{\partial}{\partial x}\left( \frac{k(x)}{\mu}\, \frac{\partial \Delta p}{\partial x} \right) + q(x,t), \qquad \Delta p(x,0) = 0 .
$$

The code solves for $\Delta p$ (variable `overpressure`) and adds $p_0$ only for
plotting.

## Step 2 – Uniform rock: a diffusion equation

If $k$ does not depend on $x$, it can be taken out of the derivative. Dividing by
$\phi c_t$:

$$
\frac{\partial \Delta p}{\partial t} = D\, \frac{\partial^2 \Delta p}{\partial x^2} + \frac{q}{\phi\, c_t},
\qquad
D = \frac{k}{\phi\, \mu\, c_t} .
$$

This is the heat equation with a source. $D$ is the pressure diffusivity:
$D \approx 0.99$ m²/s with our values, and the pressure signal reaches a distance
of about $\sqrt{8 D t}$, i.e. 16 km after one year. The code keeps the general
$k(x)$; this form is used to interpret results and to compare with exact solutions.

## Step 3 – Cut the aquifer into cells

Split the length $L$ into $N$ equal cells:

$$
\Delta x = \frac{L}{N}, \qquad x_i = \left(i + \tfrac12\right)\Delta x, \qquad i = 0, 1, \dots, N-1 .
$$

Cell $i$ lies between the faces $x_{i-1/2}$ and $x_{i+1/2}$. In 1D the aquifer is
a strip of width $W$ and thickness $H$, so every face has the area $A = W H$ and
every cell the volume $V = A\,\Delta x$. The unknown of cell $i$ is its average
overpressure $\Delta p_i(t)$.

(Defaults: $L = 50$ km, $N = 251$, $\Delta x = 199$ m, $W = 5$ km, $H = 100$ m.)

## Step 4 – Integrate the equation over one cell (finite volumes)

Multiply the equation by $A$ and integrate over cell $i$, from $x_{i-1/2}$ to
$x_{i+1/2}$. The integral of a derivative is the difference of its values at the
two ends:

$$
\underbrace{\phi\, c_t\, A\, \Delta x}_{S_i}\; \frac{d \Delta p_i}{dt}
= A \left[ \frac{k}{\mu} \frac{\partial \Delta p}{\partial x} \right]_{x_{i+1/2}}
- A \left[ \frac{k}{\mu} \frac{\partial \Delta p}{\partial x} \right]_{x_{i-1/2}}
+ \underbrace{\int_{x_{i-1/2}}^{x_{i+1/2}} q\, A\, dx}_{Q_i(t)}
$$

- $S_i = \phi c_t V$ is the **storage** of the cell [m³/Pa]: how much extra fluid
  fits in when its pressure rises by 1 Pa.
- Fluid moves towards lower pressure (Darcy's law). The first bracket is the fluid
  entering through the right face, the second (with its minus sign) the fluid
  entering through the left face, in m³/s.
- $Q_i(t)$ is the volume injected into cell $i$ per second.

So far nothing is approximated: this is the exact volume balance of one cell.

## Step 5 – Flow through a face: transmissibility and harmonic mean

We approximate the face flow from the two neighbouring cell values. Between the
centres of cells $i$ and $i+1$ the fluid crosses two half cells in a row, with
permeabilities $k_i$ and $k_{i+1}$. Call the pressure at the face $\Delta p_f$. The
same flow $F$ passes through both halves:

$$
F = \frac{k_i A}{\mu}\, \frac{\Delta p_f - \Delta p_i}{\Delta x/2} = \frac{k_{i+1} A}{\mu}\, \frac{\Delta p_{i+1} - \Delta p_f}{\Delta x/2} .
$$

Eliminating $\Delta p_f$ gives

$$
F = T_{i+1/2}\, \left( \Delta p_{i+1} - \Delta p_i \right),
\qquad
T_{i+1/2} = \frac{k_{i+1/2}\, A}{\mu\, \Delta x},
\qquad
k_{i+1/2} = \frac{2\, k_i\, k_{i+1}}{k_i + k_{i+1}} .
$$

$T_{i+1/2}$ is the **transmissibility** of the face [m³/(Pa·s)] and $k_{i+1/2}$
the **harmonic mean**: if one side is almost impermeable (a fault), the face is
almost impermeable too. (Code: `finite_volume.harmonic_mean`, `face_transmissibility`.)

## Step 6 – The well as a source

A well at position $x_w$ lies in one cell $i_w$. We spread its injection over that
cell:

$$
q(x,t) =
\begin{cases}
Q(t)/V & \text{in cell } i_w \\
0 & \text{elsewhere}
\end{cases}
\qquad\Longrightarrow\qquad
Q_i(t) = \int q\, A\, dx =
\begin{cases}
Q(t) & i = i_w \\
0 & i \ne i_w
\end{cases}
$$

The well rate is given as a CO₂ mass rate $\dot m$; dividing by the CO₂ density
$\rho_{\mathrm{CO_2}}$ gives the volume rate:

$$
Q(t) = \frac{\dot m}{\rho_{\mathrm{CO_2}}} \quad \text{while the pump runs, otherwise } 0 .
$$

Example: 1 Mt/yr gives $\dot m = 31.7$ kg/s, $Q = 0.045$ m³/s and, in a cell of
$V = 9.96\cdot10^{7}$ m³, $q = 4.5\cdot10^{-10}$ 1/s. (Code: `Well.volume_rate`, `model.source`.)

## Step 7 – One equation per cell

Inserting steps 5 and 6 into step 4:

$$
S_i\, \frac{d \Delta p_i}{dt} = T_{i+1/2}\, \left( \Delta p_{i+1} - \Delta p_i \right) - T_{i-1/2}\, \left( \Delta p_i - \Delta p_{i-1} \right) + Q_i(t) .
$$

This is a system of $N$ ordinary differential equations: *storage × rate of
change = inflow from the right + inflow from the left + well*. (Code:
`model.net_inflow`.)

**Check against the original equation.** For uniform $k$, all $T$ are equal to
$k A/(\mu \Delta x)$. Dividing by $S_i = \phi c_t A \Delta x$:

$$
\frac{d \Delta p_i}{dt} = D\, \frac{\Delta p_{i+1} - 2\,\Delta p_i + \Delta p_{i-1}}{\Delta x^2} + \frac{q_i}{\phi\, c_t},
$$

the standard central-difference approximation of step 2. Its error shrinks like
$\Delta x^2$; `validate.py` measures an order of 2.00.

## Step 8 – Boundary conditions

The first and last cells have a face to the outside.

- **Closed** (sealed edge, no flow, $\partial p/\partial x = 0$): that face carries no flow, $T = 0$.
- **Open** (the aquifer continues and keeps the natural pressure, $p = p_0$, i.e.
  $\Delta p = 0$ outside): the outside lies half a cell away from the cell centre, so
  the boundary face adds the inflow
  $$
  T_b\, \left(0 - \Delta p_0\right), \qquad T_b = \frac{k_0\, A}{\mu\, \Delta x/2}
  $$
  for cell 0, and the same with $k_{N-1}$ and $\Delta p_{N-1}$ for the last cell.

## Step 9 – Matrix form

Collect all cells in a vector $\Delta\mathbf{p} = (\Delta p_0, \dots, \Delta p_{N-1})$:

$$
\mathbf{S}\, \frac{d \Delta\mathbf{p}}{dt} = -\mathbf{L}\, \Delta\mathbf{p} + \mathbf{Q}(t),
$$

with $\mathbf{S} = \mathrm{diag}(S_i)$ and the tridiagonal **flow matrix**

$$
L_{ii} = T_{i-1/2} + T_{i+1/2}\ (+\,T_b \text{ at open ends}), \qquad
L_{i,i+1} = -T_{i+1/2}, \qquad
L_{i,i-1} = -T_{i-1/2} .
$$

$(\mathbf{L}\,\Delta\mathbf{p})_i$ is the net flow *out of* cell $i$. (Code: `model.flow_matrix`.)

## Step 10 – Time stepping

Use time steps $t_n = n\, \Delta t$ and write $\Delta\mathbf{p}^n \approx \Delta\mathbf{p}(t_n)$.
Approximate the time derivative by a difference quotient and evaluate the right-hand
side as a weighted average of the old and the new time, with weight $\theta$:

$$
\mathbf{S}\, \frac{\Delta\mathbf{p}^{n+1} - \Delta\mathbf{p}^{n}}{\Delta t}
= -\mathbf{L} \left[ \theta\, \Delta\mathbf{p}^{n+1} + (1-\theta)\, \Delta\mathbf{p}^{n} \right] + \theta\, \mathbf{Q}^{n+1} + (1-\theta)\, \mathbf{Q}^{n} .
$$

Moving the unknown $\Delta\mathbf{p}^{n+1}$ to the left:

$$
\left( \frac{\mathbf{S}}{\Delta t} + \theta\, \mathbf{L} \right) \Delta\mathbf{p}^{n+1}
= \left( \frac{\mathbf{S}}{\Delta t} - (1-\theta)\, \mathbf{L} \right) \Delta\mathbf{p}^{n} + \theta\, \mathbf{Q}^{n+1} + (1-\theta)\, \mathbf{Q}^{n} .
$$

| $\theta$ | method | per step | error in time |
|---|---|---|---|
| 0 | explicit Euler | one matrix–vector product | $\propto \Delta t$ |
| 1 | backward (implicit) Euler | solve a linear system | $\propto \Delta t$ |
| ½ | Crank–Nicolson | solve a linear system | $\propto \Delta t^2$ |

The matrix on the left never changes, so the code factorises it once (sparse LU)
and reuses it every step. (Code: `solvers.py`.)

## Step 11 – Explicit Euler: why the time step must be small

With $\theta = 0$, cell by cell:

$$
\Delta p_i^{n+1} = \Delta p_i^{n} + \frac{\Delta t}{S_i} \left[ T_{i+1/2} \left( \Delta p_{i+1}^{n} - \Delta p_i^{n} \right) - T_{i-1/2} \left( \Delta p_i^{n} - \Delta p_{i-1}^{n} \right) + Q_i^{n} \right] .
$$

Sorting by old values:

$$
\Delta p_i^{n+1} = \left( 1 - \frac{\Delta t\, (T_{i-1/2} + T_{i+1/2})}{S_i} \right) \Delta p_i^{n}
+ \frac{\Delta t\, T_{i+1/2}}{S_i}\, \Delta p_{i+1}^{n}
+ \frac{\Delta t\, T_{i-1/2}}{S_i}\, \Delta p_{i-1}^{n}
+ \frac{\Delta t}{S_i}\, Q_i^{n} .
$$

The weights of the three old values add up to 1. If all of them are
non-negative, the new value is an average of old values and cannot overshoot; this
holds when

$$
\Delta t \le \frac{S_i}{T_{i-1/2} + T_{i+1/2}} \quad \text{for every cell } i .
$$

For uniform rock this becomes

$$
\Delta t \le \frac{\phi\, c_t\, \mu\, \Delta x^2}{2\, k} = \frac{\Delta x^2}{2 D} \quad (1\text{D}),
\qquad
\Delta t \le \frac{\Delta x^2}{4 D} \quad (2\text{D}, \ \Delta x = \Delta y) .
$$

With our grid: 0.23 days in the interior and 0.155 days in the end cells of an open
aquifer (where $T_b = 2T$). The implicit methods have no such limit and use 18-day
steps. (Code: `max_stable_time_step`.)

## Step 12 – Check: no fluid is created or lost

With closed boundaries, sum the cell equations of step 7 over all cells. Every
interior face appears twice with opposite signs (what leaves one cell enters its
neighbour), so the sum telescopes:

$$
\sum_i S_i\, \frac{d \Delta p_i}{dt} = \sum_i Q_i(t) .
$$

The extra fluid stored in the aquifer equals the fluid injected by the wells. The
discrete schemes keep this exactly (up to rounding errors); `tests/test_solvers.py` checks it for all three
methods.

## Step 13 – Linearity: superposition and the maximum safe rate

The equation is linear in $\Delta p$ and $q$, and $\Delta p(x,0) = 0$. Therefore:

- **Scaling:** multiplying the injection rate by $\alpha$ multiplies $\Delta p$ everywhere by $\alpha$.
- **Superposition:** the overpressure of several wells is the sum of the
  overpressures of each well alone (checked in the tests).

  Consequence: one simulation with rate $Q$ is enough to find the largest rate that
  stays below the allowed rise $\Delta p_{\max}$ (step 14):

$$
Q_{\max} = Q\, \frac{\Delta p_{\max}}{\max_{i,\,n} \Delta p_i^{\,n}} .
$$

(This formula is not in the code yet.)

## Step 14 – Safety check

Add the natural pressure back:

$$
p_i^{\,n} = p_0 + \Delta p_i^{\,n} .
$$

The pressure must stay below 90 % of the fracture pressure (fracture gradient
$G_{\mathrm{frac}}$, safety factor $s$) in every cell at every time step. The 90 % is
the US EPA rule for CO₂ injection wells, 40 CFR 146.88(a) (Federal Register 75, 77230, 2010):

$$
p_{\max} = s\, G_{\mathrm{frac}}\, z = 22.95\ \mathrm{MPa},
\qquad
\max_{i,\,n}\ p_i^{\,n} \le p_{\max}
\iff
\max_{i,\,n}\ \Delta p_i^{\,n} \le \Delta p_{\max} = p_{\max} - p_0 = 7.5\ \mathrm{MPa} .
$$

If this fails, the viewer turns the title red and draws the cells above the limit in red.

## Step 15 – The same in 2D (map view)

In 2D the flow term gets a second direction:

$$
\phi\, c_t\, \frac{\partial p}{\partial t} = \frac{\partial}{\partial x}\left( \frac{k}{\mu} \frac{\partial p}{\partial x} \right) + \frac{\partial}{\partial y}\left( \frac{k}{\mu} \frac{\partial p}{\partial y} \right) + q .
$$

Steps 1–14 repeat with cells $(i, j)$ of size $\Delta x \times \Delta y \times H$:

$$
S_{ij} = \phi\, c_t\, \Delta x\, \Delta y\, H,
\qquad
T^{x}_{i+1/2,\,j} = \frac{k_{i+1/2,\,j}\, \Delta y\, H}{\mu\, \Delta x},
\qquad
T^{y}_{i,\,j+1/2} = \frac{k_{i,\,j+1/2}\, \Delta x\, H}{\mu\, \Delta y},
$$

$$
S_{ij}\, \frac{d \Delta p_{ij}}{dt} =
T^{x}_{i+1/2,j} (\Delta p_{i+1,j} - \Delta p_{ij}) - T^{x}_{i-1/2,j} (\Delta p_{ij} - \Delta p_{i-1,j})
+ T^{y}_{i,j+1/2} (\Delta p_{i,j+1} - \Delta p_{ij}) - T^{y}_{i,j-1/2} (\Delta p_{ij} - \Delta p_{i,j-1})
+ Q_{ij}(t) .
$$

Each cell now exchanges fluid with four neighbours, and $\mathbf{L}$ has five
non-zero diagonals instead of three. (Code: `model_2d.py`.)

---

## Summary: where each step lives in the code

| step | result | code |
|---|---|---|
| 1 | solve for $\Delta p = p - p_0$ | `overpressure` in `solvers.simulate` |
| 3 | grid $\Delta x$, $x_i$ | `finite_volume.cell_centres` |
| 4 | storage $S_i = \phi c_t V$ | `model.storage` |
| 5 | $T_{i+1/2}$, harmonic mean | `finite_volume.face_transmissibility`, `harmonic_mean` |
| 6 | well volume rate $Q = \dot m / \rho_{\mathrm{CO_2}}$ | `Well.volume_rate`, `model.source` |
| 7 | net inflow of each cell | `model.net_inflow` |
| 8 | open / closed boundaries | `boundary_transmissibility` |
| 9 | flow matrix $\mathbf{L}$ | `model.flow_matrix`, `finite_volume.assemble_flow_matrix` |
| 10 | θ-method | `ExplicitEuler`, `BackwardEuler`, `CrankNicolson` |
| 11 | stability limit | `max_stable_time_step` |
| 14 | safety check | `AquiferParameters.max_safe_pressure`, `plots.safety_status` |
