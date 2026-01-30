#!/usr/bin/env python3
"""Jet-in-crossflow case with a Lagrangian vapor bubble and three-fluid surface tension.

The setup follows the scaling used in the reference snippet provided by the user: the
CFL number is fixed from the acoustic speed in the air free-stream, the grid spacing is
based on the physical jet diameter, and outputs are scheduled by time step counts.
"""

import json
import math

# -----------------------------------------------------------------------------
# Baseline thermodynamic and flow parameters (physical units)
# -----------------------------------------------------------------------------
pA = 140000.0
rho_air = 1.76
gam_air = 1.4
c1 = math.sqrt(gam_air * pA / rho_air)

pL = pA
velJ = 37.9  # jet velocity (y-direction)
velA = 0.3 * c1  # crossflow velocity (x-direction)
rho_hex = 660.0
rho_vapor = 2.5
mu_air = 1.85e-5
mu_hex = 0.30e-3
mu_vapor = 1.5e-5

# Dense vapor void fraction injected with the jet
vapor_frac = 0.4

# Reference Eulerian bubble parameters
R0ref = 10.0e-6
pv = 2300.0
Ca = (pL - pv) / (rho_hex * velJ**2)
sigma_hex_air = 0.018
We = rho_hex * velJ**2 * R0ref / sigma_hex_air
Re_inv = mu_hex / (rho_hex * velJ * R0ref)

# -----------------------------------------------------------------------------
# Grid / time setup
# -----------------------------------------------------------------------------
djet = 50e-6
x_domain_beg = -10 * djet
x_domain_end = 50 * djet
y_domain_beg = 0.0
y_domain_end = 26.25 * djet
x_length = x_domain_end - x_domain_beg
y_length = y_domain_end - y_domain_beg
Ny = int(1778 / 2)
Nx = int(7112 / 2)
dx = x_length / Nx

time_end = 3.0e-6
cfl = 0.1

dt = cfl * dx / c1
Nt = int(math.ceil(time_end / dt))

# ensure at least one time step
if Nt < 1:
    Nt = 1

# frequency for data output in terms of number of time steps (target ~50 outputs)
save_interval = max(Nt // 50, 1)

# Regularisation
eps = 1.0e-6

jet_patch_length_x = 40 * dx
jet_patch_xc = x_domain_beg + 0.5 * jet_patch_length_x
jet_patch_yc = y_domain_beg + 0.5 * djet

print(
    json.dumps(
        {
            # Logistics -------------------------------------------------------
            "run_time_info": "T",
            "x_domain%beg": x_domain_beg,
            "x_domain%end": x_domain_end,
            "y_domain%beg": y_domain_beg,
            "y_domain%end": y_domain_end,
            "m": int(Nx),
            "n": int(Ny),
            "p": 0,
            "dt": dt,
            "t_step_start": 0,
            "t_step_stop": int(Nt),
            "t_step_save": int(save_interval),
            # Simulation algorithm parameters --------------------------------
            "num_patches": 2,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 3,
            "mpp_lim": "F",
            "mixture_err": "T",
            "bubbles_euler": "T",
            "bubble_model": 2,
            "polytropic": "T",
            "polydisperse": "F",
            "R0_type": 1,
            "thermal": 3,
            "R0ref": R0ref,
            "nb": 1,
            "Ca": Ca,
            "Web": We,
            "Re_inv": Re_inv,
            "time_stepper": 2,
            "weno_order": 3,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "F",
            "mapped_weno": "T",
            "weno_avg": "F",
            "null_weights": "F",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "surface_tension": "T",
            "viscous": "T",
            "elliptic_smoothing": "T",
            "elliptic_smoothing_iters": 50,
            # Boundary conditions ---------------------------------------------
            "bc_x%beg": -16,
            "bc_x%end": -3,
            "bc_y%beg": -16,
            "bc_y%end": -3,
            "num_bc_patches": 1,
            "patch_bc(1)%dir": 2,
            "patch_bc(1)%loc": -1,
            "patch_bc(1)%geometry": 1,
            "patch_bc(1)%type": -17,
            "patch_bc(1)%centroid(1)": 0.0,
            "patch_bc(1)%length(1)": djet,
            # Formatted Database File Structures
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cf_wrt": "T",
            "lag_db_wrt": "T",
            "parallel_io": "T",
            # Inflow boundary condition
            "bc_x%vb1": velA,
            "bc_x%vb2": 0.0,
            "bc_x%vb3": 0.0,
            # Patch 1: Initial crossflow (air-dominated)
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": x_domain_beg + 0.5 * x_length,
            "patch_icpp(1)%y_centroid": y_domain_beg + 0.5 * y_length,
            "patch_icpp(1)%length_x": x_length,
            "patch_icpp(1)%length_y": y_length,
            "patch_icpp(1)%vel(1)": velA,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": pA,
            "patch_icpp(1)%alpha_rho(1)": eps * rho_hex,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%cf_val": 0,
            "patch_icpp(1)%cf_val2": 0,
            "patch_icpp(1)%cf_val3": 0,
            "patch_icpp(1)%alpha_rho(2)": eps * rho_vapor,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(1)%alpha_rho(3)": (1.0 - 2.0 * eps) * rho_air,
            "patch_icpp(1)%alpha(3)": 1.0 - 2.0 * eps,
            "patch_icpp(1)%r0": 1.0,
            "patch_icpp(1)%v0": 0.0,
            # Patch 2: Jet seeded with dense vapor bubbles
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": jet_patch_xc,
            "patch_icpp(2)%y_centroid": jet_patch_yc,
            "patch_icpp(2)%length_x": jet_patch_length_x,
            "patch_icpp(2)%length_y": djet,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": velJ,
            "patch_icpp(2)%pres": pL,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - vapor_frac - eps) * rho_hex,
            "patch_icpp(2)%alpha(1)": 1.0 - vapor_frac - eps,
            "patch_icpp(2)%alpha_rho(2)": vapor_frac * rho_vapor,
            "patch_icpp(2)%alpha(2)": vapor_frac,
            "patch_icpp(2)%alpha_rho(3)": eps * rho_air,
            "patch_icpp(2)%alpha(3)": eps,
            "patch_icpp(2)%cf_val": 1,
            "patch_icpp(2)%cf_val2": 1,
            "patch_icpp(2)%cf_val3": 0,
            "patch_icpp(2)%r0": 1.0,
            "patch_icpp(2)%v0": 0.0,
            # Fluid properties
            "fluid_pp(1)%gamma": 1.0 / (4.0 - 1.0),
            "fluid_pp(1)%pi_inf": 1.649e8,
            "fluid_pp(1)%Re(1)": 1 / mu_hex,
            "fluid_pp(2)%gamma": 1.0 / (1.09 - 1.0),
            "fluid_pp(2)%pi_inf": 0.0,
            "fluid_pp(2)%Re(1)": 1 / mu_vapor,
            "fluid_pp(3)%gamma": 1.0 / (gam_air - 1.0),
            "fluid_pp(3)%pi_inf": 0.0,
            "fluid_pp(3)%Re(1)": 1 / mu_air,
            # Surface tension --------------------------------------------------
            "sigma": sigma_hex_air,
            "sigma_2": sigma_hex_air,
            "sigma_3": sigma_hex_air,
        }
    )
)
