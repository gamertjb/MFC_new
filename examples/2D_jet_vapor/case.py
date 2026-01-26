#!/usr/bin/env python3
"""Input deck for a vaporizing vertical jet in a crossflow."""

import json
import math

# -----------------------------------------------------------------------------
# Gas (air) properties
# -----------------------------------------------------------------------------
pA = 101325.0
rhoA = 0.87
gamma_air = 1.4
muA = 2.30e-5
cvA = 717.5
pi_inf_A = 0.0
qv_A = 0.0
qvp_A = 0.0

c_air = math.sqrt(gamma_air * pA / rhoA)

# -----------------------------------------------------------------------------
# Liquid (hexane) properties
# -----------------------------------------------------------------------------
rhoL = 660.0
muL = 0.30e-3
cvL = 2200.0
pi_inf_L = 1.649e8
gammaL_field = 1.0 / (4.0 - 1.0)
qv_L = -3.3e5
qvp_L = 0.0

# -----------------------------------------------------------------------------
# Vapor properties
# -----------------------------------------------------------------------------
rhoV = 2.5
muV = 1.5e-5
gam_v = 1.09
cvV = 1600.0
pi_inf_V = 0.0
qv_V = 3.3e5
qvp_V = 0.0

gammaV_field = 1.0 / (gam_v - 1.0)

# -----------------------------------------------------------------------------
# Flow configuration
# -----------------------------------------------------------------------------
djet = 0.259e-3
uA = 60.0
qFA = 40.0
uJ = uA * math.sqrt(qFA * rhoA / rhoL)

# Domain extents (trimmed version of the original case)
x_beg = -5.0 * djet
x_end = 25.0 * djet
y_beg = 0.0
y_end = 20.0 * djet

Lx = x_end - x_beg
Ly = y_end - y_beg

# -----------------------------------------------------------------------------
# Grid resolution
# Goal: keep approximately the same physical dx (~4.0 microns) as the reference
# -----------------------------------------------------------------------------
dx_target = 4.0e-6  # m, from the original high-res case
Nx = max(1, int(round(Lx / dx_target)))
Ny = max(1, int(round(Ly / dx_target)))

dx = Lx / Nx  # actual realized spacing after rounding

# Time stepping
cfl = 0.2
time_end = 1.0e-5
dt = cfl * dx / c_air

# Output cadence expressed in physical time so adaptive stepping still saves often
num_output_frames = 100
t_save = time_end / num_output_frames

# Regularisation and surface tension
_eps = 1.0e-5
sigma_hex_air = 0.018
sigma_hex_vapor = 0.018
sigma_air_vapor = 0.0

# Convenient geometry helpers
crossflow_center_x = 0.5 * (x_beg + x_end)
crossflow_center_y = 0.5 * (y_beg + y_end)
crossflow_length_x = x_end - x_beg
crossflow_length_y = y_end - y_beg

jet_length_x = djet
jet_length_y = 1.5 * djet
jet_center_x = 0.0
jet_center_y = y_beg + 0.5 * jet_length_y

# Lagrangian bubble sizing (increase to make the initial bubble larger)
bubble_charwidth = 2.0 * djet

# Eulerian bubble sizing for the initial vapor bubble in the jet
bubble_radius = 0.05 * djet
bubble_center_x = jet_center_x
bubble_center_y = jet_center_y

print(
    json.dumps(
        {
            # Logistics -------------------------------------------------------
            "run_time_info": "T",
            "x_domain%beg": x_beg,
            "x_domain%end": x_end,
            "y_domain%beg": y_beg,
            "y_domain%end": y_end,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "cfl_adap_dt": "T",
            "n_start": 0,
            #"t_step_start": 0,
            "t_stop": time_end,
            "t_save": t_save,
            "cfl_target": cfl,
            # Numerics --------------------------------------------------------
            "num_patches": 3,
            "model_eqns": 3,
            "alt_soundspeed": "F",
            "num_fluids": 3,
            "mpp_lim": "T",
            "mixture_err": "T",
            "relax": "T",
            "relax_model": 6,
            "palpha_eps": _eps,
            "ptgalpha_eps": _eps,
            "time_stepper": 3,
            "weno_order": 1,
            "weno_eps": 1.0e-16,
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "mapped_weno": "F",
            "null_weights": "T",
            "mp_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "surface_tension": "T",
            "viscous": "T",
            "elliptic_smoothing": "F",
            #"elliptic_smoothing_iters": 50,
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
            "patch_bc(1)%centroid(1)": jet_center_x,
            "patch_bc(1)%length(1)": jet_length_x,
            "bc_x%vb1": uA,
            "bc_x%vb2": 0.0,
            "bc_x%vb3": 0.0,
            # Output ----------------------------------------------------------
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cf_wrt": "T",
            "parallel_io": "T",
            "lag_params%charwidth": bubble_charwidth,
            # Phase/alpha ordering: alpha(1)=liquid hexane, alpha(2)=air, alpha(3)=hexane vapor
            # Color functions for surface tension:
            #   cf_val  -> liquid/air interface (sigma = sigma_hex_air)
            #   cf_val2 -> liquid/vapor interface (sigma_2 = sigma_hex_vapor)
            #   cf_val3 -> air/vapor interface (sigma_3 = sigma_air_vapor)
            # Patch 1: crossflow (mostly air) ---------------------------------
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": crossflow_center_x,
            "patch_icpp(1)%y_centroid": crossflow_center_y,
            "patch_icpp(1)%length_x": crossflow_length_x * 10.0,
            "patch_icpp(1)%length_y": crossflow_length_y * 10.0,
            "patch_icpp(1)%vel(1)": uA,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": pA,
            "patch_icpp(1)%alpha_rho(1)": _eps * rhoL,
            "patch_icpp(1)%alpha_rho(2)": (1.0 - 2.0 * _eps) * rhoA,
            "patch_icpp(1)%alpha_rho(3)": _eps * rhoV,
            "patch_icpp(1)%alpha(1)": _eps,
            "patch_icpp(1)%alpha(2)": 1.0 - 2.0 * _eps,
            "patch_icpp(1)%alpha(3)": _eps,
            "patch_icpp(1)%cf_val": 0,
            "patch_icpp(1)%cf_val2": 0,
            "patch_icpp(1)%cf_val3": 0,
            # Patch 2: liquid jet with a small vapor core ---------------------
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": jet_center_x,
            "patch_icpp(2)%y_centroid": jet_center_y,
            "patch_icpp(2)%length_x": jet_length_x,
            "patch_icpp(2)%length_y": jet_length_y,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": uJ,
            "patch_icpp(2)%pres": pA,
            "patch_icpp(2)%alpha_rho(1)": (1.0 - 2.0 * _eps) * rhoL,
            "patch_icpp(2)%alpha_rho(2)": _eps * rhoA,
            "patch_icpp(2)%alpha_rho(3)": _eps * rhoV,
            "patch_icpp(2)%alpha(1)": 1.0 - 2.0 * _eps,
            "patch_icpp(2)%alpha(2)": _eps,
            "patch_icpp(2)%alpha(3)": _eps,
            "patch_icpp(2)%cf_val": 1,
            "patch_icpp(2)%cf_val2": 0,
            "patch_icpp(2)%cf_val3": 0,
            # Patch 3: vapor bubble inside the liquid jet ---------------------
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%alter_patch(2)": "T",
            "patch_icpp(3)%x_centroid": bubble_center_x,
            "patch_icpp(3)%y_centroid": bubble_center_y,
            "patch_icpp(3)%radius": bubble_radius,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": uJ,
            "patch_icpp(3)%pres": pA,
            "patch_icpp(3)%alpha_rho(1)": 0.0,
            "patch_icpp(3)%alpha_rho(2)": 0.0,
            "patch_icpp(3)%alpha_rho(3)": rhoV,
            "patch_icpp(3)%alpha(1)": 0.0,
            "patch_icpp(3)%alpha(2)": 0.0,
            "patch_icpp(3)%alpha(3)": 1.0,
            "patch_icpp(3)%cf_val": 0,
            "patch_icpp(3)%cf_val2": 1,
            "patch_icpp(3)%cf_val3": 1,
            # Fluid properties -------------------------------------------------
            "fluid_pp(1)%gamma": gammaL_field,
            "fluid_pp(1)%pi_inf": pi_inf_L,
            "fluid_pp(1)%cv": cvL,
            "fluid_pp(1)%qv": qv_L,
            "fluid_pp(1)%qvp": qvp_L,
            "fluid_pp(1)%Re(1)": 1.0 / muL,
            "fluid_pp(2)%gamma": 1.0 / (gamma_air - 1.0),
            "fluid_pp(2)%pi_inf": pi_inf_A,
            "fluid_pp(2)%cv": cvA,
            "fluid_pp(2)%qv": qv_A,
            "fluid_pp(2)%qvp": qvp_A,
            "fluid_pp(2)%Re(1)": 1.0 / muA,
            "fluid_pp(3)%gamma": gammaV_field,
            "fluid_pp(3)%pi_inf": pi_inf_V,
            "fluid_pp(3)%cv": cvV,
            "fluid_pp(3)%qv": qv_V,
            "fluid_pp(3)%qvp": qvp_V,
            "fluid_pp(3)%Re(1)": 1.0 / muV,
            # Surface tension --------------------------------------------------
            "sigma": sigma_hex_air,
            "sigma_2": sigma_hex_vapor,
            "sigma_3": sigma_air_vapor,
        }
    )
)
