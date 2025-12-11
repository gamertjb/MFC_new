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
R_air = 287.0
k_air = 0.0262
cpA = cvA + R_air
Mv_air = 28.97e-3

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
D_v = 1.0e-5

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
cfl = 0.95
time_end = 3.0e-4
dt = 5.0 * cfl * dx / c_air

# Output cadence expressed in physical time so adaptive stepping still saves often
num_output_frames = 100
t_save = time_end / num_output_frames

# Regularisation and surface tension
_eps = 1.0e-6
sigma_hex_air = 0.018
sigma_hex_vapor = 0.018

# Convenient geometry helpers
crossflow_center_x = 0.5 * (x_beg + x_end)
crossflow_center_y = 0.5 * (y_beg + y_end)
crossflow_length_x = x_end - x_beg
crossflow_length_y = y_end - y_beg

jet_length_x = djet
jet_length_y = djet
jet_center_x = 0.0
# Keep the jet patch comfortably away from the domain boundaries so that the
# vapor core (bubble) starts fully contained inside the liquid jet and does not
# touch any walls.
jet_margin = 0.1 * jet_length_y
jet_center_y = y_beg + 0.5 * jet_length_y + jet_margin

# Reference scales for the Lagrangian bubble model
c0_ref = math.sqrt((pi_inf_L + (1.0 / gammaL_field + 1.0) * pA) / rhoL)
rho0_ref = rhoL
T0_ref = 300.0
x0_ref = djet

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
            "num_patches": 2,
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
            "lag_db_wrt": "T",
            # Lagrangian bubble model ----------------------------------------
            "bubbles_lagrange": "T",
            "bubble_model": 2,
            "lag_params%nBubs_glb": 1,
            "lag_params%solver_approach": 2,
            "lag_params%cluster_type": 2,
            "lag_params%pressure_corrector": "T",
            "lag_params%smooth_type": 1,
            "lag_params%heatTransfer_model": "T",
            "lag_params%massTransfer_model": "F",
            "lag_params%epsilonb": 1.0,
            "lag_params%charwidth": djet,
            "lag_params%valmaxvoid": 0.9,
            "lag_params%c0": c0_ref,
            "lag_params%rho0": rho0_ref,
            "lag_params%T0": T0_ref,
            "lag_params%x0": x0_ref,
            "lag_params%Thost": T0_ref,
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
            "patch_icpp(2)%cf_val2": 1,
            # Fluid properties -------------------------------------------------
            "fluid_pp(1)%gamma": gammaL_field,
            "fluid_pp(1)%pi_inf": pi_inf_L,
            "fluid_pp(1)%cv": cvL,
            "fluid_pp(1)%qv": qv_L,
            "fluid_pp(1)%qvp": qvp_L,
            "fluid_pp(1)%Re(1)": 1.0 / muL,
            "fluid_pp(1)%mul0": muL,
            "fluid_pp(1)%ss": sigma_hex_air,
            "fluid_pp(1)%pv": 2.0e4,
            "fluid_pp(1)%gamma_v": gam_v,
            "fluid_pp(1)%M_v": 86.18e-3,
            "fluid_pp(1)%mu_v": muV,
            "fluid_pp(1)%k_v": 0.014,
            "fluid_pp(1)%cp_v": cvV + 96.5,
            "fluid_pp(1)%D_v": D_v,
            "fluid_pp(2)%gamma": 1.0 / (gamma_air - 1.0),
            "fluid_pp(2)%pi_inf": pi_inf_A,
            "fluid_pp(2)%cv": cvA,
            "fluid_pp(2)%qv": qv_A,
            "fluid_pp(2)%qvp": qvp_A,
            "fluid_pp(2)%Re(1)": 1.0 / muA,
            "fluid_pp(2)%gamma_v": gamma_air,
            "fluid_pp(2)%M_v": Mv_air,
            "fluid_pp(2)%mu_v": muA,
            "fluid_pp(2)%k_v": k_air,
            "fluid_pp(2)%cp_v": cpA,
            "fluid_pp(2)%D_v": 2.0e-5,
            "fluid_pp(3)%gamma": gammaV_field,
            "fluid_pp(3)%pi_inf": pi_inf_V,
            "fluid_pp(3)%cv": cvV,
            "fluid_pp(3)%qv": qv_V,
            "fluid_pp(3)%qvp": qvp_V,
            "fluid_pp(3)%Re(1)": 1.0 / muV,
            "fluid_pp(3)%gamma_v": gam_v,
            "fluid_pp(3)%M_v": 86.18e-3,
            "fluid_pp(3)%mu_v": muV,
            "fluid_pp(3)%k_v": 0.014,
            "fluid_pp(3)%cp_v": cvV + 96.5,
            "fluid_pp(3)%D_v": D_v,
            # Surface tension --------------------------------------------------
            "sigma": sigma_hex_air,
            "sigma_2": sigma_hex_vapor,
        }
    )
)
