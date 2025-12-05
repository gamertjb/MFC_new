#!/usr/bin/env python3
"""Lagrangian vapor bubble carried by a liquid jet into a subsonic crossflow."""

import json
import math

# -----------------------------------------------------------------------------
# Reference scales
# -----------------------------------------------------------------------------
djet = 2.5e-4  # m
x0 = djet
rho0 = 660.0  # kg/m^3 (liquid reference density)
c0 = 1200.0  # m/s (liquid reference sound speed)
p0 = rho0 * c0 * c0
T0 = 300.0  # K

# -----------------------------------------------------------------------------
# Fluid properties
# -----------------------------------------------------------------------------
rho_liq = rho0
rho_vap = 2.5
rho_air = 1.2

gamma_liq = 4.0
pi_inf_liq = 1.65e8
mu_liq = 0.30e-3

sigma_liq = 0.018  # N/m, used by the Lagrangian microphysics
pv_liq = 2500.0  # Pa

gamma_air = 1.4
mu_air = 1.8e-5

gamma_vap = 1.09
pi_inf_vap = 0.0
mu_vap = 1.5e-5

# Gas/vapor transport properties for the Lagrangian model
cp_g = 1.0e3
cp_v = 2.1e3
k_g = 0.025
k_v = 0.02
MW_g = 28.0
MW_v = 18.0

diff_vapor = 2.5e-5

# -----------------------------------------------------------------------------
# Flow configuration
# -----------------------------------------------------------------------------
u_cross = 60.0  # m/s subsonic crossflow velocity (x-direction)
u_jet = 80.0  # m/s jet injection speed (y-direction)

x_beg = -6.0 * djet
x_end = 18.0 * djet
y_beg = 0.0
y_end = 15.0 * djet

Lx = x_end - x_beg
Ly = y_end - y_beg

# -----------------------------------------------------------------------------
# Grid resolution and time stepping
# -----------------------------------------------------------------------------
dx_target = 0.08 * djet
Nx = max(1, int(round(Lx / dx_target)))
Ny = max(1, int(round(Ly / dx_target)))

dx = Lx / Nx

cfl = 0.5
# nondimensionalize dt using x0/c0 scaling
dt = cfl * (dx / c0) * (c0 / x0)

# end time and output cadence (nondimensional)
t_stop = 3.0e-6 * c0 / x0  # 3 microseconds of physical time
t_save = t_stop / 50.0

# -----------------------------------------------------------------------------
# Geometry helpers
# -----------------------------------------------------------------------------
jet_length_x = djet
jet_length_y = 4.0 * djet
jet_center_x = 0.0
jet_center_y = y_beg + 0.5 * jet_length_y

bubble_entry_x = jet_center_x / x0
bubble_entry_y = (jet_center_y + 0.25 * jet_length_y) / x0
bubble_radius = 0.10  # fraction of jet diameter (already nondimensional)

# -----------------------------------------------------------------------------
# Patch volume fractions (nondimensionalized by rho0)
# -----------------------------------------------------------------------------
_eps = 1.0e-6
alpha_liq_cross = _eps
alpha_vap_cross = _eps
alpha_air_cross = 1.0 - 2.0 * _eps

alpha_liq_jet = 1.0 - 2.0 * _eps
alpha_vap_jet = _eps
alpha_air_jet = _eps

print(
    json.dumps(
        {
            # Logistics -------------------------------------------------------
            "run_time_info": "T",
            "x_domain%beg": x_beg / x0,
            "x_domain%end": x_end / x0,
            "y_domain%beg": y_beg / x0,
            "y_domain%end": y_end / x0,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "dt": dt,
            "cfl_adap_dt": "T",
            "cfl_target": cfl,
            "n_start": 0,
            "t_stop": t_stop,
            "t_save": t_save,
            # Numerics --------------------------------------------------------
            "num_patches": 2,
            "model_eqns": 3,
            "num_fluids": 3,
            "alt_soundspeed": "F",
            "mpp_lim": "T",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "viscous": "T",
            "surface_tension": "T",
            "weno_order": 5,
            "weno_eps": 1.0e-16,
            "mapped_weno": "T",
            "weno_Re_flux": "T",
            "weno_avg": "T",
            "null_weights": "T",
            "mp_weno": "F",
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
            "patch_bc(1)%centroid(1)": jet_center_x / x0,
            "patch_bc(1)%length(1)": jet_length_x / x0,
            "bc_x%vb1": u_cross / c0,
            "bc_x%vb2": 0.0,
            # Output ----------------------------------------------------------
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "parallel_io": "T",
            "lag_db_wrt": "T",
            # Patch 1: crossflow (mostly air) ---------------------------------
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * (x_beg + x_end) / x0,
            "patch_icpp(1)%y_centroid": 0.5 * (y_beg + y_end) / x0,
            "patch_icpp(1)%length_x": (x_end - x_beg) * 10.0 / x0,
            "patch_icpp(1)%length_y": (y_end - y_beg) * 10.0 / x0,
            "patch_icpp(1)%vel(1)": u_cross / c0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": 101325.0 / p0,
            "patch_icpp(1)%alpha_rho(1)": alpha_liq_cross * rho_liq / rho0,
            "patch_icpp(1)%alpha_rho(2)": alpha_vap_cross * rho_vap / rho0,
            "patch_icpp(1)%alpha_rho(3)": alpha_air_cross * rho_air / rho0,
            "patch_icpp(1)%alpha(1)": alpha_liq_cross,
            "patch_icpp(1)%alpha(2)": alpha_vap_cross,
            "patch_icpp(1)%alpha(3)": alpha_air_cross,
            "patch_icpp(1)%cf_val": 0,
            "patch_icpp(1)%cf_val2": 0,
            # Patch 2: liquid jet carrying the seeded bubble ------------------
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": jet_center_x / x0,
            "patch_icpp(2)%y_centroid": jet_center_y / x0,
            "patch_icpp(2)%length_x": jet_length_x / x0,
            "patch_icpp(2)%length_y": jet_length_y / x0,
            "patch_icpp(2)%vel(1)": 0.0,
            "patch_icpp(2)%vel(2)": u_jet / c0,
            "patch_icpp(2)%pres": 101325.0 / p0,
            "patch_icpp(2)%alpha_rho(1)": alpha_liq_jet * rho_liq / rho0,
            "patch_icpp(2)%alpha_rho(2)": alpha_vap_jet * rho_vap / rho0,
            "patch_icpp(2)%alpha_rho(3)": alpha_air_jet * rho_air / rho0,
            "patch_icpp(2)%alpha(1)": alpha_liq_jet,
            "patch_icpp(2)%alpha(2)": alpha_vap_jet,
            "patch_icpp(2)%alpha(3)": alpha_air_jet,
            "patch_icpp(2)%cf_val": 1,
            "patch_icpp(2)%cf_val2": 1,
            # Lagrangian Bubbles ----------------------------------------------
            "bubbles_lagrange": "T",
            "bubble_model": 2,
            "lag_params%nBubs_glb": 1,
            "lag_params%solver_approach": 2,
            "lag_params%cluster_type": 2,
            "lag_params%pressure_corrector": "T",
            "lag_params%smooth_type": 1,
            "lag_params%heatTransfer_model": "T",
            "lag_params%massTransfer_model": "T",
            "lag_params%epsilonb": 1.0,
            "lag_params%valmaxvoid": 0.9,
            "lag_params%write_bubbles": "F",
            "lag_params%write_bubbles_stats": "F",
            "lag_params%c0": c0,
            "lag_params%rho0": rho0,
            "lag_params%T0": T0,
            "lag_params%x0": x0,
            "lag_params%diffcoefvap": diff_vapor,
            "lag_params%Thost": T0,
            "lag_params%charwidth": 1.0,
            # Fluids Physical Parameters --------------------------------------
            # Liquid jet (host medium)
            "fluid_pp(1)%gamma": 1.0 / (gamma_liq - 1.0),
            "fluid_pp(1)%pi_inf": gamma_liq * (pi_inf_liq / p0) / (gamma_liq - 1.0),
            "fluid_pp(1)%Re(1)": 1.0 / (mu_liq / (rho0 * c0 * x0)),
            "fluid_pp(1)%mul0": mu_liq,
            "fluid_pp(1)%ss": sigma_liq,
            "fluid_pp(1)%pv": pv_liq,
            "fluid_pp(1)%gamma_v": gamma_liq,
            "fluid_pp(1)%M_v": MW_v,
            "fluid_pp(1)%k_v": k_v,
            "fluid_pp(1)%cp_v": cp_v,
            # Vapor properties for three-fluid surface tension
            "fluid_pp(2)%gamma": 1.0 / (gamma_vap - 1.0),
            "fluid_pp(2)%pi_inf": pi_inf_vap,
            "fluid_pp(2)%Re(1)": 1.0 / (mu_vap / (rho0 * c0 * x0)),
            "fluid_pp(2)%gamma_v": gamma_vap,
            "fluid_pp(2)%M_v": MW_v,
            "fluid_pp(2)%k_v": k_v,
            "fluid_pp(2)%cp_v": cp_v,
            # Air crossflow / bubble gas state
            "fluid_pp(3)%gamma": 1.0 / (gamma_air - 1.0),
            "fluid_pp(3)%pi_inf": 0.0,
            "fluid_pp(3)%Re(1)": 1.0 / (mu_air / (rho0 * c0 * x0)),
            "fluid_pp(3)%gamma_v": gamma_air,
            "fluid_pp(3)%M_v": MW_g,
            "fluid_pp(3)%k_v": k_g,
            "fluid_pp(3)%cp_v": cp_g,
            # Surface tension coefficients for interfaces 1-2 and 1-3
            "sigma": sigma_liq,
            "sigma_2": sigma_liq,
        }
    )
)
