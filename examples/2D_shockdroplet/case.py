#!/usr/bin/env python3
import json
import math

# Liang et al. (2020) Case 1 style setup:
# planar air shock + water droplet containing a centered vapor cavity.

# -------------------------------
# Reference geometry
# -------------------------------
D0 = 2.3e-3  # droplet diameter [m]
R0 = D0 / 2
rc_over_R0 = 0.25
rc = rc_over_R0 * R0

# -------------------------------
# Thermodynamic states
# -------------------------------
# Air (fluid 3)
gam_a = 1.4
pi_a = 0.0
p_sat = 2339.0
p_a0 = p_sat
rho_a0 = 1.2 * (p_a0 / 101325.0)

# Water liquid (fluid 1) - stiffened gas approximation used in other MFC droplet examples
rho_w0 = 1000.0
gam_w = 6.12
pi_w = 3.43e8
p_w0 = p_sat

# Water vapor (fluid 2)
# Saturation pressure is prescribed for the cavity at initialization.
rho_v0 = 0.017
gam_v = 1.33
pi_v = 0.0

# Surface tension [N/m]
sigma_wg = 0.072  # water-air
sigma_wv = 0.072  # water-vapor

# -------------------------------
# Incident shock in air
# -------------------------------
Min = 1.47

ps_over_p0 = 1.0 + 2.0 * gam_a / (gam_a + 1.0) * (Min**2 - 1.0)
rhos_over_rho0 = ((gam_a + 1.0) * Min**2) / ((gam_a - 1.0) * Min**2 + 2.0)

p_as = ps_over_p0 * p_a0
rho_as = rhos_over_rho0 * rho_a0

# Shock speed and post-shock velocity in lab frame
c_a0 = math.sqrt(gam_a * (p_a0 + pi_a) / rho_a0)
Ms = math.sqrt((gam_a + 1.0) / (2.0 * gam_a) * (ps_over_p0 - 1.0) + 1.0)
ss = Ms * c_a0
u_as = (2.0 * c_a0 / (gam_a + 1.0)) * (Min**2 - 1.0) / Min
u_as = max(0.0, u_as)

# -------------------------------
# Domain and numerics
# -------------------------------
# Use 160 cells across D0 (inside target 100-200)
cells_per_D = 160
Lx = 12.0 * D0
Ly = 3.0 * D0
Nx = int(cells_per_D * Lx / D0)
Ny = int(cells_per_D * Ly / D0)

dx = Lx / Nx
cfl = 0.05

# conservative max speed estimate
c_as = math.sqrt(gam_a * (p_as + pi_a) / rho_as)
c_w = math.sqrt(gam_w * (p_w0 + pi_w) / rho_w0)
max_speed = max(ss, abs(u_as) + c_as, c_w)
dt = 0.2 * cfl * dx / max_speed

# run until shock traverses droplet and cavity collapse/jet onset window
# (about 2.5 D0 / shock speed)
t_stop = 2.5 * D0 / ss
t_save = t_stop / 100.0

eps = 1.0e-4

x_beg = -3.0 * D0
x_end = 9.0 * D0
y_beg = 0.0
y_end = 3.0 * D0

# centered droplet (axisymmetric centerline at y=0)
x0 = 0.0
y0 = 0.0

print(
    json.dumps(
        {
            # Logistics
            "run_time_info": "T",
            # Computational Domain Parameters
            "x_domain%beg": x_beg,
            "x_domain%end": x_end,
            "y_domain%beg": y_beg,
            "y_domain%end": y_end,
            "m": Nx,
            "n": Ny,
            "p": 0,
            "cyl_coord": "T",
            "dt": dt,
            "cfl_adap_dt": "F",
            "cfl_target": cfl,
            "n_start": 0,
            "t_stop": t_stop,
            "t_save": t_save,
            # Simulation algorithm
            "num_patches": 4,
            "model_eqns": 3,
            "num_fluids": 3,
            "alt_soundspeed": "F",
            "mpp_lim": "T",
            "mixture_err": "F",
            "relax": "T",
            "relax_model": 6,
            "palpha_eps": eps,
            "ptgalpha_eps": eps,
            "time_stepper": 1,
            "weno_order": 3,
            "weno_eps": 1.0e-12,
            "mapped_weno": "F",
            "riemann_solver": 2,
            "wave_speeds": 1,
            "avg_state": 2,
            "surface_tension": "T",
            # BCs
            "bc_x%beg": -6,
            "bc_x%end": -6,
            "bc_y%beg": -2,
            "bc_y%end": -6,
            # Output
            "format": 1,
            "precision": 2,
            "prim_vars_wrt": "T",
            "cons_vars_wrt": "T",
            "alpha_wrt": "T",
            "parallel_io": "T",
            # Patch 1: pre-shock air background
            "patch_icpp(1)%geometry": 3,
            "patch_icpp(1)%x_centroid": 0.5 * (x_beg + x_end),
            "patch_icpp(1)%y_centroid": 0.5 * (y_beg + y_end),
            "patch_icpp(1)%length_x": x_end - x_beg,
            "patch_icpp(1)%length_y": y_end - y_beg,
            "patch_icpp(1)%vel(1)": 0.0,
            "patch_icpp(1)%vel(2)": 0.0,
            "patch_icpp(1)%pres": p_a0,
            "patch_icpp(1)%alpha_rho(1)": eps * rho_w0,
            "patch_icpp(1)%alpha_rho(2)": eps * rho_v0,
            "patch_icpp(1)%alpha_rho(3)": (1.0 - 2.0 * eps) * rho_a0,
            "patch_icpp(1)%alpha(1)": eps,
            "patch_icpp(1)%alpha(2)": eps,
            "patch_icpp(1)%alpha(3)": 1.0 - 2.0 * eps,
            # Patch 2: shocked air slab initialized upstream
            "patch_icpp(2)%geometry": 3,
            "patch_icpp(2)%alter_patch(1)": "T",
            "patch_icpp(2)%x_centroid": -2.0 * D0,
            "patch_icpp(2)%y_centroid": 0.5 * (y_beg + y_end),
            "patch_icpp(2)%length_x": 2.0 * D0,
            "patch_icpp(2)%length_y": y_end - y_beg,
            "patch_icpp(2)%vel(1)": u_as,
            "patch_icpp(2)%vel(2)": 0.0,
            "patch_icpp(2)%pres": p_as,
            "patch_icpp(2)%alpha_rho(1)": eps * rho_w0,
            "patch_icpp(2)%alpha_rho(2)": eps * rho_v0,
            "patch_icpp(2)%alpha_rho(3)": (1.0 - 2.0 * eps) * rho_as,
            "patch_icpp(2)%alpha(1)": eps,
            "patch_icpp(2)%alpha(2)": eps,
            "patch_icpp(2)%alpha(3)": 1.0 - 2.0 * eps,
            # Patch 3: liquid droplet
            "patch_icpp(3)%geometry": 2,
            "patch_icpp(3)%alter_patch(1)": "T",
            "patch_icpp(3)%x_centroid": x0,
            "patch_icpp(3)%y_centroid": y0,
            "patch_icpp(3)%radius": R0,
            "patch_icpp(3)%vel(1)": 0.0,
            "patch_icpp(3)%vel(2)": 0.0,
            "patch_icpp(3)%pres": p_w0,
            "patch_icpp(3)%alpha_rho(1)": (1.0 - 2.0 * eps) * rho_w0,
            "patch_icpp(3)%alpha_rho(2)": eps * rho_v0,
            "patch_icpp(3)%alpha_rho(3)": eps * rho_a0,
            "patch_icpp(3)%alpha(1)": 1.0 - 2.0 * eps,
            "patch_icpp(3)%alpha(2)": eps,
            "patch_icpp(3)%alpha(3)": eps,
            # Patch 4: centered vapor cavity inside droplet
            "patch_icpp(4)%geometry": 2,
            "patch_icpp(4)%alter_patch(3)": "T",
            "patch_icpp(4)%x_centroid": x0,
            "patch_icpp(4)%y_centroid": y0,
            "patch_icpp(4)%radius": rc,
            "patch_icpp(4)%vel(1)": 0.0,
            "patch_icpp(4)%vel(2)": 0.0,
            "patch_icpp(4)%pres": p_sat,
            "patch_icpp(4)%alpha_rho(1)": eps * rho_w0,
            "patch_icpp(4)%alpha_rho(2)": (1.0 - 2.0 * eps) * rho_v0,
            "patch_icpp(4)%alpha_rho(3)": eps * rho_a0,
            "patch_icpp(4)%alpha(1)": eps,
            "patch_icpp(4)%alpha(2)": 1.0 - 2.0 * eps,
            "patch_icpp(4)%alpha(3)": eps,
            # Fluid properties
            "fluid_pp(1)%gamma": 1.0 / (gam_w - 1.0),
            "fluid_pp(1)%pi_inf": gam_w * pi_w / (gam_w - 1.0),
            "fluid_pp(2)%gamma": 1.0 / (gam_v - 1.0),
            "fluid_pp(2)%pi_inf": gam_v * pi_v / (gam_v - 1.0),
            "fluid_pp(3)%gamma": 1.0 / (gam_a - 1.0),
            "fluid_pp(3)%pi_inf": gam_a * pi_a / (gam_a - 1.0),
            # Surface tensions (fluid 1 with fluids 2 and 3)
            "sigma": sigma_wv,
            "sigma_2": sigma_wg,
            "sigma_3": 1.0e-6,
            # here sigma_3 represents near-zero surface tension between fluids 2 and 3
        }
    )
)
