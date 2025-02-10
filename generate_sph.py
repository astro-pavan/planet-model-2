from planet import planet
import woma
import numpy as np
import os
import shutil

R_earth = 6.371e6  # m
M_earth = 5.9724e24  # kg

m = 5.6
m_water = m * 0.2
m_mantle = (m - m_water) * 0.66
m_core = m - m_mantle - m_water

sub_neptune = planet('sub-neptune', m_core, m_mantle, m_water, 0, 226)

n_particles = 1e6

mars = woma.Planet(
    name='Mars',
    A1_mat_layer=["ANEOS_iron", "ANEOS_forsterite"],
    A1_T_rho_type=["adiabatic", "adiabatic"],
    P_s=1e5,
    T_s=300,
    A1_R_layer=[0.25 * R_earth, 0.5 * R_earth],
    verbosity=0
)

mars.gen_prof_L2_find_M_given_R_R1(M_max=1 * M_earth, verbosity=0)
impactor = mars

M_total = sub_neptune.M + impactor.M
n_particles_target = int(n_particles * (sub_neptune.M / M_total))
n_particles_impactor = int(n_particles * (impactor.M / M_total))

target_particles, target_A1_h, target_A2_pos = sub_neptune.make_particle_planet(n_particles)
impactor_particles = woma.ParticlePlanet(impactor, n_particles_impactor, verbosity=1)

A1_pos_i, A1_vel_i = woma.impact_pos_vel_b_v_c_t(
    b=0.3,
    v_c=1.1,
    t=1800,
    R_t=sub_neptune.R,
    R_i=impactor.R,
    M_t=sub_neptune.M,
    M_i=impactor.M,
    units_v_c="v_esc",
    units_b="b"
)

impactor_particles.A2_pos += A1_pos_i
impactor_particles.A2_vel += A1_vel_i

A1_COM = np.sum(target_A2_pos.T * target_particles.A1_m, axis=1) / np.sum(target_particles.A1_m)

pos = np.array(target_A2_pos) - A1_COM

A1_vel_COM = (sub_neptune.M * A1_vel_i) / M_total

target_particles_A2_vel = np.zeros_like(target_A2_pos) - A1_vel_COM
impactor_particles.A2_vel -= A1_vel_COM

target_particles.A1_mat[target_particles.A1_mat == 304] = 901

os.mkdir('../impact-3')
os.mkdir('../impact-3/snapshots')
shutil.copyfile('impact-parameters.yml', '../impact-3/impact-parameters.yml')
shutil.copyfile('slurm-submit.sh', '../impact-3/slurm-submit.sh')

import h5py
with h5py.File(f'../impact-3/initial_conditions.hdf5', "w") as f:
    woma.save_particle_data(
        f,
        np.append(target_A2_pos, impactor_particles.A2_pos, axis=0),
        np.append(target_particles_A2_vel, impactor_particles.A2_vel, axis=0),
        np.append(target_particles.A1_m, impactor_particles.A1_m),
        np.append(target_A1_h, impactor_particles.A1_h),
        np.append(target_particles.A1_rho, impactor_particles.A1_rho),
        np.append(target_particles.A1_P, impactor_particles.A1_P),
        np.append(target_particles.A1_u, impactor_particles.A1_u),
        np.append(target_particles.A1_mat, impactor_particles.A1_mat_id),
        boxsize=100 * R_earth,
        file_to_SI=woma.Conversions(M_earth, R_earth, 1),
    )