import pathlib
from collections import namedtuple

import numpy as np
from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('TkAgg')

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

Simulation = namedtuple("Simulation", ["gradation", "results_directory", "color"])


def get_particle_volume(distribution_file):
    radii = np.genfromtxt(distribution_file, dtype=float)
    return 4*np.pi*np.sum(radii**3)/3


def fines_model(D, radius, height, D_th, d_boundary, K, nl):
    D0 = (1-2*(height + radius)/radius/height*d_boundary)*D_th
    D0[D0 < 0] = 0
    p = K*(D - D0)**nl
    p[D - D0 < 0] = 0
    return p


def residual_D_th(D_th, data, material_parameters):
    D_th = abs(D_th)
    res = 0
    for (D, p_sim, r, h, _) in data:
        p_model = fines_model(D, r, h, D_th, 0, material_parameters[1], material_parameters[2])
        res += np.sum((p_sim - p_model)**2)
    return res


def residual_common(par, data):
    par = abs(par)
    res = 0
    for (D, p_sim, r, h, _) in data:
        max_idx = np.argmax(p_sim)
        rel_d = D - D[max_idx] + 0.025
        p_model = fines_model(rel_d, r, h, par[0], 0, par[1], par[2])
        res += np.sum((p_sim - p_model)**2)
    return res


def main():
    gradation_dir = pathlib.Path("~/DEMsim/size_distributions/kcv").expanduser()
    results_dir = pathlib.Path("~/DEMresults/kcv/monosize").expanduser()
    simulations_mono = [
        Simulation(gradation=gradation_dir / "mono/kcv_2mm_0_1kg.dat",
                   results_directory=results_dir / "0_1kg_2mm/fines/1/compaction", color="b"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_1kg.dat",
                   results_directory=results_dir/"0_1kg_2mm/fines/2/compaction", color="b"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_1kg.dat",
                   results_directory=results_dir/"0_1kg_2mm/fines/3/compaction", color="b"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_2kg.dat",
                   results_directory=results_dir/"0_2kg_2mm/fines/1/compaction", color="r"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_2kg.dat",
                   results_directory=results_dir/"0_2kg_2mm/fines/2/compaction", color="r"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_2kg.dat",
                   results_directory=results_dir/"0_2kg_2mm/fines/3/compaction", color="r"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_5kg.dat",
                   results_directory=results_dir/"0_5kg_2mm/fines/1/compaction", color="g"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_5kg.dat",
                  results_directory=results_dir/"0_5kg_2mm/fines/2/compaction", color="g"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_5kg.dat",
                   results_directory=results_dir/"0_5kg_2mm/fines/3/compaction", color="g"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_5kg.dat",
                   results_directory=results_dir/"0_5kg_2mm/fines/4/compaction", color="g"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_0_5kg.dat",
                   results_directory=results_dir/"0_5kg_2mm/fines/5/compaction", color="g"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_1_0kg.dat",
                   results_directory=results_dir/"1_0kg_2mm/fines/1/compaction", color="k"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_1_0kg.dat",
                   results_directory=results_dir/"1_0kg_2mm/fines/2/compaction", color="k"),
        Simulation(gradation=gradation_dir/"mono/kcv_2mm_1_0kg.dat",
                   results_directory=results_dir/"1_0kg_2mm/fines/3/compaction", color="k"),
    ]
    data_sets = []
    for sim in simulations_mono:
        particle_volume = get_particle_volume(sim.gradation)
        surface_position_file = sim.results_directory / "surface_positions.dou"
        surface_data = np.genfromtxt(surface_position_file, delimiter=",")
        radius = surface_data[0, 2]
        height = surface_data[:, -2]
        cylinder_volume = np.pi*radius**2*height
        density = particle_volume/cylinder_volume
        idx = np.argmax(density)
        surface_force_file = sim.results_directory /"surface_forces.dou"
        surface_forces = np.genfromtxt(surface_force_file, delimiter=",")
        compaction_pressure = surface_forces[:, -2]/np.pi/radius**2/1e6
        plt.figure(0)
        plt.plot(density, compaction_pressure, sim.color, lw=2)
        data_sets.append((density[:idx], compaction_pressure[:idx], radius, height[:idx], sim.color))

        max_idx = np.argmax(compaction_pressure)
        plt.figure(1)
        plt.plot(density - density[max_idx] + 0.025, compaction_pressure, sim.color, lw=2)

    parameters_old = np.array([[0.005, 3166059e6, 4.2]])
    parameters_new = 0*parameters_old
    res = 10
    while res > 0:
        parameters_new = fmin(residual_common, parameters_old, args=(data_sets,), maxfun=1e4, maxiter=1e4)
        res = np.sum(abs(parameters_old - parameters_new))
        parameters_old = parameters_new
        print("\tsuggested parameters", parameters_new)
    parameters_new = abs(parameters_new)
    print(parameters_new)

    D1 = np.linspace(0, 0.025, 1000)
    p_model = fines_model(D1, 1 + 0*D1, 1 + 0*D1, parameters_new[0], 0, parameters_new[1], parameters_new[2])

    plt.plot(D1, p_model, '--k', lw=4)

    boundary_ratio = np.zeros(len(data_sets))
    start_D = 0*boundary_ratio
    for i, data_set in enumerate(data_sets):
        density = data_set[0]
        r = data_set[2]
        height = data_set[3]
        c = data_set[4]
        boundary_ratio[i] = 2*(height[0] + r)/r/height[0]
        start_D[i] = parameters_new[0] - 0.025 + density[-1]
        plt.figure(2)

        plt.plot(boundary_ratio[i], start_D[i], c + 'x', mew=2, ms=12)
    db, D0 = np.polyfit(boundary_ratio, start_D, 1)
    db /= -D0
    print(D0, db)
    plt.plot(np.linspace(40, 150), D0 - db*np.linspace(40, 150), '--k', lw=4)
    """
    for (D_sim, _, r, h, c) in data_sets:
        p_model = fines_model(D_sim, r, h, *parameters_new)
        plt.plot(D_sim, p_model, '--' + c, lw=2)
    """

    for data_set in data_sets:
        density = data_set[0]
        r = data_set[2]
        height = data_set[3]
        c = data_set[4]
        p = fines_model(density, r, height, D0, db, parameters_new[1], parameters_new[2])
        plt.figure(0)
        plt.plot(density, p, '--' + c, lw=2)

    plt.figure(0)
    plt.ylim(0, 8)
    plt.xlabel("Relative density [-]")
    plt.ylabel("Pressure [MPa]")
    plt.tight_layout()
    plt.savefig("fines_all.png")
    plt.show()

if __name__ == '__main__':
    main()