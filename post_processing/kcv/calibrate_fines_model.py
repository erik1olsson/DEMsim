import pathlib
from collections import namedtuple

import numpy as np
from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib

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


def residual(par, data):
    par = abs(par)
    res = 0
    for (D, p_sim, r, h, _) in data:
        p_model = fines_model(D, r, h, par[0], par[1], par[2], par[3])
        res += np.sum((p_sim - p_model)**2)
    return res


def main():
    gradation_dir = pathlib.Path("~/DEMsim/size_distributions/kcv").expanduser()
    results_dir = pathlib.Path("~/DEMsim/results/kcv").expanduser()
    simulations = [
        Simulation(gradation=gradation_dir / "mono/kcv_2mm_0_1kg.dat",
                   results_directory=results_dir / "monosize_0_1kg_2mm/fines/1/compaction", color="b"),
        # Simulation(gradation=gradation_dir/"fuller/kcv_fines_0_1kg.dat",
        #            results_directory=results_dir/"uniaxial_0_1kg/fines/4/compaction", color="r"),
        # Simulation(gradation=gradation_dir/"fuller/kcv_fines_0_2kg.dat",
        #             results_directory=results_dir/"uniaxial_0_2kg/fines/1/compaction", color="b"),
        # Simulation(gradation=gradation_dir/"fuller/kcv_fines_0_2kg.dat",
        #           results_directory=results_dir/"uniaxial_0_2kg/fines/2/compaction", color="b"),
        # Simulation(gradation=gradation_dir/"fuller/kcv_fines_0_2kg.dat",
        #            results_directory=results_dir/"uniaxial_0_2kg/fines/5/compaction", color="b"),
        # Simulation(gradation=gradation_dir/"fuller/kcv_fines_0_5kg.dat",
        #            results_directory=results_dir/"uniaxial_0_5kg/fines/compaction", color="g"),
        # Simulation(gradation=gradation_dir/"fuller/kcv_fines_1_0kg.dat",
        #            results_directory=results_dir/"uniaxial_1_0kg/fines/compaction", color="k")

    ]
    data_sets = []
    for sim in simulations:
        particle_volume = get_particle_volume(sim.gradation)
        print(sim.results_directory)
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
        print(radius, height[0])
        data_sets.append((density[:idx], compaction_pressure[:idx], radius, height[:idx], sim.color))

    parameters_old = np.array([[0.5, 1e-3, 3166059e6, 4.2]])
    parameters_new = 0*parameters_old
    res = 10
    while res > 0:
        parameters_new = fmin(residual, parameters_old, args=(data_sets,), maxfun=1e4, maxiter=1e4)
        res = np.sum(abs(parameters_old - parameters_new))
        parameters_old = parameters_new
        print("\tsuggested parameters", parameters_new)
    parameters_new = abs(parameters_new)
    print(parameters_new)

    for (D_sim, _, r, h, c) in data_sets:
        p_model = fines_model(D_sim, r, h, *parameters_new)
        plt.plot(D_sim, p_model, '--' + c, lw=2)

    plt.xlabel("Relative density [-]")
    plt.ylabel("Pressure [MPa]")
    plt.tight_layout()
    plt.savefig("fines_all.png")
    plt.show()

if __name__ == '__main__':
    main()