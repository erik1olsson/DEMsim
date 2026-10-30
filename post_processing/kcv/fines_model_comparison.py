from collections import namedtuple
import pathlib

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from calibrate_fines_model import get_particle_volume

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

Case = namedtuple("Case", ["full_path", "continuum_path", "gradation", "color", "full_simulations",
                           "continuum_simulations"])


def get_pressure_density_relationship(results_path, particle_volume):
    surface_position_file = results_path / "surface_positions.dou"
    surface_data = np.genfromtxt(surface_position_file, delimiter=",")
    radius = surface_data[0, 2]
    height = surface_data[:, -2]
    cylinder_volume = np.pi*radius**2*height
    density = particle_volume/cylinder_volume
    surface_force_file = results_path / "surface_forces.dou"
    surface_forces = np.genfromtxt(surface_force_file, delimiter=",")
    compaction_pressure = surface_forces[:, -2]/np.pi/radius**2/1e6
    length = np.min((len(density), len(compaction_pressure)))
    return np.column_stack((density[:length], compaction_pressure[:length]))


def main():
    gradation_path = pathlib.Path("~/DEMsim/size_distributions/kcv/mono").expanduser()
    results_path = pathlib.Path("~/DEMresults/kcv/monosize").expanduser()
    cases = [
        # Case(full_path=results_path / "uniaxial_0_1kg/fines/compaction",
        #      continuum_path=results_path / "uniaxial_0_1kg/fines_continuum/compaction",
        #      gradation=gradation_path / "fuller/kcv_fines_0_1kg.dat",
        #      color="b"),

        Case(full_path=results_path / "0_1kg_2mm/fines/",
             continuum_path=results_path / "0_1kg_2mm/fines_continuum/",
             gradation=gradation_path / "kcv_2mm_0_1kg.dat",
             color="b", full_simulations=[1,2,3], continuum_simulations=[1]),
        Case(full_path=results_path/"0_2kg_2mm/fines/",
             continuum_path=results_path/"0_2kg_2mm/fines_continuum/",
             gradation=gradation_path/"kcv_2mm_0_2kg.dat",
             color="r", full_simulations=[1, 2, 3], continuum_simulations=[1]),
        Case(full_path=results_path/"0_5kg_2mm/fines/",
             continuum_path=results_path/"0_5kg_2mm/fines_continuum/",
             gradation=gradation_path/"kcv_2mm_0_5kg.dat",
             color="g", full_simulations=[1, 2, 3], continuum_simulations=[1]),
        Case(full_path=results_path/"1_0kg_2mm/fines/",
             continuum_path=results_path/"1_0kg_2mm/fines_continuum/",
             gradation=gradation_path/"kcv_2mm_1_0kg.dat",
             color="k", full_simulations=[1, 2, 3], continuum_simulations=[1])

    ]

    for case in cases:
        particle_volume = get_particle_volume(case.gradation)
        print(particle_volume)
        for full_sim in case.full_simulations:
            directory = case.full_path / str(full_sim) / "compaction"
            full_results = get_pressure_density_relationship(directory, particle_volume)
            plt.plot(full_results[:, 0], full_results[:, 1], case.color, lw=2)

        for continuum_sim in case.continuum_simulations:
            directory = case.continuum_path / str(continuum_sim) / "compaction"
            continuum_results = get_pressure_density_relationship(directory, particle_volume)
            plt.plot(continuum_results[:, 0], continuum_results[:, 1], "--" + case.color, lw=2)

    # plt.xlim(0.44, 0.48)
    # plt.ylim(0, 1)
    plt.xlabel("Relative density [-]")
    plt.ylabel("Pressure [MPa]")
    plt.tight_layout()
    plt.savefig("fines_continuum_zoom.png")
    plt.show()


if __name__ == '__main__':
    main()
