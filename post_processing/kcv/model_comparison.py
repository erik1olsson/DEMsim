import pathlib

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from fines_model import get_particle_volume

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def pressure_volume_relationship(directory, particle_volume):
    surface_position_file = directory / "surface_positions.dou"
    surface_data = np.genfromtxt(surface_position_file, delimiter=",")
    radius = surface_data[0, 2]
    height = surface_data[:, -2]
    cylinder_volume = np.pi*radius**2*height
    print(cylinder_volume)
    relative_density = particle_volume/cylinder_volume
    surface_force_file = directory / "surface_forces.dou"
    surface_forces = np.genfromtxt(surface_force_file, delimiter=",")
    compaction_pressure = surface_forces[:, -2]/np.pi/radius**2
    return np.column_stack((relative_density, compaction_pressure/1e6))



def main():
    for mass, lw in zip([0.1, 1.], [2, 3]):
        for distribution, c in zip(["full", "fines"], ["b", "r"]):
            mass_str = f"{str(mass).replace('.', '_')}kg"
            distribution_filename = f"kcv_{distribution}_" + mass_str + ".dat"
            distribution_file = pathlib.Path("~/DEMsim/size_distributions/kcv").expanduser()/distribution_filename
            particle_volume = get_particle_volume(distribution_file)
            print("Particle volume", particle_volume)

            results_path = pathlib.Path(f"~/DEMresults/kcv/uniaxial_{mass_str}/{distribution}/compaction").expanduser()
            sim_data = pressure_volume_relationship(results_path, particle_volume)
            plt.plot(sim_data[:, 0], sim_data[:, 1], c, lw=lw)

            # fines_path = pathlib.Path(f"~/DEMsim/results/kcv/uniaxial_{mass_str}/fines_continuum/compaction").expanduser()
            # sim_data = pressure_volume_relationship(fines_path, particle_volume)
            # plt.plot(sim_data[:, 0], sim_data[:, 1], '--b', lw=3)

    plt.xlabel("Density [-]")
    plt.ylabel("Pressure [MPa]")
    plt.ylim(0.01, 8)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
