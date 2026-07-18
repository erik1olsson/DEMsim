import pathlib

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

from pressure_density import get_particle_volume

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    distribution = "full"
    mass = 0.1

    mass_str = f"{str(mass).replace('.', '_')}kg"
    distribution_filename = f"kcv_{distribution}_" + mass_str + ".dat"
    distribution_file = pathlib.Path("~/DEMsim/size_distributions/kcv").expanduser() / distribution_filename
    particle_volume = get_particle_volume(distribution_file)
    print("Particle volume", particle_volume)
    compaction_path = pathlib.Path(f"~/DEMresults/kcv/uniaxial_{mass_str}/{distribution}/compaction").expanduser()
    surface_position_file = compaction_path / "surface_positions.dou"
    surface_data = np.genfromtxt(surface_position_file, delimiter=",")
    radius = surface_data[0, 2]
    height = surface_data[:, -2]
    cylinder_volume = np.pi*radius**2*height
    density = particle_volume / cylinder_volume
    surface_force_file = compaction_path / "surface_forces.dou"
    surface_forces = np.genfromtxt(surface_force_file, delimiter=",")
    compaction_pressure = surface_forces[:, -2]/np.pi/radius**2

    pa = compaction_pressure/1e6
    plt.plot(density, pa, lw=2)

    fines_path = pathlib.Path(f"~/DEMsim/results/kcv/uniaxial_{mass_str}//fines_continuum/compaction").expanduser()
    surface_position_file = fines_path / "surface_positions.dou"
    surface_force_file = fines_path / "surface_forces.dou"
    surface_data = np.genfromtxt(surface_position_file, delimiter=",")
    radius = surface_data[0, 2]
    height = surface_data[:, -2]
    cylinder_volume = np.pi*radius**2*height
    density = particle_volume/cylinder_volume
    surface_forces = np.genfromtxt(surface_force_file, delimiter=",")
    compaction_pressure = surface_forces[:, -2]/np.pi/radius**2
    pa = compaction_pressure/1e6
    plt.plot(density, pa, '--', lw=2)


    plt.xlabel("Density [-]")
    plt.ylabel("Pressure [MPa]")
    plt.xlim(0.45, 0.51)
    plt.ylim(0.01, 8)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
