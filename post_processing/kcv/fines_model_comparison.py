from collections import namedtuple
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

Case = namedtuple("Case", ["full_path", "continuum_path", "gradation", "color"])


def get_pressure_density_relationship(results_path, particle_volume):
    surface_position_file = results_path / "surface_positions.dou"
    surface_data = np.genfromtxt(surface_position_file, delimiter=",")
    radius = surface_data[0, 2]
    height = surface_data[:, -2]
    cylinder_volume = np.pi*radius**2*height
    print(radius, height[0])
    density = particle_volume/cylinder_volume
    surface_force_file = results_path / "surface_forces.dou"
    surface_forces = np.genfromtxt(surface_force_file, delimiter=",")
    compaction_pressure = surface_forces[:, -2]/np.pi/radius**2/1e6
    print(results_path, density[0])
    return np.column_stack((density, compaction_pressure))


def main():
    gradation_path = pathlib.Path("~/DEMsim/size_distributions/kcv").expanduser()
    results_path = pathlib.Path("~/DEMresults/kcv").expanduser()
    cases = [
        Case(full_path=results_path / "uniaxial_0_1kg/fines/compaction",
             continuum_path=results_path / "uniaxial_0_1kg/fines_continuum/compaction",
             gradation=gradation_path / "fuller/kcv_fines_0_1kg.dat",
             color="b"),

        Case(full_path=results_path/"uniaxial_0_1kg/full/compaction",
             continuum_path=results_path/"uniaxial_0_1kg/full_continuum/compaction",
             gradation=gradation_path/"fuller/kcv_full_0_1kg.dat",
             color="r")
    ]

    for case in cases:
        particle_volume = get_particle_volume(case.gradation)
        print(particle_volume)
        full_results = get_pressure_density_relationship(case.full_path, particle_volume)
        plt.plot(full_results[:, 0], full_results[:, 1], case.color, lw=2)

        continuum_results = get_pressure_density_relationship(case.continuum_path, particle_volume)
        plt.plot(continuum_results[:, 0], continuum_results[:, 1], "--" + case.color, lw=2)

    plt.show()


if __name__ == '__main__':
    main()
