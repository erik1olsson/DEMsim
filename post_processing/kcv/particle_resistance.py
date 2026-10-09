import pathlib
from collections import namedtuple

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

Simulation = namedtuple("Simulation", ["radius", "distance", "levels", "material"])

fines_parameters = {
    "mono_6000": (0.45660371680201517, 599263.7337328821, 3.1993378400249726),
    "uniaxial_0_1kg": (0.45660371680201517, 599263.7337328821, 3.1993378400249726)
}

def main():
    cases = [
        Simulation(radius=8, distance=24, levels=[0.1, 0.5, 1, 2.5, 5.7], material="uniaxial_0_1kg"),
        Simulation(radius=8, distance=24, levels=[0.1, 0.5, 1, 2.5, 5, 7.5], material="mono_6000"),
        Simulation(radius=4, distance=24, levels=[0.1, 0.5, 1, 2.5, 5.7], material="uniaxial_0_1kg"),
        Simulation(radius=8, distance=10, levels=[0.1, 0.5, 1, 2.5, 5.7], material="uniaxial_0_1kg"),
        # Simulation(radius=8, distance=16, levels=[0.1, 0.5, 1, 2.5, 5.7], material="uniaxial_0_1kg")
    ]
    colors = ["b", "r", "g", "k", "y", "m"]

    for fig, case in enumerate(cases):
        base_directory = pathlib.Path("~/DEMresults/kcv/" + case.material).expanduser()
        sim_directory = "particle_resistance_{r}_{d}".format(r=case.radius, d=case.distance)
        # particle_data = np.genfromtxt(base_directory / sim_directory / "animation" / "particles" /
        #                              "particles_0.dou", delimiter=",")
        # particle_radii = particle_data[: ,7]
        # particle_volume = np.sum(4*3.1415*particle_radii**3/3)

        for p, c in zip(case.levels, colors):
            data = np.genfromtxt(base_directory / sim_directory / "particle_resistance_p={p}MPa".format(p=p) /
                                 "particle_3.dou",  delimiter=",")
            wall_data = np.genfromtxt(base_directory / sim_directory / "particle_resistance_p={p}MPa".format(p=p) /
                                 "surface_positions.dou",  delimiter=",")
            r = wall_data[0, 2]
            h = wall_data[0, -2]
            print(data[0, 4])
            z = data[0, 4] - data[:, 4]
            fz = data[:, -1]
            idx = np.argmax(z)
            z_loading = z[:idx]
            fz_loading = fz[:idx]
            pz_loading = fz_loading/case.radius**2/3.1415

            z_unloading = z[idx:]
            fz_unloading = fz[idx:]
            pz_unloading = fz_unloading/case.radius**2/3.1415
            plt.figure(10+fig)
            plt.plot(z_loading, fz_loading, c, lw=2)



            n = fines_parameters[case.material][2]
            K = fines_parameters[case.material][1]
            D0 = fines_parameters[case.material][0]
            cyL_vol = 3.1415*r**2*h
            # D_true = particle_volume/cyL_vol
            # D_fines = (particle_volume - 4*3.1415*(case.radius/1000)**3/3)/(cyL_vol - 4*3.1415*(case.radius/1000)**3/3)
            D1 = (p/K)**(1/n) + D0
            D_reploading = ((pz_loading)/K)**(1/n)

            # D_model = 0.06*(1 - np.exp(-1500*z_loading)) + D1
            D_model = 0.03*(1-np.exp(-1500*z_loading)) + D1
            # D_model = 10*z_loading
            p_model = K*(D_model - D0)**n
            print(D_model)
            fz_model = (p_model - p)*case.radius**2*3.1415
            plt.plot(z_loading, fz_model, '--' + c, lw=2)

            plt.figure(1)
            plt.plot(z_loading, D_reploading + D0 - D1, c, lw=2)

            D0u = D0 + 0.8*(D_reploading[-1] - D0)
            Ku = pz_loading[-1]/(D_reploading[-1] - D0u)**n
            pz_unloading[pz_unloading < 0] = 0
            Drep_unloading = (pz_unloading/Ku)**(1/n) + D0u
            plt.figure(1)
            # plt.plot(z_unloading, Drep_unloading - D1, c, lw=2)

    x = np.linspace(0, 0.008, 1000)
    y = 0.03*(1-np.exp(-1500*x))
    # y = 30*x
    plt.figure(1)
    plt.plot(x, y, '--k', lw=3)
    plt.show()


if __name__ == '__main__':
    main()
