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

Simulation = namedtuple("Simulation", ["radius", "distance", "levels"])


def model(d, Kp, n, radius, pressure):
    return Kp[pressure]*d**n*radius**2


def model_residual(par, data):
    Kp = {val: p for (val, p) in zip(par[0:3], [0.1, 0.5, 1])}
    for data_set in simulations:
        pass


def main():
    cases = [
        Simulation(radius=8, distance=24, levels=[0.1, 0.5, 1]),
        Simulation(radius=4, distance=24, levels=[0.1])
    ]
    colors = ["b", "r", "g", "k"]
    datasets = {}
    for case in cases:
        base_directory = pathlib.Path("~/DEMresults/kcv/uniaxial_0_1kg").expanduser()
        sim_directory = "particle_resistance_{r}_{d}".format(r=case.radius, d=case.distance)
        for p, c in zip(case.levels, colors):
            data = np.genfromtxt(base_directory / sim_directory / "particle_resistance_p={p}MPa".format(p=p) /
                                 "particle_3.dou",  delimiter=",")
            z = data[0, 4] - data[:, 4]
            fz = data[:, -1]


            plt.plot(z, fz/case.radius**2, c, lw=2)
    plt.show()

if __name__ == '__main__':
    main()
