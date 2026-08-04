import pathlib

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

def main():
    directory = pathlib.Path("~/DEMresults/kcv/uniaxial_0_1kg/particle_restance_8_24").expanduser()
    pressure_levels = [0.1, 0.5]
    for p in pressure_levels:
        data = np.genfromtxt(directory / "particle_resistance_p={p}MPa".format(p=p) / "particle_3.dou",
                             delimiter=",")
        z = data[0, 4] - data[:, 4]
        fz = data[:, -1]

        plt.plot(z, fz, lw=2)
    plt.show()

if __name__ == '__main__':
    main()
