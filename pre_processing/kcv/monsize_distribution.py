import pathlib

import numpy as np

import matplotlib.pyplot as plt
import matplotlib


matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{gensymb}')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    r = 2e-3
    mass = 0.5
    N = int(mass/2630*3/4/3.1415/r**3)
    print(N)

    radii = np.zeros(N) + r
    size_distribution_directory = pathlib.Path("../../size_distributions/kcv/mono")
    filename = "kcv_{r}mm_{mass}_kg".format(r=int(r*1000), mass=mass).replace(".", "_")
    filename += ".dat"
    np.savetxt(size_distribution_directory / filename, radii)


if __name__ == '__main__':
    main()
