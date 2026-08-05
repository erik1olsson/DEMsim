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
    N = 6000
    r = 1
    radii = np.zeros(N) + r
    size_distribution_directory = pathlib.Path("../../size_distributions/kcv")
    filename = "kcv_{r}mm_{N}".format(r=r, N=N).replace(".", "_")
    filename += ".dat"
    np.savetxt(size_distribution_directory / filename, radii/1e3)


if __name__ == '__main__':
    main()
