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


def main():
    r = 8
    H = 16

    a = 0.1
    p_vec = np.linspace(0, 1000, 1000)
    h = 0*p_vec
    z = np.linspace(0, H, 1000)
    for p0, c in zip([0.1, 0.5, 1, 2.5, 5.7], ["b", "r", "g", "k", "y"]):
        for i, p in enumerate(p_vec):


            szz = p*np.exp(-2*a/r*z) + p0

            n = 4.20438353e+00
            K = 3.16605936e+06
            D0 = 4.59801261e-01

            D_reploading = (szz/K)**(1/n) + D0
            plt.figure(2)
            D1 = (p0/K)**(1/n) + D0

            h[i] = np.trapz(D_reploading - D1, z)
        plt.plot(h, p_vec, '--' + c, lw=2)
    plt.xlim(0, 1)
    plt.ylim(0, 200)

    plt.show()


if __name__ == '__main__':
    main()
