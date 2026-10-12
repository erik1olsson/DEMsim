import math
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
    vf_fines = 0.25
    mass = 1
    r_fines  = 2e-3
    r_stones = 8e-3
    N_fines = int(vf_fines*3*mass/2630/4/3.1415/r_fines**3)
    N_stones = int((1-vf_fines)*3*mass/2630/4/3.1415/r_stones**3)

    print(N_fines, N_stones)

    radii = np.zeros(N_stones + N_fines)
    radii[:N_stones] = r_stones
    radii[N_stones:] = r_fines

    filename = "kcv_{r1}mm-{r2}mm_vf={vf}_mass={mass}kg".format(r1=int(r_fines*1000), r2=int(r_stones*1000), vf=vf_fines,
                                                          mass=mass).replace(".", "_")

    filename += ".dat"
    size_distribution_directory = pathlib.Path("../../size_distributions/kcv/bimodal")
    np.savetxt(size_distribution_directory / filename, radii)



if __name__ == '__main__':
    main()
