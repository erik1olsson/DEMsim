import pathlib

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


def main():
    DEMpath = pathlib.Path('~/DEMresults/kcv/monosize/0_1kg_2mm/particle_response/1/compaction').expanduser()
    force_data = np.genfromtxt(DEMpath / "surface_forces.dou", delimiter=',')
    surface_data = np.genfromtxt(DEMpath/"surface_positions.dou", delimiter=',')

    force = force_data[:, -2]
    displacement = surface_data[:, -2]
    # displacement -= displacement[0]

    plt.plot(-displacement, force, 'r', lw=2)

    DEMpath = pathlib.Path('~/DEMresults/kcv/monosize/0_1kg_2mm/particle_response_continuum/particle_compaction').expanduser()
    force_data = np.genfromtxt(DEMpath/"surface_forces.dou", delimiter=',')
    surface_data = np.genfromtxt(DEMpath/"surface_positions.dou", delimiter=',')

    force = force_data[:, -2]
    displacement = surface_data[:, -2]
    # displacement -= displacement[0]

    plt.plot(-displacement, force, '--r', lw=2)

    plt.show()


if __name__ == '__main__':
    main()
