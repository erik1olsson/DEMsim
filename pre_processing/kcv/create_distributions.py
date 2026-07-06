import os
import pathlib

import numpy as np
import sympy as sp

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
    density = 2650
    target_mass = 1.

    rmax_val = 8
    rmin_val = 0.5

    cutoff_radius = 2.5

    r, x, rmax, rmin, rc, B = sp.symbols(['r', 'x', 'rmax', 'rmin', 'rc', 'B'], positive=True)
    GV_fuller = (r/rmax)**0.5

    GV_fines = GV_fuller.subs(r, rc)
    GV_stones = (GV_fuller - GV_fines)/(1-GV_fines)


    GN_stones = (GV_stones/r**3 + 3*sp.integrate(GV_stones.subs(r, x)/x**4, [x, rc, r]))
    A = 1/GN_stones.subs(r, rmax)
    GN_stones *= A

    f_fuller = sp.lambdify((r, rmax), GV_fuller, modules="numpy")
    f_stones = sp.lambdify((r, rmax, rc), GV_stones, modules="numpy")
    fn_stones = sp.lambdify((r, rmax, rc), GN_stones, modules="numpy")
    f_fines = sp.lambdify((rmax, rc), GV_fines, modules="numpy")


    rc_val = 0.5

    r_vals_fuller = np.linspace(rmin_val, rmax_val, 2000)
    r_vals_stones = np.linspace(rc_val, rmax_val, 2000)
    vf_fines = f_fines(rmax_val, rc_val)*0
    y_fuller = f_fuller(r_vals_fuller, rmax_val)
    y_stones = f_stones(r_vals_stones, rmax_val, rc_val)
    yn_stones = fn_stones(r_vals_stones, rmax_val, rc_val)

    current_stone_mass = 0.
    radii = []
    current_total_mass = 0
    while current_total_mass < target_mass:
        val = np.random.random()
        radius = np.interp(val, yn_stones, r_vals_stones)
        current_stone_mass += 4*3.1415*(radius/1000)**3/3*density
        radii.append(radius)
        current_total_mass = current_stone_mass/(1-vf_fines)
    print(len(radii))
    radii = np.array(sorted(radii))
    plt.plot(r_vals_fuller, y_fuller, 'k',  lw=2, label="GV_fuller")
    plt.plot(r_vals_stones, y_stones, 'b', lw=2, label="GV_stones")
    plt.plot(r_vals_stones, yn_stones, 'r', lw=2, label="GN_stones")

    plt.plot(radii, np.arange(len(radii))/len(radii), 'r*')
    plt.plot(radii, np.cumsum(4*3.1415*(radii/1000)**3/3*density)/current_stone_mass, 'b*')

    plt.xlim(rmin_val, rmax_val)

    size_distribution_directory = pathlib.Path("../../size_distributions/kcv")
    if not size_distribution_directory.is_dir():
        os.makedirs(size_distribution_directory)
    filename = "kcv_full_{mass}kg.dat".format(mass=str(target_mass).replace(".", "_"))
    np.savetxt(size_distribution_directory / filename, radii/1e3)

    large_stones = radii[radii > cutoff_radius]
    fines = radii[radii <= cutoff_radius]
    print(large_stones.shape)

    filename = "kcv_large_{mass}kg.dat".format(mass=str(target_mass).replace(".", "_"))
    np.savetxt(size_distribution_directory / filename, large_stones/1e3)

    filename = "kcv_fines_{mass}kg.dat".format(mass=str(target_mass).replace(".", "_"))
    np.savetxt(size_distribution_directory / filename, fines/1e3)

    plt.legend(loc='best')
    plt.grid(True)
    plt.show()


if __name__ == '__main__':
    main()
