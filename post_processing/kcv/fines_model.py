import pathlib
from zlib import DEF_MEM_LEVEL

import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib

from fines_model_comparison import get_pressure_density_relationship

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

# b = (D1 - D0r)/(Dmax - D0r)
#bDmax - bD0r = D1 - D0r
#D0r - bD0r = D1 - bDmax

class IncrementalModel:
    def __init__(self, initial_density):
        # Base material parameters
        self.D0 = 4.59801261e-01
        self.kl = 3.16605936e+06
        self.n = 4.20438353e+00

        self.a = 0.658433977547588
        # self.a = 0.8
        self.b = 0.9

        self.D = initial_density
        self.pressure = 0
        self.legend = 0

        # Internal state variable
        self.Dmax = initial_density

        # Help variables to make calculations easier
        self.pmax = 0

        self.D0u = self.D0
        self.ku = self.kl

        self.D0r = self.D0
        self.kr = self.kl


    def update(self, delta_d):
        self.D += delta_d

        if self.D > self.Dmax:
            # plastic loading
            self.Dmax = self.D

            if self.D > self.D0:
                self.pressure = self.kl*(self.D - self.D0)**self.n
                self.pmax = self.pressure
                self.D0u = self.D0 + self.a*(self.Dmax - self.D0)
                self.ku = self.pressure/(self.D - self.D0u)**self.n
                self.legend = 0

            else:
                self.pressure = 0
        else:
            if delta_d <= 0:
                # Unloading
                if self.D > self.D0u:
                    self.pressure = self.ku*(self.D - self.D0u)**self.n
                    gamma = (self.pressure/self.pmax)**(1/(self.b*self.n))
                    if gamma < 1:
                        self.D0r = (self.D - gamma*self.Dmax)/(1 - gamma)
                    else:
                        self.D0r = self.D0u
                    self.kr = self.pmax/(self.Dmax - self.D0r)**(self.b*self.n)
                    self.legend = 1
                else:
                    self.pressure = 0
            else:
                # Reloading
                if self.D > self.D0r:
                    self.pressure = self.kr*(self.D - self.D0r)**(self.b*self.n)
                    self.legend = 2



def model(D, D0, K, nl):
    p = K*(D - D0)**nl
    p[D < D0] = 0
    return p


def residual_loading(par, density, pressure):
    par = abs(par)
    if par[0] < np.min(density):
        par[0] = np.min(density)
    imax = np.argmax(density)
    p_model = model(density[:imax], par[0], par[1], par[2])
    return np.sum((p_model- pressure[:imax])**2)


def residual_unloading(par, density, pressure):
    par = abs(par)
    par[2] = 4.20438353e+00
    imax = np.argmax(density)
    p_model = model(density[imax:], par[0], par[1], par[2])
    return np.sum((p_model- pressure[imax:])**2)


def get_particle_volume(distribution_file):
    radii = np.genfromtxt(distribution_file, dtype=float)
    return 4*np.pi*np.sum(radii**3)/3


def main():
    distribution = "fines"
    particle_volume = get_particle_volume(size_distribution_file)
    data = get_pressure_density_relationship(distribution)
    density = data[:, 0]
    print(np.max(density))
    pa = data[:, 1]/1e6
    plt.plot(density, pa, lw=2)
    loading_par = fmin(residual_loading, [0.47, 3e5, 3], args=(density, pa), maxfun=1e6, maxiter=1e6)
    unloading_par = fmin(residual_unloading, [0.49, 3e5, 3], args=(density, pa), maxfun=1e6, maxiter=1e6)
    print(loading_par)
    print(unloading_par)
    print((unloading_par[0] - loading_par[0])/(np.max(density) - loading_par[0]))
    unloading_par[2] = 4.20438353e+00
    d = np.linspace(0.45, 0.51, 1000)
    p_model = model(d, loading_par[0], loading_par[1], loading_par[2])
    plt.plot(d, p_model, '--k', lw=2)

    p_model = model(d, unloading_par[0], unloading_par[1], unloading_par[2])
    plt.plot(d, p_model, '--k', lw=2)

    inc_model = IncrementalModel(0.45)
    D1 = 0.45
    D2 = 0.505
    d_levels = np.linspace(D1, D2, 6)
    current_d = D1
    unloading_d = 0.01
    d_vec = []
    pressure = []
    model_d = []
    legend = []
    for d_level in d_levels[1:]:
        loading = np.linspace(current_d, d_level, 1000)
        for d in loading:
            previous_d = current_d
            current_d = d
            d_vec.append(current_d)
            inc_model.update(current_d - previous_d)
            model_d.append(inc_model.D)
            pressure.append(inc_model.pressure)
            legend.append(inc_model.legend)

        unloading = np.linspace(current_d, d_level - unloading_d, 1000)
        for d in unloading:
            previous_d = current_d
            current_d = d
            d_vec.append(current_d)
            inc_model.update(current_d - previous_d)
            model_d.append(inc_model.D)
            pressure.append(inc_model.pressure)
            legend.append(inc_model.legend)

    plt.xlabel("Density [-]")
    plt.ylabel("Pressure [MPa]")
    plt.xlim(0.45, 0.51)
    plt.ylim(0.01, 8)
    plt.tight_layout()
    plt.figure(1)

    model_d = np.array(model_d)
    pressure = np.array(pressure)
    legend = np.array(legend)
    plt.plot(model_d, pressure, '--k')
    plt.plot(model_d[legend == 0], pressure[legend == 0], 'b*')
    plt.plot(model_d[legend == 1], pressure[legend == 1], 'g*')
    plt.plot(model_d[legend == 2], pressure[legend == 2], 'r*')

    plt.show()


if __name__ == '__main__':
    main()
