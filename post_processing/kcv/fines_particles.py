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


log_filename = pathlib.Path("~/DEMsim/cmake-build-debug/log").expanduser()
data_lines = []
with open(log_filename, "r") as log_file:
    while line := log_file.readline():
        if line[0].isnumeric():
            data_lines.append(line)

with open("fines_data.dat", "w") as dat_file:
    for line in data_lines[:-1]:
        dat_file.write(line)


data = np.genfromtxt("fines_data.dat", delimiter=",")

p1_data = data[data[:, 0] == 1, :]
p2_data = data[data[:, 0] == 2, :]

data_set = p2_data

plt.figure(1)
plt.plot(data_set[:, 1], data_set[:, -13])

plt.figure(2)
plt.plot(data_set[:, 1])

plt.figure(3)
plt.plot(data_set[:, -13])

plt.xlabel("Time increment [-]")
plt.xlabel("Force [N]")
plt.savefig("particle_2_force.png")
plt.show()