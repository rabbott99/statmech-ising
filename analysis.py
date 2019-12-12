import numpy as np
import matplotlib.pyplot as plt

def read_output(filename):
    data = np.loadtxt(filename)
    T = data[:,0].flatten()
    # energy
    E = data[:,1::3]
    # specific heat
    C = data[:,2::3]
    # susceptibility
    chi = data[:,3::3]

    return T, E, C, chi

T, E, C, chi = read_output("out.txt")

E_mean = E[:,0]
E_errs = E[:,1]

plt.errorbar(T, E_mean, yerr=E_errs)
plt.show()
