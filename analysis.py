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

E_mean = np.mean(E, axis=1, keepdims=True)
N = E.shape[1]
E_errs = np.sqrt((N - 1.0) / N * np.sum((E - E_mean)**2, axis=1))
E_mean = E_mean.flatten()

plt.errorbar(T, E_mean, yerr=E_errs)
plt.show()
