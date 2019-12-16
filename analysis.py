import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat, unumpy
import os

class Output:
    def __init__(self, filename):
        data = np.loadtxt(filename)
        self.T   = data[:,0].flatten()
        self.E   = data[:,1]
        self.C   = data[:,2]
        self.chi = data[:,3]
        self.M   = data[:,4]

        offset = 4
        self.E_err   = data[:,1 + offset]
        self.C_err   = data[:,2 + offset]
        self.chi_err = data[:,3 + offset]
        self.M_err   = data[:,4 + offset]
        
L_values = [10, 16, 24, 36, 50]
#L_values = L_values[:3]

def get_filename(L):
    return "L{}output.txt".format(L)

def get_output(L):
    return Output(get_filename(L))

def Tc_from_chi(output):
    chi = output.chi
    idx = np.argmax(chi)
    return output.T[idx]

def extrapolate_Tc_from_chi(L_vals = L_values):
    Linv = 1 / np.array(L_vals, dtype=np.float)
    Tc = [Tc_from_chi(get_output(L)) for L in L_vals]

    _, intercept = np.polyfit(Linv, Tc, 1)
    return intercept

def chi_max(output):
    chi = output.chi
    chi_err = output.chi_err
    idx = np.argmax(chi)
    return ufloat(chi[idx], chi_err[idx])

def gamma_over_nu(Lvals = L_values):
    logL = np.log(Lvals)
    logchi = np.log([chi_max(get_output(L)).nominal_value for L in Lvals])
    slope, _ = np.polyfit(logL, logchi, 1)
    return slope

def get_M_at_Tc(L):
    output = get_output(L)
    return output.M[np.argmax(output.chi)]

def beta_over_nu(Lvals = L_values):
    logL = np.log(Lvals)
    logM = np.log([get_M_at_Tc(L) for L in Lvals])
    slope, _ = np.polyfit(logL, logM, 1)
    return -slope

def Tc_from_Cv(output):
    C = output.C
    idx = np.argmax(C)
    return output.T[idx]

def extrapolate_Tc_from_Cv(L_vals = L_values):
    Linv = 1 / np.array(L_vals, dtype=np.float)
    Tc = [Tc_from_Cv(get_output(L)) for L in L_vals]

    _, intercept = np.polyfit(Linv, Tc, 1)
    return intercept

def Cv_max(output):
    Cv = output.C
    Cv_err = output.C_err
    idx = np.argmax(Cv)
    return ufloat(Cv[idx], Cv_err[idx])

def entropy_at_idx(output, idx):
    ret = ufloat(0,0)
    T = output.T
    for i in range(1,idx + 1):
        C = ufloat(output.C[i], output.C_err[i])
        dT = T[i] - T[i - 1]
        ret += C / T[i] * dT
    return ret

def entropy(output):
    ret = []
    for i in range(len(output.T)):
        ret.append(entropy_at_idx(output, i).nominal_value)
    return np.array(ret)

def free_energy(L):
    output = get_output(L)
    S = entropy(output)
    return output.E - output.T * S * L**2

# Computes S(T = 4.5)
def entropy_limit(L):
    output = get_output(L)
    return entropy_at_idx(output, len(output.T) - 1)

if __name__ == '__main__':
    # Ensure output director exists
    output_dir = "results"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Part 4: Calculating T_c by extrapolating from susceptibility peaks
    Tc_chi = extrapolate_Tc_from_chi()
    print("T_c extrapolated from chi: {}".format(Tc_chi))
    with open(output_dir + "/Tc_from_chi.txt", "w") as f:
        print("{:.4f}".format(Tc_chi), file=f)

    # Part 5: Calculating \gamma / \nu from peaks in susceptibility
    gamma_nu = gamma_over_nu()
    print("gamma / nu from chi peaks: {}".format(gamma_nu))
    with open(output_dir + "/gamma_nu_from_chi.txt", "w") as f:
        print("{:.4f}".format(gamma_nu), file=f)

    # Part 6: Calcluating \beta / \nu from a linear fit
    beta_nu = beta_over_nu()
    print("beta / nu from M(T_c): {}".format(beta_nu))
    with open(output_dir + "/beta_nu_from_M_at_Tc.txt", "w") as f:
        print("{:.4f}".format(beta_nu), file=f)

    # Part 7: Calculating T_c from C_V
    Tc_Cv = extrapolate_Tc_from_Cv()
    print("T_c extrapolated from Cv: {}".format(Tc_Cv))
    with open(output_dir + "/Tc_from_Cv.txt", "w") as f:
        print("{:.4f}".format(Tc_Cv), file=f)

    # Part 8: Computing S(infty) for different values of L
    with open(output_dir + "/entropy_table.txt", "w") as f:
        print("\\begin{tabular}{c|", end="", file=f)
        for i in range(len(L_values)): print("c", end="", file=f)
        print("}", file=f)
        print("$L$", file=f)
        for L in L_values:
            print("& {}".format(L), file=f)
        print("\\\\", file=f)
        print("\\hline", file=f)
        print("$S(\infty)$", file=f)
        for L in L_values:
            print("S(L = {}) = {}".format(L, entropy_limit(L)))
            print("& {:S}".format(entropy_limit(L)), file=f)
        print("\\end{tabular}", file=f)


