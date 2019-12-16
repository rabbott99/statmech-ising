from analysis import *
from uncertainties import umath

# Make matplotlib work with latex taken from here:
# https://matplotlib.org/3.1.1/tutorials/text/usetex.html
from matplotlib import rc
rc('text', usetex=True)

def plot_feature(L, feature, normalize = False):
    filename = "L{}output.txt".format(L)
    output = Output(filename)
    T = output.T
    vals = getattr(output, feature)
    errs = getattr(output, feature + "_err")

    if normalize:
        vals = vals / L**2
        errs = errs / L**2

    plt.errorbar(T, vals, yerr=errs, label="L = {}".format(L))

# Plotting energy
for L in L_values:
    plot_feature(L, "E", True)
plt.title("Energy vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Energy per spin")
plt.legend()
plt.savefig("energy-plot.png")
plt.cla()

# Plotting Specific heat
for L in L_values:
    plot_feature(L, "C", False)
plt.title("Specific Heat vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Specific Heat")
plt.legend()
plt.savefig("specific-heat-plot.png")
plt.cla()

# Plotting manetization
for L in L_values:
    plot_feature(L, "M", False)
plt.title("Magnetization vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Average Magnetization")
plt.legend()
plt.savefig("magnetization-plot.png")
plt.cla()


# Plotting susceptibility
for L in L_values:
    plot_feature(L, "chi", False)
plt.title("Susceptibility vs Temperature")
plt.xlabel("Temperature")
plt.ylabel("Susceptibility")
plt.legend()
plt.savefig("susceptibility-plot.png")
plt.cla()

# Plotting chi(Tc) vs L
x = np.log(np.array(L_values))
chi_at_Tc = [umath.log(chi_max(get_output(L))) for L in L_values]
y = [x.nominal_value for x in chi_at_Tc]
yerr = [x.std_dev for x in chi_at_Tc]
plt.errorbar(x, y, yerr=yerr, fmt='o')

# Plot linear fit
slope, intercept = np.polyfit(x, y, 1)
y = slope * x + intercept
plt.plot(x, y, color='red')

plt.title("Susceptibility as a function of $\log L$")
plt.xlabel("$\log L$")
plt.ylabel("$\\log \\chi(T_c(L))$")
plt.savefig("chi-max-plot.png")
plt.cla()

# Plotting Tc(L) vs L, where Tc is taken from the susceptibility
x = 1 / np.array(L_values)
Tc = [Tc_from_chi(get_output(L)) for L in L_values]
plt.scatter(x, Tc)

# Plot linear fit
slope, intercept = np.polyfit(x, Tc, 1)
x = np.arange(0, 0.12, 0.01)
y = slope * x + intercept
plt.plot(x, y, color='red')

plt.title("$T_c$ from susceptibility vs $L^{-1}$")
plt.xlabel("$L^{-1}$")
plt.ylabel("$T_c(L)$")
plt.savefig("chi-Tc-vs-Linv-plot.png")
plt.cla()

# Plotting Tc(L) vs L, where Tc is taken from the specific heat
x = 1 / np.array(L_values)
Tc = [Tc_from_Cv(get_output(L)) for L in L_values]
plt.scatter(x, Tc)

# Plot linear fit
slope, intercept = np.polyfit(x, Tc, 1)
x = np.arange(0, 0.12, 0.01)
y = slope * x + intercept
plt.plot(x, y, color='red')

plt.title("$T_c$ from specific heat vs $L^{-1}$")
plt.xlabel("$L^{-1}$")
plt.ylabel("$T_c(L)$")
plt.savefig("Cv-Tc-vs-Linv-plot.png")
plt.cla()

# Plotting L^{-\gamma / \nu} \chi as a function of L^{1 / \nu} (T - T_c(L))
def plot_chi_scaling(L, gamma_nu, nu):
    output = get_output(L)
    T = output.T

    Tc = Tc_from_chi(output)
    scaling_var = L**(1 / nu) * (T - Tc)
    
    chi = output.chi * L**(-gamma_nu)
    chi_err = output.chi_err * L**(-gamma_nu)
    
    plt.errorbar(scaling_var, chi, yerr=chi_err, label="L = {}".format(L))
    
# Using exact gamma / nu = 1.75
for L in L_values:
    plot_chi_scaling(L, 1.75, 1)
plt.title("$L^{-\\gamma / \\nu} \chi$ as a function of $L^{1 / \\nu} (T - T_c(L))$ ($\\gamma / \\nu$ exact)")
plt.xlabel("$L^{1 / \\nu} (T - T_c(L))$")
plt.ylabel("$L^{-\\gamma / \\nu} \chi$")
plt.xlim([-25, 25])
plt.legend()
plt.savefig("chi-scaling-gamma-nu-exact.png")
plt.cla()

# Using gamma / nu fit from fitting log(chi(T_c)) vs. log(L)
gamma_nu = gamma_over_nu(L_values)
for L in L_values:
    plot_chi_scaling(L, gamma_nu, 1)
plt.title("$L^{-\\gamma / \\nu} \chi$ as a function of $L^{1 / \\nu} (T - T_c(L))$ ($\\gamma / \\nu$ from fit)")
plt.xlabel("$L^{1 / \\nu} (T - T_c(L))$")
plt.ylabel("$L^{-\\gamma / \\nu} \chi$")
plt.legend()
plt.xlim([-25, 25])
plt.savefig("chi-scaling-gamma-nu-fit.png")
plt.cla()

# Plotting L^{-\beta / \nu} M as a function of L^{1 / \nu} (T - T_c(L))
def plot_M_scaling(L, beta_nu, nu):
    output = get_output(L)
    T = output.T

    Tc = Tc_from_chi(output)
    scaling_var = L**(1 / nu) * (T - Tc)
    
    M = output.M * L**(beta_nu)
    M_err = output.M_err * L**(beta_nu)
    
    plt.errorbar(scaling_var, M, yerr=M_err, label="L = {}".format(L))

beta_nu = beta_over_nu(L_values)
for L in L_values:
    plot_M_scaling(L, beta_nu, 1)
plt.title("$L^{-\\beta / \\nu} M$ as a function of $L^{1 / \\nu} (T - T_c(L))$ ($\\beta / \\nu$ from fit)")
plt.xlabel("$L^{1 / \\nu} (T - T_c(L))$")
plt.ylabel("$L^{-\\beta / \\nu} M$")
plt.legend()
plt.savefig("M-scaling.png")
plt.cla()

# Plotting C_V(T_c) versus \log L
x = np.log(L_values)
y = [Cv_max(get_output(L)).nominal_value for L in L_values]
plt.scatter(x, y)
plt.title("$C_V$ as a function of $\log L$")
plt.xlabel("$\log L$")
plt.ylabel("$C_V(T_c(L))$")
plt.savefig("Cv-max-plot.png")
plt.cla()

# Plotting entropy vs temperature
for L in L_values:
    output = get_output(L)
    x = output.T
    y = entropy(output)
    plt.plot(x, y, label="L = {}".format(L))

plt.title("Entropy as a function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Entropy")
plt.legend()
plt.savefig("entropy-plot.png")
plt.cla()

# Plotting free energy vs temperature
for L in L_values:
    output = get_output(L)
    x = output.T
    y = free_energy(L) / L**2
    plt.plot(x, y, label="L = {}".format(L))

plt.title("Free energy per spin as a function of temperature")
plt.xlabel("Temperature")
plt.ylabel("Free energy per spin")
plt.legend()
plt.savefig("free-energy-plot.png")
