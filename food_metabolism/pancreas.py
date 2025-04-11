# https://www.researchgate.net/figure/Concentration-response-of-insulin-secretion-A-Average-insulin-secretion-rate_fig3_336211327
# https://www.sciencedirect.com/science/article/pii/S0168822713004221

def pancreas_dummy_ode(t, y):
    kcl = 0.15

    Gblood = y[0]
    dGblood = -kcl * y[0]
    dinsulin, dglucagon, dsomat = pancreatic_secretion_response_to_Gblood(Gblood)
    dGcl = +kcl * y[0]
    dydt = [dGblood, dinsulin, dglucagon, dsomat, dGcl]
    return dydt

def pancreatic_secretion_response_to_Gblood(Gblood):
    """
        provides piecewise functions approximating insulin, glucagon, and
        somatostatin secretion response relative to plasma glucose concentrations (mM).
        Note that insulin and somatostatin secretions are dependent on a
        baseline maximum secretion rate, while glucagon is based on the baseline
        secretion in the absence of glucose glucagon is based on the baseline secretion
        in the absence of glucose.
    """
    insulin_m1 = (.75 - .0) / (11 - 2.5)
    insulin_m2 = (.85 - .75) / (17 - 11)
    insulin_m3 = (1 - 0.85) / (30 - 17)
    insulin_x1 = 2.5
    insulin_x2 = -34
    insulin_x3 = -170/3

    glucagon_m1 = (0.3 - 1) / (7.5 - 0)
    glucagon_m2 = (0.55 - 0.3) / (20 - 7.5)
    glucagon_m3 = (1.3 - 0.55) / (30 - 20)
    __glucagon_m2 = (0 - 0.3) / (30 - 7.5)
    glucagon_x1 = 75/7
    glucagon_x2 = -7.5
    __glucagon_x2 = 30
    glucagon_x3 = 38/3

    somatostatin_m1 =  1 / 30
    # Gblood is mmol
    k_baseline_insulin_secretion_max = 1
    k_baseline_glucagon_secretion_empty = 1
    k_baseline_somatostatin_secretion_max = 1

    # INSULIN
    dinsulin = 0
    if Gblood >= 30:
        dinsulin +=  k_baseline_insulin_secretion_max
    elif Gblood >= 17:
        dinsulin += insulin_m3 * k_baseline_insulin_secretion_max * (Gblood - insulin_x3)
    elif Gblood >= 11:
        dinsulin += insulin_m2 * k_baseline_insulin_secretion_max * (Gblood - insulin_x2)
    elif Gblood >= 2.5:
        dinsulin += insulin_m1 * k_baseline_insulin_secretion_max * (Gblood - insulin_x1)


    # GLUCAGON
    dglucagon = 0
    # if Gblood >= 30:
    #     dglucagon += 1.3 * k_baseline_glucagon_secretion_empty
    # elif Gblood >= 20:
    #     dglucagon += glucagon_m3 * k_baseline_glucagon_secretion_empty * (Gblood - glucagon_x3)
    # elif Gblood >= 7.5:
    #     dglucagon += __glucagon_m2 * k_baseline_glucagon_secretion_empty * (Gblood - __glucagon_x2)
    # elif Gblood >= 0:
    #     dglucagon += glucagon_m1 * k_baseline_glucagon_secretion_empty * (Gblood - glucagon_x1)
     
    if Gblood <= 7.5:
        dglucagon += glucagon_m1 * k_baseline_glucagon_secretion_empty * (Gblood - glucagon_x1)
    elif Gblood <= 30:
        dglucagon += __glucagon_m2 * k_baseline_glucagon_secretion_empty * (Gblood - __glucagon_x2)

    # SOMATOSTATIN
    dsomatostatin = somatostatin_m1 * k_baseline_somatostatin_secretion_max * Gblood

    return (dinsulin, dglucagon, dsomatostatin)

def plot_hormone_secretion_from_G_conc():
    import matplotlib.pyplot as plt
    import numpy as np

    Gblood = np.linspace(0, 30, 31)
    hormones = np.array([pancreatic_secretion_response_to_Gblood(glucose) for glucose in Gblood])

    hormone_labels = ('Insulin', 'Glucagon', 'Somatostatin')
    colors = ('blue', 'green', 'red')

    plt.figure(figsize=(10, 6))
    for i in range(hormones.shape[1]):
       plt.plot(Gblood, hormones[:, i], label=hormone_labels[i], color=colors[i])
    plt.grid(color='gray', linestyle='--', linewidth=0.5)

    plt.xlabel('Blood Glucose Concentration (mM)')
    plt.ylabel('Hormone Secretion (%)')
    plt.title('Hormone Secretion in Response to Blood Glucose Levels')
    plt.legend()
    plt.show()

def plot_ode_soln():
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.integrate import solve_ivp

    t_span = (0, 6)
    y0_Gblood = 1
    y0 = [y0_Gblood, 0, 0, 0, 0]
    sol = solve_ivp(pancreas_dummy_ode, t_span=t_span, y0=y0)
    names = ("Plasma Glucose", "Insulin", "Glucagon", "Somatostatin", "Cleared Glucose")
    # mass_balance = 1 - (sol.y[0, :] + sol.y[4, :])
    # plt.plot(mass_balance)
    # plt.show()
    # plt.figure(figsize=(10, 6))
    # for i in range(len(y0)):
    #     plt.plot(sol.t, sol.y[i, :], label=names[i])
    # plt.legend()
    # plt.show()

    fig, (ax_main, ax_mass) = plt.subplots(
        nrows=2, 
        figsize=(10, 6), 
        gridspec_kw={'height_ratios': [4, 1]}, 
        sharex=True
    )

    for i in range(len(y0)):
        ax_main.plot(sol.t, sol.y[i, :], label=names[i])
    ax_main.legend()
    ax_main.set_ylabel("Concentration (Temporary Placeholder)")
    ax_main.set_title("Pancreatic Hormone Secretion Kinetics")

    mass_balance = 1 - (sol.y[0, :] + sol.y[4, :])
    ax_mass.plot(sol.t, mass_balance, color='gray')
    ax_mass.set_ylabel("Mass\nBalance")
    ax_mass.set_xlabel("Time")
    ax_mass.set_ylim([-1e-12, 1e-12])
    ax_mass.axhline(0, linestyle='--', color='black', linewidth=0.8)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # plot_hormone_secretion_from_G_conc()
    # plot_ode_soln()
