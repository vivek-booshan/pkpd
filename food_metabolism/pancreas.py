# https://www.researchgate.net/figure/Concentration-response-of-insulin-secretion-A-Average-insulin-secretion-rate_fig3_336211327
# https://www.sciencedirect.com/science/article/pii/S0168822713004221

def pancreas(Gblood):
    insulin_m1 = (.75 - .0) / (11 - 2.5)
    insulin_m2 = (.85 - .75) / (17 - 11)
    insulin_m3 = (1 - 0.85) / (30 - 17)
    insulin_x1 = 2.5
    insulin_x2 = -34
    insulin_x3 = -170/3

    glucagon_m1 = (0.3 - 1) / (7.5 - 0)
    glucagon_m2 = (0.55 - 0.3) / (20 - 7.5)
    glucagon_m3 = (1.3 - 0.55) / (30 - 20)
    glucagon_x1 = 75/7
    glucagon_x2 = -7.5
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
    if Gblood >= 30:
        dglucagon += 1.3 * k_baseline_glucagon_secretion_empty
    elif Gblood >= 20:
        dglucagon += glucagon_m3 * k_baseline_glucagon_secretion_empty * (Gblood - glucagon_x3)
    elif Gblood >= 7.5:
        dglucagon += glucagon_m2 * k_baseline_glucagon_secretion_empty * (Gblood - glucagon_x2)
    elif Gblood >= 0:
        dglucagon += glucagon_m1 * k_baseline_glucagon_secretion_empty * (Gblood - glucagon_x1)

    # SOMATOSTATIN
    dsomatostatin = somatostatin_m1 * k_baseline_somatostatin_secretion_max * Gblood

    return (dinsulin, dglucagon, dsomatostatin)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    Gblood = np.linspace(0, 30, 31)
    hormones = np.array([pancreas(glucose) for glucose in Gblood])

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

