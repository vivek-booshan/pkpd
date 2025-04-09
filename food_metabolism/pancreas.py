# https://www.researchgate.net/figure/Concentration-response-of-insulin-secretion-A-Average-insulin-secretion-rate_fig3_336211327
im1 = (.75 - .0) / (11 - 2.5)
im2 = (.85 - .75) / (17 - 11)
im3 = (1 - 0.85) / (30 - 17)
ix1 = 2.5
ix2 = -34
ix3 = -170/3

gm1 = (0.3 - 1) / (7.5 - 0)
gm2 = (0.55 - 0.3) / (20 - 7.5)
gm3 = (1.3 - 0.55) / (30 - 20)
gx1 = 75/7
gx2 = 7.5
gx3 = 38/3

sm1 =  1 / 30

def pancreas(Gblood):
# Gblood is mmol

    # INSULIN
    dinsulin = 0
    if Gblood > 30:
        dinsulin += k_baseline_max
    elif Gblood > 17:
        dinsulin += im3 * k_baseline_max * (Gblood - ix3)
    elif Gblood > 11:
        dinsulin += im2 * k_baseline_max * (Gblood - ix2)
    elif Gblood > 2.5:
        dinsulin += im1 * k_baseline_max * (Gblood - ix1)


    # GLUCAGON
    dglucagon = 0
    if Gblood > 30:
        dglucagon += 1.3 * k_baseline_glucagon_secretion_max
    elif Gblood > 20:
        dglucagon += gm3 * k_baseline_glucagon_secretion_max * (Gblood - gx3)
    elif Gblood > 7.5:
        dglucagon += gm2 * k_baseline_glucagon_secretion_max * (Gblood - gx2)
    elif Gblood > 0:
        dglucagon += gm1 * k_baseline_glucagon_secretion_max * (Gblood - gx1)

        
    # SOMATOSTATIN
    dsomatostatin = sm1 * k_baseline_somatostatin_secretion_max * Gblood

    return (dinsulin, dglucagon, dsomatostatin)
