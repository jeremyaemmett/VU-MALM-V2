""" Calculate rates needed to nudge chemical concentrations to equilibrium values in sub-saturated soil layers """

def transport(species, chd, atd, forcing_t, dt, watidx):

    U1 = chd[species]['conc']
    U2 = forcing_t['eqc'][species]

    rates = 0.0 * U1
    rates[0:(watidx-1)] = (U2[0:(watidx-1)] - U1[0:(watidx-1)]) / dt  # Rates
    atd[species] = rates

    return atd
