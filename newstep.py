import parameters as params
import numpy as np

# 3/17/2025

def newstep(species, chd, dtd, dt):

    """ Update dynamic quantities using their predicted rates of change over the current timestep
    :param species:
    :param chd:
    :param dtd:
    :param dt:
    :return:
    """

    chd[species]['conc'] += dtd[species]['prof'] * dt

    chd[species]['conc'][chd[species]['conc'] < 2e-10] = 1e-10  # Suppress numerical noise around the minimum value

    return(chd)


def adaptive_timestep(forcingt, elapsed):

    """ Update the current time-step for e.g. stability, efficiency. Alignment with write time-steps is ensured when the
    timestep would otherwise over-step these, to prevent skipping output records
    :param forcingt:
    :param elapsed:
    :return:
    """

    best_dt = params.dt

    if max(forcingt['tem']) <= 0.0: best_dt = params.write_dt # Proceed at the write time-step if all layers are frozen

    if max(forcingt['tem']) > 0.0:
        best_dt = params.dt
        # If the new dt would overshoot the next write interval, and the

    if elapsed % params.write_dt + best_dt > params.write_dt and params.write_dt - elapsed % params.write_dt > 1e-10:
        best_dt = round(params.write_dt - elapsed % params.write_dt, 2)

    #print(' best dt: ', best_dt, ' elapsed: ', elapsed)

    #stop

    return(best_dt)