""" Main code, containing the computational time-loop and major routine calls """

import matplotlib.pyplot as plt
from datetime import timedelta
import parameters as params
import initialize as init
import conversions
import numpy as np
import atmosphere
import diffusion
import forcing
import newstep
import output
import system
import plots2
import plots
import time

# 3/17/2025

thicks, depths = init.define_layers()

fd = forcing.forcing_data()  # Forcing variables
#plots.forcing_plots(fd)
chd, mmd, atd, dtd, ptd, ifd, sfd, dict_dict = init.variables()  # Chemistry, diffusion, and plant transport variables

species_list = system.list_user_chemicals()

#fig, axes = plt.subplots(1, len(species_list), figsize=(15, 5))
#cmap = plt.get_cmap('rainbow')

# Main computational time loop
dt, count, elapsed, doy = params.dt, 0, 0.0, conversions.datetime2doy(params.start_datetime)
run_time, start = (params.end_datetime - params.start_datetime).days, time.time()
watidx, iceidx, watidx_old, iceidx_old = 0, 0, 0, 0
if params.jump2plots == False: output.prepare_output(dict_dict, depths)
while elapsed < run_time and params.jump2plots == False:

    current_datetime = params.start_datetime + timedelta(days=elapsed)
    print_flag = current_datetime.hour == 0 and current_datetime.minute == 0
    write_flag = elapsed % params.write_dt == 0
    if print_flag:
        print(current_datetime)
        #print(np.dot(chd['ch4']['conc'], depths))
        #print(chd['ch4']['conc'])
    #    print(current_datetime)

    # Slice the forcing data at the current day-of-year (doy)
    forcing_t = forcing.forcing_data_t(fd, doy)

    # Ice table depth
    iceidx = 0 if np.max(forcing_t['tem']) < 0.0 else np.where(forcing_t['tem'] <= 0.0)[0][0]
    watidx = 0 if np.max(forcing_t['tem']) < 0.0 else np.where(forcing_t['sat'] == 100.0)[0][0]
    level_change = watidx != watidx_old or iceidx != iceidx_old
    watidx_old, iceidx_old = watidx, iceidx
    iceline, watline = depths[iceidx], depths[watidx]
    #print('ch4 before processes', np.dot(chd['ch4']['conc'][0:watidx], thicks[0:watidx]))
    if print_flag:
        print('ch4 before processes', np.dot(chd['ch4']['conc'], thicks))
        #print(np.dot(chd['ch4']['conc'][0:watidx], thicks[0:watidx]))
        #print(np.dot(chd['ch4']['conc'][watidx:], thicks[watidx:]))

    if iceline != 0.0: # Call processes only if the ground isn't completely frozen

        # Plot variables
        norm = elapsed / run_time
        lwidth = 3.0 if elapsed == 0 else 1.0

        # Loop through the chemical species
        for s in range(0, len(species_list)):

            species = species_list[s]

            # Fill unsaturated layers with equilibrium concentrations
            if watidx > 1 and params.instant_diffusion: chd[species]['conc'][0:watidx] = forcing_t['eqc'][species][0:watidx]

            # Calculate fill-in rates for sub-sat. layers (diff. between equilibrium and last-updated concentrations)
            if watidx > 1: atd = atmosphere.transport(species, chd, atd, forcing_t, dt, watidx)

            #stop

            #axes[s].set_title(species)
            #axes[s].set_facecolor('gainsboro')
            #axes[s].set_ylim([1.0, 0.0])
            #axes[s].set_xlabel('Concentration')
            #axes[s].set_ylabel('Depth') if s == 0 else axes[s].set_ylabel(' ')
            #axes[s].grid(True)

            ###if print_flag:
            ###    if species == 'ch4': print('ch4 before diffusion process', np.dot(chd['ch4']['conc'], thicks))
            #    axes[s].step(cd['chem'][species]['conc'], depths, linewidth=lwidth, color=cmap(norm))
                #axes[s].plot([0.0, 1.0], [iceline, iceline], linewidth=lwidth, linestyle = '--', color='gray', alpha = 0.5)
                #axes[s].plot([0.0, 1.0], [watline, watline], linewidth=lwidth, linestyle='-', color='gray', alpha = 0.5)
            # Calculate gas diffusion rates through soil (mol/m2/day for each layer)
            dtd[species]['prof'] *= 0.0
            if params.diffusion_flag: dtd, ifd = diffusion.transport(species, chd, dtd, ifd, forcing_t, dt, level_change)
            ###if species == 'ch4': print('ch4 after diffusion process', np.dot(chd['ch4']['conc'], thicks))
            # Diffusion fill-in partitioning
            #if params.instant_diffusion and watidx > 1: dtd[species]['prof'][0:watidx] += atd[species][0:watidx]

            # Calculate gas transport rates through plants (mol/m2/day for each layer)
            # Plant fill-in partitioning
            if params.instant_diffusion and watidx > 1: ptd[species]['prof'][0:watidx] += 0.0 * atd[species][0:watidx]

            # Update concentration fields (mol/m3; microbeC/m3) and calculate surface fluxes (mol/m2/day)
            subsat_sum_1 = np.dot(chd[species]['conc'][0:(watidx-1)], thicks[0:(watidx-1)])
            #if species == 'ch4' and print_flag:
            #    print('subsat_sum_1: ', subsat_sum_1)
            #    print(dtd[species]['prof'][0:(watidx-1)])
            ###if species == 'ch4': print('ch4 before newstep', np.dot(chd['ch4']['conc'], thicks))
            chd = newstep.newstep(species, chd, dtd, dt)
            ###if species == 'ch4': print('ch4 after newstep', np.dot(chd['ch4']['conc'], thicks))
            subsat_sum_2 = np.dot(chd[species]['conc'][0:(watidx-1)], thicks[0:(watidx-1)])
            if species == 'ch4' and print_flag: print('subsat_sum_2: ', subsat_sum_2)
            #if species == 'co2':
                #print(' ')
                #print(atd[species][0:(watidx - 1)])
                #print(dtd[species]['prof'][0:(watidx-1)])
            #sfd[species] = np.dot(forcing_t['eqc'][species][0:(watidx-1)] - chd[species]['conc'][0:(watidx-1)],
            sfd[species] = (subsat_sum_2 - subsat_sum_1) / dt if level_change == False else 0.0
            #if species == 'ch4' and print_flag:
            #    print(sfd[species], ifd[species])

    if write_flag: output.write_output(count, current_datetime, chd, dtd, ifd, sfd)

    # Update the elapsed time, day-of-year, and time-step
    count += 1
    dt = newstep.adaptive_timestep(forcing_t, elapsed)
    elapsed, doy = round(elapsed + dt, 2), (doy + dt) % 365

    #stop

end = time.time()
print('Runtime: ', end - start)

# dict_name, dict_dict, fd, log_flag, fsize, n_dimensions, flip_y_flag, symmetric_y_flag, resampling_flag
test = plots2.plots('chd', dict_dict, fd, True, 48, 2, False, False, False)
test = plots2.plots('ifd', dict_dict, fd, True, 48, 1, True, True, True)
test = plots2.plots('sfd', dict_dict, fd, True, 48, 1, True, True, True)

#plt.show()

