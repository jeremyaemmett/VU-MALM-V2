import parameters as params
import numpy as np
import system
import csv

# 2/27/2025_2


def read_row(row, row_num, var_num, n_layers, profiles, n_dimensions):

    if n_dimensions == 1: n_layers = 1
    for column_idx in range((var_num-1)*n_layers, var_num*n_layers):
        val = row[column_idx]
        val = val.replace('[', '').replace("'", '').replace(']', '')  # Clean-up non-numerical characters
        val = float(val)
        profiles[row_num, column_idx-(var_num-1)*n_layers] = val

    return(profiles)


def write_header(file, depths, keys):

    [file.write(depths[j] + ',') if j < len(depths)-1 else file.write(depths[j]) for j in range(0, len(depths))]
    file.write('\n')
    [file.write(keys[j] + ',') if j < len(keys) - 1 else file.write(keys[j]) for j in range(0, len(keys))]
    file.write('\n\n')

    return


def prepare_output(dict_dict, depths):

    depths = list(depths.astype(str))

    keys = list(dict_dict.keys())
    filenames = [keys[d] + '.txt' for d in range(0, len(keys))]

    for i in range(0, len(filenames)):

        fkeys = list(dict_dict[keys[i]].keys())

        with open(params.main_directory + 'output/' + filenames[i], "w") as file:
            write_header(file, depths, fkeys)

        file.close()

    return()


def write_output(count, c_datetime, chd, dtd, ifd, sfd):

    species_list = system.list_user_chemicals()

    with open(params.main_directory + 'output/chd.txt', "a") as file:
        file.write(f"{count},{c_datetime},{[','.join(map(str, chd[i]['conc'])) for i in species_list]}\n")

    with open(params.main_directory + 'output/ifd.txt', "a") as file:
        ifluxes = [ifd[i] for i in species_list]
        file.write(f"{count},{c_datetime},{','.join(map(str, ifluxes))}\n")

    with open(params.main_directory + 'output/sfd.txt', "a") as file:
        sfluxes = [sfd[i] for i in species_list]
        file.write(f"{count},{c_datetime},{','.join(map(str, sfluxes))}\n")


def read_output(filename, n_dimensions):

    # Get the number of time records, data key names, and depths (if applicable), from the file
    header_length = 3
    lines_in_file = open(filename, 'r').readlines()
    n_records = len(lines_in_file) - header_length
    with open(filename) as f:
        reader = csv.reader(f)
        depths = np.array(next(reader)).astype(float)
        reader = csv.reader(f)
        keys = list(next(reader))

    profiles = {}  # Dictionary to contain time record data for each data key name

    if n_dimensions == 2:  # Prepare a 2D array
        for i in range(0, len(keys)): profiles[keys[i]] = np.zeros((n_records, len(depths)))

    if n_dimensions == 1:  # Prepare a 1D array
        for i in range(0, len(keys)): profiles[keys[i]] = np.zeros((n_records, 1))

    with open(filename) as f:

        # Skip the header
        reader = csv.reader(f)
        for i in range(0, 3): next(reader)

        # Loop through the file rows (each row is a time record)
        r, date_times = 0, []
        for row in reader:
            date_times.append(row[1])
            row = row[2:]  # Extract the columns with data

            for c in range(0, len(keys)):
                profiles[keys[c]] = read_row(row, r, 1 + c, len(depths), profiles[keys[c]], n_dimensions)

            r = r + 1

        if n_dimensions == 1:
            for i in range(0, len(keys)): profiles[keys[i]] = profiles[keys[i]].T

    return profiles, depths, date_times