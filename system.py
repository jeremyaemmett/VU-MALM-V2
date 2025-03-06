""" System folder/file navigation, modification, and reading tasks """

import parameters as params
import conversions as conv
import numpy as np
import csv
import os

# 3/2/2025


def files_with_substring(folder, substring):

    return [folder + '/' + f for f in os.listdir(folder) if all(substr in f for substr in substring)]


def csv2array(filename):

    with open(filename, newline='') as csvfile: data = list(csv.reader(csvfile))

    return(np.array(data))


def datedprofiles2arrays(label, depths):

    root_files = files_with_substring(params.main_directory + 'data/' + params.site + '/sitedata', ['roo', 'AI'])
    days, profiles = [], []
    for f in root_files:
        d, v = csv2array(f)[1:, 0], csv2array(f)[1:, 1]
        d, v = [float(d[di]) for di in range(0, len(d))], [float(v[vi]) for vi in range(0, len(v))]
        xnew, root_profile_interp = conv.interp2grid(d, v, depths)
        days.append(conv.ddmmyyyy2doy(f.split('/')[-1].split('_')[1]))
        profiles.append(root_profile_interp)

    return(np.array(days), np.array(profiles))


def list_user_chemicals():

    all_chemicals = []
    reaction_keys = list(params.reactions.keys())
    for i in range(0, len(reaction_keys)):

        reaction_chemicals = list(params.reactions[reaction_keys[i]]['cons'].keys()) + \
                             list(params.reactions[reaction_keys[i]]['prod'].keys())
        all_chemicals += reaction_chemicals

    return np.sort(list(set(all_chemicals)))


def list_reaction_q10s():

    reactions = list(params.reactions.keys())
    print(reactions)
    for i in range(0, len(reactions)):
        print(reactions[i], params.reactions[reactions[i]]['q10']['val'])

    return


def list_reactions_involving_chemicals():

    unique_chemicals = list_user_chemicals()

    reactions_involving_chemicals = {}
    for j in range(0, len(unique_chemicals)):

        chemical = unique_chemicals[j]

        cons_reactions, prod_reactions = [], []
        reaction_keys = list(params.reactions.keys())

        for i in range(0, len(reaction_keys)):

            cons_chemicals = list(params.reactions[reaction_keys[i]]['cons'].keys())
            prod_chemicals = list(params.reactions[reaction_keys[i]]['prod'].keys())

            if chemical in cons_chemicals: cons_reactions.append(reaction_keys[i])
            if chemical in prod_chemicals: prod_reactions.append(reaction_keys[i])

        reactions_dictionary = {'cons': cons_reactions, 'prod': prod_reactions}
        reactions_involving_chemicals[chemical] = reactions_dictionary

    return reactions_involving_chemicals
