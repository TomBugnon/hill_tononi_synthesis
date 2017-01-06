#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Simulation of figures 3 & 4 in Hill-Tononi paper.
# Author: Pablo Martinez Cañada (pablomc@ugr.es)

import nest
import nest.topology as tp
import numpy as np
import time
import math

import figure_4_plot
reload(figure_4_plot)

import figure_3_plot
reload(figure_3_plot)

sim_fig_3 = True
sim_fig_4 = False

if sim_fig_4:

    Params = {
        'Np': 40, # cells in the primary visual area
        'Ns': 30, # cells in the secondary visual area
        'visSize': 8.0, # visual angle (degrees)
        'ret_rate': 30.0, # average firing rate of retina ganglion cells (spikes^(-1))
        'ret_amplitude': 0.0, # amplitude of the sinusoidal poisson generator
                              # used for retina ganglion cells (spikes^(-1))
        'temporal_frequency': 0.0, # frequency of the generator (Hz)
        'threads': 12, # threads to use in NEST simulation
        #'intervals': [1000.0, 1000.0, 7500.0],  # keiko
        #'intervals': [500.0],  # Intervals (in ms) of the waking,transition
        'intervals': [500.0, 500.0, 500.0, 500.0, 3000.0],  # Intervals (in ms) of the waking,transition
                                           # and sleep modes
        'resolution': 1.0 # Simulation step (in ms)
    }

    # Run simulation of figure 4
    figure_4_plot.simulation(Params)



if sim_fig_3:

    root_data_folder = '/data/nsdm/'

    network ='network_full_keiko'
    # network = 'network_full_leonardo'

    # scramble network connections? only works with network_full_leonardo!
    scramble = True
    # scramble = False

    # structured_input = True
    structured_input = False
    ret_rate = 100.0

    if structured_input:
        # vertical
        lambda_dg = 2.0
        input_flag = True
        data_folder = '%s/vertical_rate%d_%s_%s/' % (root_data_folder, int(ret_rate), network, 'scrambled' if scramble else 'intact')
    else:
        lambda_dg = -1.0
        input_flag = False
        data_folder = '%s/random_rate%d_%s_%s/' % (root_data_folder, int(ret_rate), network, 'scrambled' if scramble else 'intact')

    Params = {
        'network' : network,
        'Np': 40,
        'Ns': 30,
        'visSize': 8.0,
        'ret_rate': ret_rate,#20.0,
        'ret_amplitude': 0.0, # random
        'temporal_frequency': 2.0, # (Hz)
        'spatial_frequency' : 0.5, # (cpd)
        'threads': 12,
        #'intervals': [100.0, 250.0, 650.0],  # original
        #'intervals': [5000.0],  # keiko
        'intervals': [2000.0],  # leonardo
        'resolution': 1.0,
        'phi_dg': 0.0,  # vertical
        #'phi_dg': 0.5*np.pi, # horizontal

        'scrambled' : scramble, # scramble the connectivity: no invariance for horizontal/vertical stimulus
        'load_connections_from_file': './scrambled_connection_Tp_Cortex.mat',

        'lambda_dg': lambda_dg,  # visSize / number_of_lines
        'input_flag': input_flag,

        'data_folder': data_folder
    }

    # Run simulation of figure 3
    figure_3_plot.simulation(Params)


