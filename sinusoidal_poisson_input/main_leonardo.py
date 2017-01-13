#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Simulation of figures 3 & 4 in Hill-Tononi paper.
# Author: Pablo Martinez Cañada (pablomc@ugr.es)
# Modified by Keiko Fujii and Leonardo S. Barbosa

import nest
import nest.topology as tp
import numpy as np
import time
import math

import figure_4_plot
reload(figure_4_plot)

import figure_3_plot_leonardo
reload(figure_3_plot_leonardo)

import test_scrambled_intact
reload(test_scrambled_intact)

sim_fig_3 = True
sim_fig_4 = False

#p_ratio = 1.
p_ratio = 2.

# vis_size = [10, 7]
vis_size = [40, 30]

for structured_input in [True, False]:
#for structured_input in [True]:
#for structured_input in [False]:
    for scramble in [True, False]:
    # for scramble in [True]:
    # for scramble in [False]:

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


            root_folder = '/data/nsdm'

            #network ='network_full_keiko'
            #network ='network_full_keiko2'
            network = 'network_full_leonardo'
            #network = 'network_full_leonardo2'

            # scramble network connections? only works with network_full_leonardo!
            # scramble = True
            # scramble = False

            # structured_input = True
            # structured_input = False
            ret_rate = 100.0

            edge_wrap = True
            #edge_wrap = False
            net_folder = '/%s_%s_edge_wrap_%d_Np_%d_Ns_%d_p_ratio_%d' % \
                          (network,
                           'scrambled' if scramble else 'intact',
                           1*edge_wrap, vis_size[0], vis_size[1], p_ratio)

            if structured_input:
                # vertical
                lambda_dg = 2.0
                #lambda_dg = 8.0
                input_flag = True
                data_folder = '/vertical_rate%d' % (int(ret_rate))
            else:
                lambda_dg = -1.0
                # input_flag = False
                input_flag = True
                data_folder = '/%s_rate%d' % (
                    'random' if input_flag else 'spontaneous',
                    int(ret_rate))

            Params = {
                'dump_connections' : False, # Takes a lot of disk space and time! half gigabyte...
                'load_connections' : False,   # Load connections from files GENERATED IN A PREVIOUS RUN
                'show_main_figure' : False,
                'start_membrane_potential' : 120.0,
                'end_membrane_potential' : 130.0,
                'show_V4_num_conn_figure' : True,
                'show_V4_connectivity_figure' : False,
                'show_center_connectivity_figure' :False,
                'network' : network,
                #'Np': 8,
                #'Ns': 4,
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

                'edge_wrap' : edge_wrap,

                'scrambled' : scramble, # scramble the connectivity: no invariance for horizontal/vertical stimulus

                'lambda_dg': lambda_dg,  # visSize / number_of_lines
                'input_flag': input_flag,

                'root_folder': root_folder,
                'net_folder': net_folder,
                'data_folder': data_folder,

                # for debugin
                'p_ratio': p_ratio,
                #'dry_run': True
                'dry_run': True

            }

            # Run simulation of figure 3
            figure_3_plot_leonardo.simulation(Params)
            # test_scrambled_intact.simulation(Params)


