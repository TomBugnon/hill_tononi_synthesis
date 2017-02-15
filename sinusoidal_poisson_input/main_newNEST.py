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

import figure_4_plot_tom
reload(figure_4_plot_tom)

import figure_3_plot_newNEST
reload(figure_3_plot_newNEST)


sim_fig_3 = True

# vis_size = [10, 7]
vis_size = [40, 30]

for run in range(1,2):
#for run in range(1,21):

    # Which set of parameters for synaptic depression and NMDA synapse: old defaults (9/1/16) or New defaults (1/20.16)
    # for synapse_old in [True, False]:
    #     for NMDA_old in [True, False]:
    for synapse_old in [True]:
        for NMDA_old in [False]:

            if sim_fig_3:

                root_folder = '/Users/Tom/Documents/projects/network_simulations/HEAD_version/newNEST_singleccx'
                # root_folder = '/Users/Tom/Desktop/garbage'

                network = 'network_full_newNEST'

                ret_rate = 100.0

                edge_wrap = True
                net_folder = '/%s_edge_wrap_%d_Np_%d_Ns_%d_NMDA_%s_synapse_%s' % \
                              (network,
                               1*edge_wrap, vis_size[0], vis_size[1],
                               'old' if NMDA_old else 'new',
                               'old' if synapse_old else 'new')

                # vertical
                lambda_dg = 2.0
                #lambda_dg = 8.0
                input_flag = True
                data_folder = '/vertical_rate%d_run%d' % (int(ret_rate), run)


                Params = {
                    'show_main_figure' : False,
                    'start_membrane_potential' : 120.0,
                    'end_membrane_potential' : 130.0,
                    'show_center_connectivity_figure' : False,
                    'save_recorders' : True,
                    #'save_recorders' : False,
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

                    'lambda_dg': lambda_dg,  # visSize / number_of_lines
                    'input_flag': input_flag,

                    'root_folder': root_folder,
                    'net_folder': net_folder,
                    'data_folder': data_folder,

                    # for debugin
                    'dry_run': False,
                    #'dry_run': True

                    'plot_all_regions' : True,

                    'synapse_old' : synapse_old,
                    'NMDA_old' : NMDA_old,

                    'figure_title': ('Figure 3 from H&T2005. %s synaptic depression params, %s NMDA params' %
                                     ('old' if synapse_old else 'new',
                                      'old' if NMDA_old else 'new'))

                }

                # Run simulation of figure 3
                figure_3_plot_newNEST.simulation(Params)

