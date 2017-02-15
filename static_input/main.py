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

import static_figure
reload(static_figure)

# import test_scrambled_intact
# reload(test_scrambled_intact)


# sim_fig_3 = True
# sim_fig_4 = False

#p_ratio = 1.
p_ratio = 2.

# vis_size = [10, 7]
vis_size = [40, 30]

for run in range(1,2):


    root_folder = '/Users/Tom/Documents/projects/network_simulations/HEAD_version/full/'
    root_folder = root_folder + 'static_images/'

    # root_folder = '/Users/Tom/Desktop/garbage'

    #network ='network_full_keiko'
    #network ='network_full_keiko2'
    # network = 'network_full_leonardo'
    network = 'network_full_tom'
    #network = 'network_full_leonardo2'

    # scramble network connections? only works with network_full_leonardo!
    # scramble = True
    # scramble = False

    # structured_input = True
    # structured_input = False


    ret_rate = 100.0

    synapse_keiko = True
    NMDA_keiko = True


    edge_wrap = True
    #edge_wrap = False
    net_folder = '/%s_edge_wrap_%d_Np_%d_Ns_%d_p_ratio_%d_NMDA_%s_synapse_%s' % \
                 (network,
                  1*edge_wrap, vis_size[0], vis_size[1], p_ratio,
                  'keiko' if NMDA_keiko else 'default',
                  'keiko' if synapse_keiko else 'default')

    new_image = True #If false, loads image_folder+image_name. if True, uses new_image_elements to randomly create new image and saves it in image_folder+image_name
    image_folder = '/Users/Tom/Documents/projects/network_simulations/'
    image_name = ''
    new_image_elements = {'horizontal':{'size':3, 'number':1}, 'vertical':{'size':2, 'number':1}, 'cross':{'size':2, 'number':1}} #Elements of the created image if none is specified

    data_folder = '/%s_rate%d_run%d' % (image_name,int(ret_rate), run)


    Params = {
        'dump_connections' : False, # Takes a lot of disk space and time! half gigabyte...
        'load_connections' : False,   # Load connections from files GENERATED IN A PREVIOUS RUN
        'show_main_figure' : False,
        'start_membrane_potential' : 120.0,
        'end_membrane_potential' : 140.0,
        'show_V4_num_conn_figure' : False,
        'show_V4_connectivity_figure' : False,
        'show_center_connectivity_figure' : False,
        'save_recorders' : True,
        #'save_recorders' : False,
        'network' : network,
        #'Np': 8,
        #'Ns': 4,
        'Np': 40,
        'Ns': 30,
        'visSize': 8.0,

        'threads': 12,
        #'intervals': [100.0, 250.0, 650.0],  # original
        #'intervals': [5000.0],  # keiko
        'initiation_time' : 100.,
        'intervals': [100.0],  # leonardo

        'resolution': 1.0,

        'edge_wrap' : edge_wrap,

        'root_folder': root_folder,
        'net_folder': net_folder,
        'data_folder': data_folder,

        # for debugin
        'p_ratio': p_ratio,
        'dry_run': False,
        #'dry_run': True

        'plot_all_regions' : True,
        'plot_static_all_regions' : True,

        'synapse_keiko' : synapse_keiko,
        'NMDA_keiko' : NMDA_keiko,

        #input image:
        'image_folder' : image_folder,
        'image_name' : image_name,
        'new_image' : new_image,
        'new_image_elements' : new_image_elements,

        #retina params:
        'ret_rate': ret_rate,  # rate = luminance * ret_rate + baseline_ret_rate
        'ret_rate_baseline':20

    }

    # Run simulation of figure 3
    static_figure.simulation(Params)
    # test_scrambled_intact.simulation(Params)


