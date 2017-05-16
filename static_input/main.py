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


sim_fig_3 = True

# vis_size = [10, 7]
vis_size = [40, 30]
# vis_size = [80, 60]

for run in range(1,2):

    root_folder = '/Users/Tom/Documents/projects/network_simulations/HEAD_version/full/'
    root_folder = root_folder + 'static_images/'

    network = 'network_full_newNEST_modify'

    NMDA_old = False
    synapse_old = True


    specific_name = 'all_horiz_intralaminar_1.5'

    max_ret_rate = 400.0

    edge_wrap = True
    net_folder = '/%s_edge_wrap_%d_Np_%d_Ns_%d_NMDA_%s_synapse_%s' % \
                  (network,
                   1*edge_wrap, vis_size[0], vis_size[1],
                   'old' if NMDA_old else 'new',
                   'old' if synapse_old else 'new')

    new_image = False #If false, loads image_folder+image_name. if True, uses new_image_elements to randomly create new image and saves it in image_folder+image_name
    image_folder = root_folder + 'images/'
    image_name = 'horiz_vert_cross' # Don't add the extension (.pickle)
    new_image_elements = {'horizontal':{'size':(15,4), 'number':1}, 'vertical':{'size':(4,15), 'number':1}, 'cross':{'size':(9,3), 'number':0}} #Elements of the created image if none is specified

    data_folder = '/%s_rate%d_run%d_%s' % (image_name,int(max_ret_rate), run, specific_name)


    Params = {
        'dump_connections' : False, # Takes a lot of disk space and time! half gigabyte...
        'load_connections' : False,   # Load connections from files GENERATED IN A PREVIOUS RUN
        'show_main_figure' : False,
        'start_membrane_potential' : 0.0,
        'end_membrane_potential' : 150.0,
        'show_V4_num_conn_figure' : False,
        'show_V4_connectivity_figure' : False,
        'show_center_connectivity_figure' : False,

        'save_recorders' : False,
        #'save_recorders' : False,
        'network' : network,
        #'Np': 8,
        #'Ns': 4,
        # 'Np': 80,
        # 'Ns': 60,
        # 'visSize': 16.,
        'Np': 40,
        'Ns': 30,
        'visSize': 8.,

        'threads': 12,
        #'intervals': [100.0, 250.0, 650.0],  # original
        #'intervals': [5000.0],  # keiko
        'initiation_time' : 0.,
        'intervals': [200.0],

        'resolution': 1.0,

        'edge_wrap' : edge_wrap,

        'root_folder': root_folder,
        'net_folder': net_folder,
        'data_folder': data_folder,

        # for debugin
        'dry_run': False,
        #'dry_run': True

        'plot_all_regions' : False,
        'plot_topo_all_regions' : True,

        'synapse_old' : synapse_old,
        'NMDA_old' : NMDA_old,

        #input image:
        'image_folder' : image_folder,
        'image_name' : image_name,
        'new_image' : new_image,
        'new_image_elements' : new_image_elements,

        #retina params:
        'ret_rate': max_ret_rate,  # rate(i) = ret_rate if luminance = 1 else min_rate
        'ret_rate_baseline':50., #Unused if automatic_min_rate = True, in which case min_rate = ret_rate_baseline
        'automatic_min_rate':True, #If true, determines min_rate ('off' pixels) so that mean rate of whole retina  is mean_rate
        'mean_rate':70.


    }

    # Run simulation of figure 3
    static_figure.simulation(Params)
    # test_scrambled_intact.simulation(Params)


