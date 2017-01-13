#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Simulation of figure 3 of Hill-Tononi paper
# Author: Pablo Martinez Cañada (pablomc@ugr.es)
#
# Evoked activity: the stimulus is ON during the second interval.
# I used the same parameters of Figure 4 but I had to change gKL to 0.8
# to get a little stronger evoked response.

import nest
import nest.topology as tp
import numpy as np
import numpy.random as rd
import os
import matplotlib.pyplot as plt

import plotting
reload(plotting)

import pickle
import os.path
import scipy.io
import pylab
from nest import raster_plot

import math

import shutil


def phaseInit(pos, lam, alpha):
    '''Initializer function for phase of drifting grating nodes.

       pos  : position (x,y) of node, in degree
       lam  : wavelength of grating, in degree
       alpha: angle of grating in radian, zero is horizontal

       Returns number to be used as phase of sinusoidal Poisson generator.
    '''
    return 360.0 / lam * (math.cos(alpha) * pos[0] + math.sin(alpha) * pos[1])


def simulation(Params):


    #! =================
    #! Import network
    #! =================

    # NEST Kernel and Network settings
    nest.ResetKernel()
    nest.ResetNetwork()
    nest.SetKernelStatus({"local_num_threads": Params['threads'],'resolution': Params['resolution']})
    nest.SetStatus([0],{'print_time': True})

    # initialize random seed
    import time
    msd = int(round(time.time() * 1000))
    nest.SetKernelStatus({'grng_seed' : msd})
    nest.SetKernelStatus({'rng_seeds' : range(msd+Params['threads']+1, msd+2*Params['threads']+1)})


    import importlib
    network = importlib.import_module(Params['network'])
    reload(network)
    models, layers, conns  = network.get_Network(Params)
    #import network_full_keiko
    #reload(network_full_keiko)
    # models, layers, conns  = network_full_keiko.get_Network(Params)

    # Create models
    for m in models:
            nest.CopyModel(m[0], m[1], m[2])

    # Create layers, store layer info in Python variable
    for l in layers:
            exec '%s = tp.CreateLayer(l[1])' % l[0]

    # Create connections, need to insert variable names
    for c in conns:
            eval('tp.ConnectLayers(%s,%s,c[2])' % (c[0], c[1]))


    # Prepare for file IO
    import glob
    # --- Set folder information
    data_folder = Params['data_folder']
    if not os.path.isdir(data_folder):
        os.makedirs(data_folder)

    # --- To save spike data, set pairs of population id and its name
    population_name = [ {'population': Retina_layer, 'name': 'Retina'},
                        {'population': Vp_vertical, 'name': 'Vp_v'},
                        {'population': Vp_horizontal, 'name': 'Vp_h'},
                        {'population': Rp_layer, 'name': 'Rp'},
                        {'population': Tp_layer, 'name': 'Tp'},
                        {'population': Vs_vertical, 'name': 'Vs_v'},
                        {'population': Vs_horizontal, 'name': 'Vs_h'}]


    if Params.has_key('load_connections_from_file') and Params['load_connections_from_file']:

        # Preparation
        scramble_populations = [(Vp_vertical, 'Vp_vertical'),
                                (Vp_horizontal, 'Vp_horizontal')]
        scramble_layers = ['L23_exc',
                           'L23_inh',
                           'L4_exc',
                           'L4_inh',
                           'L56_exc',
                           'L56_inh']
        #scramble_layers = ['L4_exc']

        # Get min &max index of each layer
        h_min_idx = {}
        h_max_idx = {}
        v_min_idx = {}
        v_max_idx = {}
        target_neurons = []
        for tmp_model in scramble_layers:
            tmp_h = nest.GetLeaves(Vp_horizontal, properties={'model': tmp_model}, local_only=True)[0]
            tmp_v = nest.GetLeaves(Vp_vertical, properties={'model': tmp_model}, local_only=True)[0]
            h_min_idx[tmp_model] = min(tmp_h)
            h_max_idx[tmp_model] = max(tmp_h)
            v_min_idx[tmp_model] = min(tmp_v)
            v_max_idx[tmp_model] = max(tmp_v)
            target_neurons = target_neurons + range(h_min_idx[tmp_model], h_max_idx[tmp_model]+1) + range(v_min_idx[tmp_model], v_max_idx[tmp_model]+1)

        # Save intact network information
        for p in range(0, len(population_name), 1):

            population = population_name[p]['population']
            p_name = population_name[p]['name']
            filename_AMPA = data_folder + 'connection_' + p_name + '_AMPA_syn' + '_intact.dat'
            filename_NMDA = data_folder + 'connection_' + p_name + '_NMDA_syn' + '_intact.dat'
            filename_GABAA = data_folder + 'connection_' + p_name + '_GABA_A_syn' + '_intact.dat'
            filename_GABAB = data_folder + 'connection_' + p_name + '_GABA_B_syn' + '_intact.dat'
            tp.DumpLayerConnections(population, 'AMPA_syn', filename_AMPA)
            tp.DumpLayerConnections(population, 'NMDA_syn', filename_NMDA)
            tp.DumpLayerConnections(population, 'GABA_A_syn', filename_GABAA)
            tp.DumpLayerConnections(population, 'GABA_B_syn', filename_GABAB)

        # Reset network
        nest.ResetNetwork()

        '''
        # The following code works, but takes time longer

        # Get the list of connection_file
        file_list = glob.glob(data_folder+'connection_*_intact.dat')

        # Remove a file if the size is zero
        remove_idx = []
        for file_idx in range(0,len(file_list)):
            intact_filename = file_list[file_idx]
            fsize=os.path.getsize(intact_filename)
            if fsize == 0:
                remove_idx = remove_idx + [file_idx]

        remove_idx.sort()
        remove_idx.reverse()
        for i in remove_idx:
            del file_list[i]
        '''

        # TODO : put GABA_A, GABA_B and NMDA connection files
        file_list = []
        file_list.append({'filename': data_folder + 'connection_Vp_h_AMPA_syn_intact.dat',
                           'synapse': 'AMPA'})
        file_list.append({'filename': data_folder + 'connection_Vp_v_AMPA_syn_intact.dat',
                           'synapse': 'AMPA'})


        # Do the following procedure for all connection files
        for file_idx in range(0,len(file_list)):

            # Set filenames
            intact_filename = file_list[file_idx]['filename']
            receptors = nest.GetDefaults('ht_neuron')['receptor_types']
            syn_model = receptors[ file_list[file_idx]['synapse'] ]

            scrambled_filename = intact_filename.replace('intact', 'scrambled')
            print(intact_filename)

            # Get original(intact) connectivity data
            src_network = np.loadtxt(open(intact_filename,'rb'))
            np_pre = src_network[:, 0].astype(int)
            np_post = src_network[:, 1].astype(int)
            np_w = src_network[:, 2]
            np_d = src_network[:, 3]

            if Params['scrambled']:

                # Preserve the original structure if
                # -- pre neruons are not in target populations (i.e. scramble_layers)
                # -- OR
                # -- post neurons are not in target populations(i.e. scramble_layers)
                preserved_rows = np.where( ~np.in1d(np_pre,target_neurons) | ~np.in1d(np_post,target_neurons) )[0]
                preserved_pre = np_pre[preserved_rows]
                preserved_post = np_post[preserved_rows]
                preserved_w = np_w[preserved_rows]
                preserved_d = np_d[preserved_rows]

                # If len(preserved_rows)==len(np_pre), that means all neurons do not belong to scramble target areas.
                # If len(preserved_rows) < len(np_pre), that means there are some neurons which need to be scrambled.
                if len(preserved_rows) > len(np_pre):
                    print('ERROR: preserved_rows should not be larger than np_pre')

                elif len(preserved_rows) == len(np_pre):
                    scrambled_pre = preserved_pre.tolist()
                    scrambled_post = preserved_post.tolist()
                    scrambled_w = preserved_w.tolist()
                    scrambled_d = preserved_d.tolist()

                else:  # --- len(preserved_rows) < len(np_pre)

                    scrambled_pre = []
                    scrambled_post = []
                    scrambled_w = []
                    scrambled_d = []

                    for tmp_model_pre in scramble_layers:

                        for tmp_model_post in scramble_layers:

                            # Get row index such that
                            # --- pre_neuron is in "tmp_model_pre"
                            # --- AND
                            # --- post_neuron is in "tmp_model_post"
                            bool_pre_h = np.in1d(np_pre, range(h_min_idx[tmp_model_pre], h_max_idx[tmp_model_pre]+1))
                            bool_pre_v = np.in1d(np_pre, range(v_min_idx[tmp_model_pre], v_max_idx[tmp_model_pre]+1))
                            bool_post_h = np.in1d(np_post, range(h_min_idx[tmp_model_post], h_max_idx[tmp_model_post]+1))
                            bool_post_v = np.in1d(np_post, range(v_min_idx[tmp_model_post], v_max_idx[tmp_model_post]+1))

                            tmp_rows_pre_h = np.where(bool_pre_h & (bool_post_h|bool_post_v))[0]
                            tmp_rows_pre_v = np.where(bool_pre_v & (bool_post_h|bool_post_v))[0]

                            # Get connectivity data such that pre neuron is in "tmp_model_pre"
                            # --- pre, w and d information should be kept.
                            # --- post index should be scrambled.
                            tmp_pre = np_pre[np.append(tmp_rows_pre_h, tmp_rows_pre_v)].tolist()
                            tmp_w = np_w[np.append(tmp_rows_pre_h, tmp_rows_pre_v)].tolist()
                            tmp_d = np_d[np.append(tmp_rows_pre_h, tmp_rows_pre_v)].tolist()
                            tmp_post = []

                            num_pre_h = len(tmp_rows_pre_h)
                            num_pre_v = len(tmp_rows_pre_v)

                            # --- pre : population = horizontal, model = "tmp_model_pre"
                            # --- post: population = horizontal(1/2)+vertical(1/2), model = "tmp_model_post"
                            if num_pre_h > 0:
                                # Assign the same number of connections for horizontal population and vertical population
                                num_scrambled_h = int(round(num_pre_h / 2))
                                num_scrambled_v = num_pre_h - num_scrambled_h

                                # Choose post neuron index randomly
                                scrambled_h_idx = rd.randint(low =h_min_idx[tmp_model_post],
                                                             high=h_max_idx[tmp_model_post],
                                                             size=[num_scrambled_h, 1])
                                scrambled_v_idx = rd.randint(low =v_min_idx[tmp_model_post],
                                                             high=v_max_idx[tmp_model_post],
                                                             size=[num_scrambled_v, 1])

                                # append scrambled post index
                                tmp_post = tmp_post + np.append(scrambled_h_idx, scrambled_v_idx).tolist()

                            # --- pre : population = vertical, model = "tmp_model_pre"
                            # --- post: population = horizontal(1/2)+vertical(1/2), model = "tmp_model_post"
                            if num_pre_v > 0:
                                # Assign the same number of connections for horizontal population and vertical population
                                num_scrambled_h = int(round(num_pre_v / 2))
                                num_scrambled_v = num_pre_v - num_scrambled_h

                                # Choose post neuron index randomly
                                scrambled_h_idx = rd.randint(low =h_min_idx[tmp_model_post],
                                                             high=h_max_idx[tmp_model_post],
                                                             size=[num_scrambled_h, 1])
                                scrambled_v_idx = rd.randint(low =v_min_idx[tmp_model_post],
                                                             high=v_max_idx[tmp_model_post],
                                                             size=[num_scrambled_v, 1])

                                # append scrambled post index
                                tmp_post = tmp_post + np.append(scrambled_h_idx, scrambled_v_idx).tolist()

                            scrambled_pre = scrambled_pre + tmp_pre
                            scrambled_post = scrambled_post + tmp_post
                            scrambled_w = scrambled_w + tmp_w
                            scrambled_d = scrambled_d +tmp_d

                    # append preserved connection data
                    scrambled_pre = scrambled_pre + preserved_pre.tolist()
                    scrambled_post = scrambled_post + preserved_post.tolist()
                    scrambled_w = scrambled_w + preserved_w.tolist()
                    scrambled_d = scrambled_d + preserved_d.tolist()

                    # Save scrambled_data
                    scrambled_all_data = np.zeros([len(scrambled_pre), 4])
                    scrambled_all_data[:, 0] = scrambled_pre
                    scrambled_all_data[:, 1] = scrambled_post
                    scrambled_all_data[:, 2] = scrambled_w
                    scrambled_all_data[:, 3] = scrambled_d
                    np.savetxt(scrambled_filename, scrambled_all_data, fmt='%.6f')

                # Connect
                con_dict = {'rule': 'one_to_one'}
                syn_dict = {"model": "ht_synapse",
                            'receptor_type': syn_model,
                            'weight': scrambled_w,
                            'delay': scrambled_d,
                            }
                nest.Connect(scrambled_pre, scrambled_post, con_dict, syn_dict)


            else:
                # just convert data type(ndarray->list) and connect based on the original data
                pre = np_pre.tolist()
                post = np_post.tolist()
                w = np_w.tolist()
                d = np_d.tolist()

                # Connect
                con_dict = {'rule': 'one_to_one'}
                syn_dict = {"model": "ht_synapse",
                            'receptor_type': 1,
                            'weight': w,
                            'delay': d,
                            }
                nest.Connect(pre, post, con_dict, syn_dict)

    # nest.DisconnectOneToOne(tp_node, tgt_map[0], {"synapse_model": "AMPA_syn"})
    #nest.Disconnect([tp_node], tgt_map, 'one_to_one', {"synapse_model": "AMPA_syn"})

    # Get target nodes for the vertical population
    # tp_nodes = nest.GetLeaves(Tp_layer, local_only=True)[0]

    if Params.has_key('show_V4_num_conn_figure') and Params['show_V4_num_conn_figure']:

        horizontal_nodes = nest.GetLeaves(Vp_horizontal, properties={'model': 'L4_exc'}, local_only=True)[0]
        vertical_nodes = nest.GetLeaves(Vp_vertical, properties={'model': 'L4_exc'}, local_only=True)[0]

        n_conns_hor = []
        for (idx, horizontal_node) in enumerate(horizontal_nodes):
            tgt_map = []
            this_conns = nest.GetConnections([horizontal_node], horizontal_nodes, synapse_model='AMPA_syn')
            tgt_map.extend([conn[1] for conn in this_conns])
            n_conns_hor.append(len(tgt_map))

        plt.figure()
        plt.hist(n_conns_hor, range(0, max(n_conns_hor + [30])))
        plt.title('# of connections Vp(h) L4pyr -> Vp(h) L4Pyr')

        n_conns_hor = []
        for (idx, horizontal_node) in enumerate(horizontal_nodes):
            tgt_map = []
            this_conns = nest.GetConnections([horizontal_node], vertical_nodes, synapse_model='AMPA_syn')
            tgt_map.extend([conn[1] for conn in this_conns])
            n_conns_hor.append(len(tgt_map))

            # nest.DisconnectOneToOne(tp_node, tgt_map[0], {"synapse_model": "AMPA_syn"})
            #nest.Disconnect([tp_node], tgt_map, 'one_to_one', {"synapse_model": "AMPA_syn"})

        plt.figure()
        plt.hist(n_conns_hor, range(0, max(n_conns_hor + [30])))
        plt.title('# of connections Vp(h) L4pyr -> Vp(v) L4Pyr')

        n_conns_ver = []
        for (idx, vertical_node) in enumerate(vertical_nodes):
            tgt_map = []
            this_conns = nest.GetConnections([vertical_node], vertical_nodes, synapse_model='AMPA_syn')
            tgt_map.extend([conn[1] for conn in this_conns])
            n_conns_ver.append(len(tgt_map))

        plt.figure()
        plt.hist(n_conns_ver, range(0, max(n_conns_ver + [30])))
        plt.title('# of connections Vp(v) L4pyr -> Vp(v) L4Pyr')

        n_conns_ver = []
        for (idx, vertical_node) in enumerate(vertical_nodes):
            tgt_map = []
            this_conns = nest.GetConnections([vertical_node], horizontal_nodes, synapse_model='AMPA_syn')
            tgt_map.extend([conn[1] for conn in this_conns])
            n_conns_ver.append(len(tgt_map))

        plt.figure()
        plt.hist(n_conns_ver, range(0, max(n_conns_ver + [30])))
        plt.title('# of connections Vp(v) L4pyr -> Vp(h) L4Pyr')

        plt.show()

    # Check connections

    # Connections from Retina to TpRelay
    # tp.PlotTargets(tp.FindCenterElement(Retina_layer), Tp_layer)

    if Params.has_key('show_V4_connectivity_figure') and Params['show_V4_connectivity_figure']:

        Vp_hor_gids = tp.GetElement(Vp_horizontal, [0,0])
        n_Vp_hor = len(Vp_hor_gids)

        f = []
        for idx in range(n_Vp_hor):
            f.append(plt.figure())

        positions = range(0,41,10)
        positions[-1] -= 1
        for xi in range(len(positions)):
            for yi in range(len(positions)):
                print("Position [%d,%d] : %d" %(xi,yi,yi*(len(positions))+xi+1))
                x = positions[xi]
                y = positions[yi]
                Vp_hor_gids = tp.GetElement(Vp_horizontal, [x,y])
                Vp_hor_status = nest.GetStatus(Vp_hor_gids)
                for (idx, n) in enumerate(Vp_hor_status):
                    if n['Tau_theta'] == 2.0:
                        print idx
                        try:
                            f[idx].add_subplot(len(positions), len(positions), yi*(len(positions))+xi+1)
                            tp.PlotTargets([Vp_hor_gids[idx]], Vp_horizontal, 'L4_exc', 'AMPA_syn', f[idx])
                        except:
                            print('%i bad' % Vp_hor_gids[idx])
        plt.show()

    # Connections from TpRelay to L4pyr in Vp (horizontally tuned)
    #topo.PlotTargets(topo.FindCenterElement(Tp), Vp_h, 'L4pyr', 'AMPA')
    #pylab.title('Connections TpRelay -> Vp(h) L4pyr')
    #pylab.show()

    # Connections from TpRelay to L4pyr in Vp (vertically tuned)
    #topo.PlotTargets(topo.FindCenterElement(Tp), Vp_v, 'L4pyr', 'AMPA')
    #pylab.title('Connections TpRelay -> Vp(v) L4pyr')
    #pylab.show()

    '''
    # pablo
    # Create vertical grating
    for n in nest.GetLeaves(Retina_layer)[0]:
            retina_0 = (nest.GetLeaves(Retina_layer)[0])[0]
            col = (n-retina_0)/Params['Np']

            cells_per_degree = Params['Np']/Params['visSize']
            cells_per_cycle = cells_per_degree/Params['spatial_frequency']

            nest.SetStatus([n], { "phase": col * 360/(cells_per_cycle-1) })
    '''
    ### keiko
    if Params['lambda_dg'] >= 0:
        [nest.SetStatus([n], {"phase": phaseInit(tp.GetPosition([n])[0],
                                                 Params["lambda_dg"],
                                                 Params["phi_dg"])})
        for n in nest.GetLeaves(Retina_layer)[0]]
    else:
        # Leonardo: Random retina input
        [nest.SetStatus([n], {"phase": phaseInit(tp.GetPosition([n])[0],
                                                 np.pi * np.random.rand(),
                                                 np.pi * np.random.rand())})
         for n in nest.GetLeaves(Retina_layer)[0]]


    # --------------------------------------------------------------------#
    # ---------- SET IB NEURONS ----------------------------------------- #
    # --------------------------------------------------------------------#

    # 30% of Cortex L56 excitatory neurons are intrinsically bursting(IB) neuron.
    # That is achieved by setting pacemaker current I_h.
    # So select 30% of L56_exc neuron, and change h_g_peak from 0.0 to 1.0.
    # (Other cortical neuron do not have I_h, thus h_g_peak=0.0)

    L56_vertical_idx = [nd for nd in nest.GetLeaves(Vp_vertical)[0] if nest.GetStatus([nd], 'model')[0]=='L56_exc']
    L56_horizontal_idx = [nd for nd in nest.GetLeaves(Vp_horizontal)[0] if nest.GetStatus([nd], 'model')[0]=='L56_exc']

    num_neuron = len(L56_vertical_idx)
    num_ib = int(num_neuron*0.3)

    ridx_vertical = np.random.randint(num_neuron, size=(1,num_ib))[0]
    ridx_horizontal = np.random.randint(num_neuron, size=(1,num_ib))[0]

    for i in range(1,num_ib,1):
        nest.SetStatus([L56_vertical_idx[ridx_vertical[i]]], {'h_g_peak': 1.0})
        nest.SetStatus([L56_horizontal_idx[ridx_horizontal[i]]], {'h_g_peak': 1.0})



    # initiate network activity
    #nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'rate': Params['ret_rate']})
    nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'rate': Params['ret_rate']})
    nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': 0.0})
    nest.Simulate(500.0)


    #! =================
    #! Recording devices
    #! =================

    nest.CopyModel('multimeter', 'RecordingNode',
            params = {'interval'   : Params['resolution'],
            #'record_from': ['V_m'],
            'record_from': ['V_m',
                            'I_syn_AMPA',
                            'I_syn_NMDA',
                            'I_syn_GABA_A',
                            'I_syn_GABA_B',
                            'g_AMPA',
                            'g_NMDA',
                            'g_GABAA',
                            'g_GABAB',
                            'I_NaP',
                            'I_KNa',
                            'I_T',
                            'I_h'],
            'record_to'  : ['memory'],
            'withgid'    : True,
            'withtime'   : True})

    recorders = []
    '''
    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Tp_layer  , 'Tp_inh'),
                              (Rp_layer  , 'Rp'),
                              (Vp_vertical, 'L23_exc'),
                              (Vp_horizontal, 'L23_exc'),
                              (Vp_vertical, 'L23_inh'),
                              (Vp_horizontal, 'L23_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc'),
                              (Vp_vertical, 'L4_inh'),
                              (Vp_horizontal, 'L4_inh'),
                              (Vp_vertical, 'L56_exc'),
                              (Vp_horizontal, 'L56_exc'),
                              (Vp_vertical, 'L56_inh'),
                              (Vp_horizontal, 'L56_inh'),
                              (Vs_vertical, 'L23_exc'),
                              (Vs_horizontal, 'L23_exc'),
                              (Vs_vertical, 'L23_inh'),
                              (Vs_horizontal, 'L23_inh'),
                              (Vs_vertical, 'L4_exc'),
                              (Vs_horizontal, 'L4_exc'),
                              (Vs_vertical, 'L4_inh'),
                              (Vs_horizontal, 'L4_inh'),
                              (Vs_vertical, 'L56_exc'),
                              (Vs_horizontal, 'L56_exc'),
                              (Vs_vertical, 'L56_inh'),
                              (Vs_horizontal, 'L56_inh')]:
    '''
    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Tp_layer  , 'Tp_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_vertical, 'L4_inh'),
                              (Vp_horizontal, 'L4_exc'),
                              (Vp_vertical, 'L23_exc'),
                              (Vp_horizontal, 'L23_exc'),
                              (Vp_vertical, 'L56_exc'),
                              (Rp_layer, 'Rp')]:
        rec = nest.Create('RecordingNode')
        recorders.append([rec,population,model])
        if (model=='Retina'):
            nest.SetStatus(rec,{'record_from': ['rate']})
        tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
        nest.Connect(rec, tgts)

    #! =================
    #! Spike detector
    #! =================
    detectors = []
    '''
    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Tp_layer  , 'Tp_inh'),
                              (Rp_layer  , 'Rp'),
                              (Vp_vertical, 'L23_exc'),
                              (Vp_horizontal, 'L23_exc'),
                              (Vp_vertical, 'L23_inh'),
                              (Vp_horizontal, 'L23_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc'),
                              (Vp_vertical, 'L4_inh'),
                              (Vp_horizontal, 'L4_inh'),
                              (Vp_vertical, 'L56_exc'),
                              (Vp_horizontal, 'L56_exc'),
                              (Vp_vertical, 'L56_inh'),
                              (Vp_horizontal, 'L56_inh'),
                              (Vs_vertical, 'L23_exc'),
                              (Vs_horizontal, 'L23_exc'),
                              (Vs_vertical, 'L23_inh'),
                              (Vs_horizontal, 'L23_inh'),
                              (Vs_vertical, 'L4_exc'),
                              (Vs_horizontal, 'L4_exc'),
                              (Vs_vertical, 'L4_inh'),
                              (Vs_horizontal, 'L4_inh'),
                              (Vs_vertical, 'L56_exc'),
                              (Vs_horizontal, 'L56_exc'),
                              (Vs_vertical, 'L56_inh'),
                              (Vs_horizontal, 'L56_inh')]:
        '''

    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Tp_layer  , 'Tp_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc')]:
        rec = nest.Create('spike_detector', params={"withgid": True, "withtime": True})
        #rec = nest.Create('spike_detector')
        detectors.append([rec,population,model])
        tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
        if model == 'Retina':
            for t in tgts:
                try:
                    nest.Connect([t], rec)
                    print('connected %d' % t)
                except:
                    print('%d did not work' % t)
        else:
            nest.Connect(tgts, rec)


    #! ====================
    #! Simulation
    #! ====================
    '''
    # change gKL to 0.8 in all populations (necessary to get a little stronger evoked response)
    for l in layers:
            sim_elements = l[1]['elements']
            for m in np.arange(0,np.size(sim_elements),1):

                    if(np.size(sim_elements)==1):
                        sim_model = sim_elements
                    else:
                        sim_model = sim_elements[m]

                    exec("la = %s" % l[0])
                    pop = [nd for nd in nest.GetLeaves(la)[0] if nest.GetStatus([nd], 'model')[0]==sim_model]
                    if (l[0]!='Retina_layer'):
                            for cell in pop:
                                    nest.SetStatus([cell], {'g_KL':0.8})
    '''
    # Simulate
    for t in Params['intervals']:

        #if (t == 250.0):  # Stimulus ON
        #    # if (t == 1500.0):  # Stimulus ON
        #    nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': 45.0})
        #else:  # Stimulus OFF
        #    nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': 0.0})

        if Params['input_flag']==True:
            nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': Params['ret_rate']})
        else:
            nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': 0.0})

        nest.Simulate(t)

    '''
    data_folder = Params['data_folder']
    if not os.path.isdir(data_folder):
        os.makedirs(data_folder)
    '''

    #! ====================
    #! Plot Results
    #! ====================

    if Params.has_key('show_main_figure') and Params['show_main_figure']:
        print "plotting..."

        rows = 9
        cols = 2

        fig = plt.figure(num=None, figsize=(13, 24), dpi=100)
        fig.subplots_adjust(hspace=0.4)

        # Plot A: membrane potential rasters

        recorded_models = [(Retina_layer,'Retina'),
                            (Vp_vertical,'L23_exc'),
                            (Vp_vertical,'L4_exc'),
                            (Vp_vertical,'L56_exc'),
                            (Rp_layer,'Rp'),
                            (Tp_layer,'Tp_exc')]

        #plotting.potential_raster(fig,recorders,recorded_models,0,Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,0)
        plotting.potential_raster(fig,recorders,recorded_models,0,Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,0)
        #starting_neuron = 800+1
        #plotting.potential_raster(fig,recorders,recorded_models,starting_neuron,Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,0)

        plt.title('Evoked')

        # Plot B: individual intracellular traces

        recorded_models =[(Vp_vertical,'L4_exc'),
                          (Vp_vertical,'L4_inh')]

        #plotting.intracellular_potentials(fig, recorders, recorded_models, 21, rows, cols, 6) #original
        # keiko
        total_time = 0.0
        for t in Params['intervals']:
            total_time += t

        #draw_neuron = (Params['Np']*Params['Np']/2)
        #plotting.intracellular_potentials(fig, recorders, recorded_models, draw_neuron, rows, cols, 6, total_time)
        plotting.intracellular_potentials(fig, recorders, recorded_models, 21, rows, cols, 6, total_time)
        #plotting.intracellular_potentials(fig, recorders, recorded_models, 820, rows, cols, 6, total_time)

        # Plot C: topographical activity of the vertical and horizontal layers


        if Params.has_key('start_membrane_potential') and  Params.has_key('end_membrane_potential'):
            start = Params['start_membrane_potential']
            stop = Params['end_membrane_potential']
        else:
            start = 130.0
            stop = 140.0

        recorded_models = [(Vp_vertical,'L23_exc')]
        labels = ["Vertical"]
        plotting.topographic_representation(fig,
                                            recorders,
                                            recorded_models,
                                            labels,
                                            Params['Np'],
                                            np.sum(Params['intervals']),
                                            Params['resolution'],
                                            rows,
                                            cols,
                                            start,
                                            stop,
                                            8,
                                            0)

        recorded_models = [(Vp_horizontal,'L23_exc')]

        labels = ["Horizontal"]

        plotting.topographic_representation(fig,recorders,recorded_models,labels,Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,start,stop,8,1)

        fig.savefig(data_folder + 'figure3.png', dpi=100)
        plt.show()


    # Plot D: movie

    #labels = ["Evoked_Vp_L23_Vertical","Evoked_Vp_L23_Horizontal"]
    #recorded_models = [(Vp_vertical,'L23_exc'),(Vp_horizontal,'L23_exc')]
    #plotting.makeMovie(fig,recorders,recorded_models,labels,Params['Np'],np.sum(Params['intervals']),Params['resolution'])


    #! ====================
    #! Save Results
    #! ====================

    print('save recorders data')

    for rec, population, model in recorders:

        # Get name of population
        for p in range(0, len(population_name), 1):
            if population_name[p]['population'] == population:
                p_name = population_name[p]['name']

        data = nest.GetStatus(rec)[0]['events']

        if model == 'Retina':
            scipy.io.savemat(data_folder + '/recorder_' + p_name + '_' + model + '.mat',
                             mdict={'senders': data['senders'],
                                    'rate': data['rate']})
        else:
            scipy.io.savemat(data_folder + '/recorder_' + p_name + '_' + model + '.mat',
                             mdict={'senders': data['senders'],
                                    'V_m': data['V_m'],
                                    'I_syn_AMPA': data['I_syn_AMPA'],
                                    'I_syn_NMDA': data['I_syn_NMDA'],
                                    'I_syn_GABA_A': data['I_syn_GABA_A'],
                                    'I_syn_GABA_B': data['I_syn_GABA_B'],
                                    'g_AMPA': data['g_AMPA'],
                                    'g_NMDA': data['g_NMDA'],
                                    'g_GABAA': data['g_GABAA'],
                                    'g_GABAB': data['g_GABAB']} )



    print('save raster images')
    plt.close()
    for rec, population, model in detectors:
        spikes = nest.GetStatus(rec, 'events')[0]

        # Get name of population
        for p in range(0, len(population_name), 1):
            if population_name[p]['population'] == population:
                p_name = population_name[p]['name']

        if len(nest.GetStatus(rec)[0]['events']['senders']) > 3:
            raster = raster_plot.from_device(rec, hist=True)
            pylab.title( p_name + '_' + model )
            f = raster[0].figure
            f.set_size_inches(15, 9)
            f.savefig(data_folder + 'spikes_' + p_name + '_' + model + '.png', dpi=100)
            plt.close()

            # Set filename and save spike data
            filename = data_folder + 'spike_' + p_name + '_' + model + '.pickle'
            pickle.dump(spikes, open(filename, 'w'))
            scipy.io.savemat(data_folder + '/spike_' + p_name + '_' + model + '.mat', mdict={'senders': spikes['senders'], 'times': spikes['times']})

            #filename_AMPA = data_folder + 'connection_' + p_name + '_AMPA_syn' + '.dat'
            #tp.DumpLayerConnections(population, 'AMPA_syn', filename_AMPA)

            '''
            filename_AMPA = data_folder + 'connection_' + p_name + '_AMPA_syn' + '.dat'
            filename_NMDA = data_folder + 'connection_' + p_name + '_NMDA_syn' + '.dat'
            filename_GABAA = data_folder + 'connection_' + p_name + '_GABA_A_syn' + '.dat'
            filename_GABAB = data_folder + 'connection_' + p_name + '_GABA_B_syn' + '.dat'
            tp.DumpLayerConnections(population, 'AMPA_syn', filename_AMPA)
            tp.DumpLayerConnections(population, 'NMDA_syn', filename_NMDA)
            tp.DumpLayerConnections(population, 'GABA_A_syn', filename_GABAA)
            tp.DumpLayerConnections(population, 'GABA_B_syn', filename_GABAB)
            '''
    '''
    for p in range(0, len(population_name), 1):

        population = population_name[p]['population']
        p_name = population_name[p]['name']
        filename_nodes = data_folder + '/gid_' + p_name + '.dat'

        tp.DumpLayerNodes(population, filename_nodes)
    '''

    network_script = Params['network'] + '.py'
    shutil.copy2(network_script, Params['data_folder'] + network_script)

    print('end')