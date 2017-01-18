#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Simulation of figure 3 of Hill-Tononi paper
# Author: Pablo Martinez Cañada (pablomc@ugr.es)
#
# Evoked activity: the stimulus is ON during the second interval.
# I used the same parameters of Figure 4 but I had to change gKL to 0.8
# to get a little stronger evoked response.

from __future__ import print_function
import nest
import nest.topology as tp
import numpy as np
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

def file_len(fname):
    i = -1
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def phaseInit(pos, lam, alpha):
    '''Initializer function for phase of drifting grating nodes.

       pos  : position (x,y) of node, in degree
       lam  : wavelength of grating, in degree
       alpha: angle of grating in radian, zero is horizontal

       Returns number to be used as phase of sinusoidal Poisson generator.
    '''
    return 360.0 / lam * (math.cos(alpha) * pos[0] + math.sin(alpha) * pos[1])


def simulation(Params):

    root_folder = Params['root_folder']
    if not os.path.isdir(root_folder):
        os.makedirs(root_folder)

    net_folder = root_folder + Params['net_folder']
    if not os.path.isdir(net_folder):
        os.makedirs(net_folder)

    conn_folder = net_folder + '/connections'
    if not os.path.isdir(conn_folder):
        os.makedirs(conn_folder)

    data_folder = net_folder + Params['data_folder']
    if not os.path.isdir(data_folder):
        os.makedirs(data_folder)

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

    # TODO this should come from network file
    synapse_models = ['AMPA_syn', 'NMDA_syn', 'GABA_A_syn', 'GABA_B_syn', 'reticular_projection']
    # synapse_models = ['AMPA_syn', 'NMDA_syn', 'GABA_A_syn', 'GABA_B_syn']

    # Create models
    print('Creating models...', end="")
    for m in models:
        print('.', end="")
        nest.CopyModel(m[0], m[1], m[2])
    print(' done.')

    # Create layers, store layer info in Python variable
    print('Creating layers...', end="")
    for l in layers:
        print('.', end="")
        exec '%s = tp.CreateLayer(l[1])' % l[0] in globals(), locals()
    print(' done.')

    # now either create brand new connections, or load connectivity files from
    # connections_file
    connections_file =  conn_folder + '/connections.pickle'

    # assume alread runned before, and it is going to load connections and layers
    if Params.has_key('load_connections') and Params['load_connections']:

        if not os.path.isfile(connections_file):
            raise IOError('File with connections was never created: %s' % connections_file)

        print('Loading connections...', end="")
        with open(connections_file) as f:
            (pre, post_intact, post_scrambled, w, d, syn) = pickle.load(f)
        print(' done.')

        if Params['scrambled']:
            post = post_scrambled
        else:
            post = post_intact

        print('Connecting neurons...', end="")
        con_dict = {'rule': 'one_to_one'}
        for (pi, ti, wi, di, syni) in zip(pre, post, w, d, syn):
            print('.', end="")
            # print(syni)

            syn_dict = {"model": syni,
                        'weight': wi,
                        'delay': di
                        }
            nest.Connect(pi, ti, con_dict, syn_dict)
        print(' done.')

        # print('Connection retina...', end="")
        # # Create connections, need to insert variable names
        # for c in conns:
        #     if c[0]=='Retina_layer':
        #         eval('tp.ConnectLayers(%s,%s,c[2])' % (c[0], c[1]))
        # print('done.')

    else:
        # Create connections, need to insert variable names
        print('Creating connectiions...', end="")
        for c in conns:
            if Params['Np'] < 40:
                print('%s -> %s' % (c[0], c[1]))
                c[2]['allow_oversized_mask'] = True
            else:
                print('.', end="")
            eval('tp.ConnectLayers(%s,%s,c[2])' % (c[0], c[1]))
        print(' done.')

    if Params.has_key('dump_connections') and Params['dump_connections']:

        import glob

        # First dump all nodes and connections to files
        for l in layers:
            print('dumping connections for %s' % l[0])
            eval("tp.DumpLayerNodes(%s, conn_folder + '/nodes_%s.dump')" % (l[0], l[0]) )
            for synapse in synapse_models:
                 eval("tp.DumpLayerConnections(%s, '%s', conn_folder + '/connections_%s_%s.dump')" % (l[0], synapse, l[0], synapse) )

        # Now load all nodes, and mark which physical layer they belong to
        all_nodes_files = glob.glob(conn_folder + "/nodes_*.dump")
        if len(all_nodes_files) < 1:
            raise ValueError('No files save for this configuration. First '
                             'run with load_connection = False')

        pop_nodes = [None] * len(all_nodes_files)
        pop_names = [None] * len(all_nodes_files)
        perm_pop_nodes = [None] * len(all_nodes_files)
        for (idx, pop_file) in enumerate(all_nodes_files):

            population = pop_file[pop_file.find('nodes_')+6:-5]
            pop_names[idx] = population
            print('Loading nodes for population %s' % population)
            pop_nodes[idx] = np.loadtxt(pop_file)

            # identify physical layers
            n_layer = 0
            n_nodes = pop_nodes[idx].shape[0]
            layers = np.zeros((n_nodes,1))
            for i_line in range(n_nodes):
                if pop_nodes[idx][i_line, 1] == -3.9 and pop_nodes[idx][i_line, 2] == 3.9:
                    n_layer += 1
                layers[i_line] = n_layer
            pop_nodes[idx] = np.append(pop_nodes[idx], layers, 1)

        # Finally, load connections and save a big vector for each one of the 5
        # parameters needed by nest.Connect: pre, post, w, d, and synapse_model
        # Save unmodified and scrambled versions in network.pickle

        all_conn_files = glob.glob(conn_folder + "/connections_*.dump")

        pre = []
        post_intact = []
        post_scrambled = []
        w = []
        d = []
        syn = []
        for (idx, pop_conn_file) in enumerate(all_conn_files):

            population = pop_conn_file[pop_conn_file.find('connections_')+12:-5]
            print('Scrambling connections for %s (%d out of %d files)' % (population, idx+1, len(all_conn_files)))

            synapse_model = [s for i,s in enumerate(synapse_models) if s in population][0]

            print('Loading file')
            pop_conns = np.loadtxt(pop_conn_file)
            tn_lines = pop_conns.shape[0]

            if tn_lines > 0:

                syn.append(synapse_model)

                pre.append(pop_conns[:,0].astype(int))
                post_intact.append(pop_conns[:,1].astype(int))
                w.append(pop_conns[:,2])
                d.append(pop_conns[:,3].astype(int))

                # Do the actual scrfambling
                # TODO this should be optimized, it is way too slow
                if any(w == population[0:2] for w in ['Vs', 'Vp']):
                    print('Scrambling...', end="")
                    for i_line in range(tn_lines):
                        if not i_line % round(tn_lines/10.):
                            print('%d%%, ' % (round(100.*i_line/tn_lines) if i_line > 0 else 0), end="")

                        for (i_pop, this_pop) in enumerate(pop_nodes):
                            t_idx = np.where(this_pop[:,0] == pop_conns[i_line, 1])[0]
                            if t_idx.size:
                                break
                        if not t_idx.size:
                            raise ValueError('Could not find tartget %s for node %s' % (pop_conns[i_line, 1], pop_conns[i_line, 0]))
                        if len(t_idx) > 1:
                            raise ValueError('Target %s apears in more than one population' % (pop_conns[i_line, 1], pop_conns[i_line, 0]))
                        new_t_idx = t_idx[0]

                        # For the moment, only scramble connections within the visual cortex
                        if any(w == pop_names[i_pop][0:2] for w in ['Vs', 'Vp']):
                            possible_targets = pop_nodes[i_pop][pop_nodes[i_pop][:,3] == pop_nodes[i_pop][new_t_idx,3], 0]
                            pop_conns[i_line, 1] = np.random.permutation(possible_targets)[0]

                print('done.')

                post_scrambled.append(pop_conns[:,1].astype(int))

        # save results
        print('Saving results...', end="")
        with open(connections_file, 'w') as f:
            pickle.dump((pre, post_intact, post_scrambled, w, d, syn), f)
        print(' done.')

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


    if Params.has_key('show_center_connectivity_figure') and Params['show_center_connectivity_figure']:

        # this code only gives one element, and you have no control over the
        # model it comes from...
        # tp.PlotTargets(tp.FindCenterElement(Vp_horizontal), Vp_horizontal, 'L56_exc', 'AMPA_syn')
        # tp.PlotTargets(tp.FindCenterElement(Vp_horizontal), Vp_vertical, 'L56_exc', 'AMPA_syn')
        # tp.PlotTargets(tp.FindCenterElement(Vp_vertical), Vp_horizontal, 'L56_exc', 'AMPA_syn')
        # tp.PlotTargets(tp.FindCenterElement(Vp_vertical), Vp_vertical, 'L56_exc', 'AMPA_syn')

        # this way I can get all nodes in the positions, and later filter per Model using GetStatus
        # IMPORTANT: when layer has even number of neurons per line/column
        #            using center (0,0) makes the function return 4 neurons
        #            per postiion since they all have the same distance.
        #            Using (-0,1, 0,1) solves this problem.
        #            When using odd number perhaps this has to change...

        for src_label in ['Tp_layer', 'Vp_vertical', 'Vp_horizontal']:
            for tgt_label in ['Tp_layer', 'Vp_vertical', 'Vp_horizontal']:

                if src_label == 'Tp_layer' and src_label == tgt_label:
                    continue

                print('source population %s' % src_label)
                print('target population %s' % tgt_label)
                n_plot = 0
                f = plt.figure()
                f.canvas.set_window_title('AMPA Connections: %s -> %s' % (src_label, tgt_label))
                all_center = eval('tp.FindNearestElement(%s, (-0.1, 0.1), True)[0]' % src_label)
                for n in all_center:
                    s = nest.GetStatus([n])
                    p = tp.GetPosition([n])
                    # print('%s : (%2.2f, %2.2f)' % (s[0]['model'], p[0][0], p[0][1]))
                    m1 = s[0]['model'].name
                    print('Source neuron model %s' % m1)
                    if m1.find('_exc') > -1:
                        if src_label.find('Tp_') > -1:
                            print( 'found Tp_')
                            # target has to be one of Vp_ or Vs_ models
                            target_models = ['L23_exc', 'L4_exc', 'L56_exc']
                        elif src_label.find('Vp_') > -1 or src_label.find('Vs') > -1:
                            # if target is Vp_ of Vs_ too, then one of those models
                            if tgt_label.find('Vp_') > -1:
                                target_models = ['L23_exc', 'L4_exc', 'L56_exc']
                            elif tgt_label.find('Tp_') > -1:
                                # otherwise it has to be a Tp_target
                                target_models = ['Tp_exc']
                            else:
                                raise ValueError('Invalide target %s for source model: %s' % (tgt_label, src_label))
                        else:
                            raise ValueError('Invalide source model: %s' % (src_label))

                        for m2 in target_models:
                            print('Target neuron model %s' % m2)
                            try:
                                get_targets_command = 'tp.GetTargetNodes([n], %s, m2, "AMPA_syn")[0]' % tgt_label
                                print(get_targets_command)
                                targets = eval(get_targets_command)
                                if len(targets) > 0:
                                    n_plot += 1
                                    f.add_subplot(3,5,n_plot)
                                    eval('tp.PlotTargets([n], %s, m2, "AMPA_syn", f)' % tgt_label)
                                    plt.title('%s -> %s' % (m1, m2), fontsize=9)
                            except:
                                print('didnt work')

                f.savefig(data_folder + '/connectivity_%s_%s.png' % (src_label, tgt_label), dpi=100)

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
                        print(idx)
                        try:
                            f[idx].add_subplot(len(positions), len(positions), yi*(len(positions))+xi+1)
                            tp.PlotTargets([Vp_hor_gids[idx]], Vp_horizontal, 'L4_exc', 'AMPA_syn', f[idx])
                        except:
                            print('%i bad' % Vp_hor_gids[idx])
        plt.show()

    if Params.has_key('dry_run') and Params['dry_run']:
        print('Only generation, loading and ploting network. No actual simulation done.')
        return

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
    print('Connecting recorders...', end="")
    nest.CopyModel('multimeter', 'RecordingNode',
            params = {'interval'   : Params['resolution'],
            'record_from': ['V_m'
                            # Put back when plotting synaptic currents, otherwise makes everything slower for no reason
                            # 'I_syn_AMPA',
                            # 'I_syn_NMDA',
                            # 'I_syn_GABA_A',
                            # 'I_syn_GABA_B',
                            #'g_AMPA',
                            #'g_NMDA',
                            #'g_GABAA',
                            #'g_GABAB',
                            #'I_NaP',
                            #'I_KNa',
                            #'I_T',
                            #'I_h'
                            ],
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
        print('.', end="")
        rec = nest.Create('RecordingNode')
        recorders.append([rec,population,model])
        if (model=='Retina'):
            nest.SetStatus(rec,{'record_from': ['rate']})
        tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
        #nest.Connect(rec, tgts, None, {'delay': recorders_delay})
        nest.Connect(rec, tgts)
    print('done.')

    #! =================
    #! Spike detector
    #! =================
    print('Connecting detectors...', end="")
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
                              (Vp_vertical, 'L23_exc'),
                              (Vp_horizontal, 'L23_exc'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc')]:
        print('.', end="")
        rec = nest.Create('spike_detector', params={"withgid": True, "withtime": True})
        #rec = nest.Create('spike_detector')
        detectors.append([rec,population,model])
        tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
        if model == 'Retina':
            for t in tgts:
                try:
                    # nest.Connect([t], rec, None, {'delay':  recorders_delay})
                    nest.Connect([t], rec)
                    print('connected %d' % t)
                except:
                    print('%d did not work' % t)
        else:
            # nest.Connect(tgts, rec, None, {'delay': recorders_delay})
            nest.Connect(tgts, rec)
    print('done.')


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



    #! ====================
    #! Plot Results
    #! ====================

    print("Creating figure 3...")

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
    plotting.potential_raster(fig,recorders,recorded_models,0,Params['Np'],
                              np.sum(Params['intervals']),
                              Params['resolution'],
                              rows,cols,0)
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

    fig.savefig(data_folder + '/figure3.png', dpi=100)

    if Params.has_key('show_main_figure') and Params['show_main_figure']:
        plt.show()


    # Plot D: movie

    #labels = ["Evoked_Vp_L23_Vertical","Evoked_Vp_L23_Horizontal"]
    #recorded_models = [(Vp_vertical,'L23_exc'),(Vp_horizontal,'L23_exc')]
    #plotting.makeMovie(fig,recorders,recorded_models,labels,Params['Np'],np.sum(Params['intervals']),Params['resolution'])


    #! ====================
    #! Save Results
    #! ====================

    print('save recorders data')
    # Set folder
    #rootdir = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/'
    #expdir = 'random/'
    # expdir = 'random_full/'
    # expdir = 'structured_full/'
    # data_folder = rootdir + expdir

    # To save spike data, set pairs of population id and its name
    population_name = [ {'population': Retina_layer, 'name': 'Retina'},
                        {'population': Vp_vertical, 'name': 'Vp_v'},
                        {'population': Vp_horizontal, 'name': 'Vp_h'},
                        {'population': Rp_layer, 'name': 'Rp'},
                        {'population': Tp_layer, 'name': 'Tp'},
                        {'population': Vs_vertical, 'name': 'Vs_v'},
                        {'population': Vs_horizontal, 'name': 'Vs_h'}]

    if Params.has_key('save_recorders') and Params['save_recorders']:
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
                                        'V_m': data['V_m']
                                        #'I_syn_AMPA': data['I_syn_AMPA'],
                                        #'I_syn_NMDA': data['I_syn_NMDA'],
                                        #'I_syn_GABA_A': data['I_syn_GABA_A'],
                                        #'I_syn_GABA_B': data['I_syn_GABA_B'],
                                        #'g_AMPA': data['g_AMPA'],
                                        #'g_NMDA': data['g_NMDA'],
                                        #'g_GABAA': data['g_GABAA'],
                                        #'g_GABAB': data['g_GABAB']
                                        } )



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
            f.savefig(data_folder + '/spikes_' + p_name + '_' + model + '.png', dpi=100)
            plt.close()

            # Set filename and save spike data
            filename = data_folder + '/spike_' + p_name + '_' + model + '.pickle'
            pickle.dump(spikes, open(filename, 'w'))
            scipy.io.savemat(data_folder + '/spike_' + p_name + '_' + model + '.mat', mdict={'senders': spikes['senders'], 'times': spikes['times']})

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
    shutil.copy2(network_script, data_folder + '/' + network_script)
    shutil.copy2(network_script, conn_folder + '/' + network_script)

    print('end')