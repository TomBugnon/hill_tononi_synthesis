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

import static_input_utils

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


    # TODO this should come from network file
    synapse_models = ['AMPA_syn', 'NMDA_syn', 'GABA_A_syn', 'GABA_B_syn', 'reticular_projection']
    # synapse_models = ['AMPA_syn', 'NMDA_syn', 'GABA_A_syn', 'GABA_B_syn']

    # Create models
    print('Creating models...')
    for m in models:
        # print('.', end="")
#        print(m), print('\n')
        nest.CopyModel(m[0], m[1], m[2])
    print(' done.')

    # Create layers, store layer info in Python variable
    print('Creating layers...')
    for l in layers:
        # print('.', end="")
        exec '%s = tp.CreateLayer(l[1])' % l[0] in globals(), locals()
    print(' done.')

    # Create connections, need to insert variable names
    for c in conns:
            eval('tp.ConnectLayers(%s,%s,c[2])' % (c[0], c[1]))


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







    # --------------------------------------------------------------------#
    # ---------- SET RETINA STATIC FIRING RATE -------------------------- #
    # --------------------------------------------------------------------#

    if not Params['new_image']:
        luminances = np.load(Params['image_folder']+Params['image_name'])
    else:
        luminances = static_input_utils.create_default_image(Params['new_image_elements'], Params['Np'])
        if not Params['image_name']=='':
            np.save(Params['image_folder']+Params['image_name'], luminances)


    # For each 'pixel' (retina neuron) we set the rate of the poisson generator to luminance * ret_rate + ret_rate_baseline

    luminances = static_input_util.get_luminances(Params['image_folder']+Params['image_name'])
    if not np.shape(luminances) == (Params['Np'], Params['Np']):
        error('Wrong definition for the input image %s' % (Params['image_folder']+Params['image_name']))


    for n in nest.GetLeaves(Retina_layer)[0]:
        nest.SetStatus([n], {"rate": luminances[tp.GetPosition(n)] * Params['ret_rate'] + Params['ret_rate_baseline']})


     # --------------------------------------------------------------------#
     # ---------------------------- INITIATE  ---------------------------- #
     # --------------------------------------------------------------------#
    nest.Simulate(Params['initiation_time'])



    #! =================
    #! Recording devices
    #! =================
    print('Connecting recorders...')
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
    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Rp_layer  , 'Rp'),
                              (Vp_vertical, 'L23_exc'),
                              (Vp_horizontal, 'L23_exc'),
                              (Vp_vertical, 'L23_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc'),
                              (Vp_vertical, 'L4_inh'),
                              (Vp_vertical, 'L56_exc'),
                              (Vp_horizontal, 'L56_exc'),
                              (Vs_vertical, 'L23_exc'),
                              (Vs_horizontal, 'L23_exc'),
                              (Vs_vertical, 'L4_exc'),
                              (Vs_horizontal, 'L4_exc'),
                              (Vs_vertical, 'L56_exc'),
                              (Vs_horizontal, 'L56_exc'),
                              (Vs_cross, 'L23_exc'),
                              (Vs_cross, 'L4_exc'),
                              (Vs_cross, 'L56_exc'),


                              ]:
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
    retina_failed_detectors = []
    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Tp_layer  , 'Tp_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc'),
                              (Vp_cross, 'L4_exc')]:

        rec = nest.Create('spike_detector', params={"withgid": True, "withtime": True})
        #rec = nest.Create('spike_detector')
        detectors.append([rec,population,model])
        tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
        if model == 'Retina':
            for t in tgts:
                try:
                    #nest.Connect(rec, [t])
                    nest.Connect([t], rec)
                    #print('connected %d' % t)
                except:
                    print('%d did not work' % t)
                    retina_failed_detectors.append(t)
        else:
            nest.Connect(tgts, rec)

    sg = [nd for nd in nest.GetLeaves(Retina_layer)[0] if nest.GetStatus([nd], 'model')[0] == 'Retina']






    #! ====================
    #! Simulation
    #! ====================
    nest.Simulate(Params['intervals'])




    #! ====================
    #! Plot
    #! ====================

    #Figure 3:

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







    # Plot C: topographical activity of all layers


    if Params.has_key('start_membrane_potential') and Params.has_key('end_membrane_potential'):
        start = Params['start_membrane_potential']
        stop = Params['end_membrane_potential']
    else:
        start = 130.0
        stop = 140.0

    vertical_models = [(Retina_layer, 'Retina'),
                       (Vp_vertical, 'L23_exc'),
                       (Vp_vertical, 'L4_exc'),
                       (Vp_vertical, 'L56_exc'),
                       (Vs_vertical, 'L23_exc'),
                       (Vs_vertical, 'L4_exc'),
                       (Vs_vertical, 'L56_exc')]

    horizontal_models = [(Retina_layer, 'Retina'),
                         (Vp_horizontal, 'L23_exc'),
                         (Vp_horizontal, 'L4_exc'),
                         (Vp_horizontal, 'L56_exc'),
                         (Vs_horizontal, 'L23_exc'),
                         (Vs_horizontal, 'L4_exc'),
                         (Vs_horizontal, 'L56_exc')]

    cross_models = [(Retina_layer, 'Retina'),
                    (),
                    (),
                    (),
                    (Vs_cross, 'L23_exc'),
                    (Vs_cross, 'L4_exc'),
                    (Vs_cross, 'L56_exc')]


    recorded_models = [vertical_models, horizontal_models, cross_models]
    labels = ["Vertical", 'Horizontal', 'Cross']
    plotting.all_topographic(recorders, recorded_models, labels, number_cells, simtime, resolution, start, stop)



    fig.savefig(data_folder + '/figure_all_topo.png', dpi=100)

    if Params.has_key('show_main_figure') and Params['show_main_figure']:
        plt.show()

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
    data_folder = Params['data_folder']

    if not os.path.isdir(data_folder):
        os.makedirs(data_folder)

    # To save spike data, set pairs of population id and its name
    population_name = [ {'population': Retina_layer, 'name': 'Retina'},
                        {'population': Vp_vertical, 'name': 'Vp_v'},
                        {'population': Vp_horizontal, 'name': 'Vp_h'},
                        {'population': Rp_layer, 'name': 'Rp'},
                        {'population': Tp_layer, 'name': 'Tp'},
                        {'population': Vs_vertical, 'name': 'Vs_v'},
                        {'population': Vs_horizontal, 'name': 'Vs_h'},
                        {'population': Vs_cross, 'name': 'Vs_c'}]

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
                                    'V_m': data['V_m']})#,
                                    # 'I_syn_AMPA': data['I_syn_AMPA'],
                                    # 'I_syn_NMDA': data['I_syn_NMDA'],
                                    # 'I_syn_GABA_A': data['I_syn_GABA_A'],
                                    # 'I_syn_GABA_B': data['I_syn_GABA_B'],
                                    # 'g_AMPA': data['g_AMPA'],
                                    # 'g_NMDA': data['g_NMDA'],
                                    # 'g_GABAA': data['g_GABAA'],
                                    # 'g_GABAB': data['g_GABAB']} )



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

    shutil.copy2('network_full_tom.py', Params['data_folder'] + 'network_full_keiko.py')
    shutil.copy2('figure_3_plot.py', Params['data_folder'] + 'figure_3_plot.py')
    shutil.copy2('main.py', Params['data_folder'] + 'main.py')
    #scipy.io.savemat(data_folder + '/retina_failed_detectors.mat', mdict={'detectors': retina_failed_detectors})
    scipy.io.savemat(data_folder + '/retina_spike_times.mat', mdict={'spike_times': retina_spikes})

    print('end')