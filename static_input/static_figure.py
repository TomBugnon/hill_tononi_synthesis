#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Simulation of figure 3 of Hill-Tononi paper
# Author: Pablo Martinez Cañada (pablomc@ugr.es)
#
# Evoked activity: the stimulus is ON during the second interval.
# I used the same parameters of Figure 4 but I had to change gKL to 0.8
# to get a little stronger evoked response.

from __future__ import print_function #Because pythoon 2.7 doesn't support print(..., end = ..)


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
reload(static_input_utils)


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
    print('Creating models...', end="")
    for m in models:
        print('.', end="")
#        print(m), print('\n')
        nest.CopyModel(m[0], m[1], m[2])
    print(' done.')

    # Create layers, store layer info in Python variable
    print('Creating layers...', end="")
    for l in layers:
        print('.', end="")
        exec '%s = tp.CreateLayer(l[1])' % l[0] in globals(), locals()
    print(' done.')

    # Create connections, need to insert variable names
    print('Creating connections...', end="")
    for c in conns:
        print('.', end="")
        eval('tp.ConnectLayers(%s,%s,c[2])' % (c[0], c[1]))
    print(' done.')


    # --------------------------------------------------------------------#
    # ---------- SET IB NEURONS ----------------------------------------- #
    # --------------------------------------------------------------------#

    print('Set IB neurons')
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
        nest.SetStatus([L56_vertical_idx[ridx_vertical[i]]], {'g_peak_h': 1.0})
        nest.SetStatus([L56_horizontal_idx[ridx_horizontal[i]]], {'g_peak_h': 1.0})







    # --------------------------------------------------------------------#
    # ---------- SET RETINA STATIC FIRING RATE -------------------------- #
    # --------------------------------------------------------------------#
    print('Set Retina inputs neurons')
    if not Params['new_image']: #Read
        with open(Params['image_folder']+Params['image_name']+'.pickle', 'rb') as handle:
            image_dic = pickle.load(handle)
    else:
        image_dic = static_input_utils.create_default_image(Params['new_image_elements'], Params['Np'])
        if not Params['image_name']=='':
            static_input_utils.save_new_image_dic(Params['image_folder'] + Params['image_name'], image_dic)

    luminances = image_dic['luminances']
    fig = plt.figure()
    plt.imshow(luminances, interpolation='nearest')
    fig.savefig(data_folder + '/input_image.png', dpi=100)


    if not np.shape(luminances) == (Params['Np'], Params['Np']):
        raise ValueError('Wrong definition for the input image %s' % (Params['image_folder']+Params['image_name']))

    # For each 'pixel' (retina neuron) we set the rate of the poisson generator to maxrate * ret_rate + minrate
    #if 'automatic_min_rate', minrate is determined so that mean retinal activity is mean_rate
    #Iterate over indices of image.
    if Params['automatic_min_rate']:
        N_on = np.sum(luminances)
        N_tot = (Params['Np'] ** 2)
        minrate = (N_tot * Params['mean_rate'] - N_on * Params['ret_rate'] )/ ( N_tot - N_on)
    else:
        minrate = Params['ret_rate_baseline']

    for i in range(Params['Np']):
        for j in range(Params['Np']):
            #Get node, set corresponding luminance
            currnode = tp.GetElement(Retina_layer, (j, i))[0] #NB: Because of the inconsistency of the GetElement function we have to transpose.
            if luminances[i,j] > 0.:
                nest.SetStatus([currnode], {"rate": Params['ret_rate']})
            else:
                nest.SetStatus([currnode], {"rate": minrate})



     # --------------------------------------------------------------------#
     # ---------------------------- INITIATE  ---------------------------- #
     # --------------------------------------------------------------------#
    print('Simulate %i secs' % Params['initiation_time'])
    nest.Simulate(Params['initiation_time'])



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
    print('done.')

    #! =================
    #! Spike detector
    #! =================
    print('Connecting detectors...', end="")
    detectors = []
    retina_failed_detectors = []
    for population, model in [(Retina_layer, 'Retina'),
                              (Tp_layer  , 'Tp_exc'),
                              (Tp_layer  , 'Tp_inh'),
                              (Vp_vertical, 'L4_exc'),
                              (Vp_horizontal, 'L4_exc'),
                              (Vs_cross, 'L4_exc')]:

        rec = nest.Create('spike_detector', params={"withgid": True, "withtime": True})
        #rec = nest.Create('spike_detector')
        detectors.append([rec,population,model])
        tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
        if model == 'Retina':
            for t in tgts:
                try:
                    # nest.Connect([t], rec, None, {'delay':  recorders_delay})
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
    # Simulate
    for t in Params['intervals']:

        # TODO: useless as is
        # if Params['input_flag']==True:
        #     nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': Params['ret_rate']})
        # else:
        #     nest.SetStatus(nest.GetLeaves(Retina_layer)[0], {'amplitude': 0.0})

        nest.Simulate(t)



    #! ====================
    #! Plot Results
    #! ====================

    print("Creating figure 3...")

    rows = 9
    cols = 2

    #Tom: We plot temporal data for the column to which the vertical line in the image is presented
    #nest.GetLeaves() numbers the elements starting top left and descending vertically.
    # Therefore starting_neuron (given to plotting.potential_rater) must be the index of the node in the first row.
    # starting_neuron = 0 plots the first column, starting_neuron = 40 the second, etc
    vertical_bar_column = image_dic['vertical'][1]
    starting_neuron_Vp = Params['Np']*vertical_bar_column




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
    plotting.potential_raster(fig,recorders,recorded_models,starting_neuron_Vp,Params['Np'],
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

    if Params.has_key('plot_topo_all_regions') and Params['plot_topo_all_regions']:

        if Params.has_key('start_membrane_potential') and Params.has_key('end_membrane_potential'):
            start = Params['start_membrane_potential']
            stop = Params['end_membrane_potential']
        else:
            start = 130.0
            stop = 140.0

        ####### USER INPUT

        # First column
        vertical_models = [(Retina_layer, 'Retina'),
                           (Vp_vertical, 'L23_exc'),
                           (Vp_vertical, 'L4_exc'),
                           (Vp_vertical, 'L56_exc'),
                           (Vs_vertical, 'L23_exc'),
                           (Vs_vertical, 'L4_exc'),
                           (Vs_vertical, 'L56_exc')]
        # Second column
        horizontal_models = [(Retina_layer, 'Retina'),
                             (Vp_horizontal, 'L23_exc'),
                             (Vp_horizontal, 'L4_exc'),
                             (Vp_horizontal, 'L56_exc'),
                             (Vs_horizontal, 'L23_exc'),
                             (Vs_horizontal, 'L4_exc'),
                             (Vs_horizontal, 'L56_exc')]
        # Third column
        cross_models = [(Retina_layer, 'Retina'),
                        (),
                        (),
                        (),
                        (Vs_cross, 'L23_exc'),
                        (Vs_cross, 'L4_exc'),
                        (Vs_cross, 'L56_exc')]

        # Column labels
        labels = ['Vertical', 'Horizontal', 'Cross']

        # Row labels
        areas = ['Retina', 'Primary\nL2/3', 'Primary\nL4', 'Primary\nL5/6', 'Secondary\nL2/3', 'Secondary\nL4',
                 'Secondary\nL4/6']

        # Number of neurons in each plotted layer:
        Np, Ns = Params['Np'], Params['Ns']
        n_neurons = [Np, Np, Np, Np, Ns, Ns, Ns]

        recorded_models = [vertical_models, horizontal_models, cross_models]

        ######## end USER INPUT

        fig = plotting.all_topographic(recorders, recorded_models, labels, areas,
                                            n_neurons,
                                            np.sum(Params['intervals']),
                                            Params['resolution'], start, stop)



        fig.savefig(data_folder + '/figure_all_topo.png', dpi=100)

        if Params.has_key('show_main_figure') and Params['show_main_figure']:
            plt.show()




    # Plot F: All areas


    if Params.has_key('plot_all_regions') and Params['plot_all_regions']:

        print("Creating plot of all areas")

        vertical_bar_column_Vp = image_dic['vertical'][1]
        vertical_bar_column_Vs = vertical_bar_column_Vp * Params['Ns']/Params['Np']
        starting_neuron_Vp = Params['Np'] * vertical_bar_column_Vp
        starting_neuron_Vs = Params['Ns'] * vertical_bar_column_Vs

        #### USER INPUT

        vertical_models = [(Retina_layer,'Retina'),
                        (Vp_vertical,'L23_exc'),
                        (Vp_vertical,'L4_exc'),
                        (Vp_vertical,'L56_exc'),
                        (Vs_vertical,'L23_exc'),
                        (Vs_vertical, 'L4_exc'),
                        (Vs_vertical, 'L56_exc')]

        horizontal_models = [(Retina_layer,'Retina'),
                        (Vp_horizontal,'L23_exc'),
                        (Vp_horizontal,'L4_exc'),
                        (Vp_horizontal,'L56_exc'),
                        (Vs_horizontal,'L23_exc'),
                        (Vs_horizontal, 'L4_exc'),
                        (Vs_horizontal, 'L56_exc')]

        cross_models = [(Retina_layer,'Retina'),
                        (),
                        (),
                        (),
                        (Vs_cross,'L23_exc'),
                        (Vs_cross, 'L4_exc'),
                        (Vs_cross, 'L56_exc')]

        labels = ['Vertical', 'Horizontal', 'Cross']

        # Row labels
        areas = ['Retina', 'Primary\nL2/3', 'Primary\nL4', 'Primary\nL5/6', 'Secondary\nL2/3', 'Secondary\nL4',
                 'Secondary\nL4/6']

        #TODO: make this look nicer
        # Number of neurons in each plotted layer:
        Np, Ns = Params['Np'], Params['Ns']
        n_neurons = [Np, Np, Np, Np, Ns, Ns, Ns]

        # Starting neuron in each plotted layer:
        starting_neurons = [starting_neuron_Vp, starting_neuron_Vp, starting_neuron_Vp, starting_neuron_Vp, starting_neuron_Vs, starting_neuron_Vs, starting_neuron_Vs]

        plotcols = [vertical_models, horizontal_models, cross_models]

        #### end USER INPUT


        fig = plotting.potential_raster_multiple_models(fig, recorders, plotcols, labels, areas, starting_neurons, n_neurons,
                                                  np.sum(Params['intervals']),
                                                  Params['resolution'],
                                                  0)

        fig.savefig(data_folder + '/figure_all_areas.png', dpi=100)



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
                        {'population': Vs_horizontal, 'name': 'Vs_h'},
                        {'population': Vs_cross, 'name': 'Vs_c'}]

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
                                        'times': data['times'],
                                        'V_m': data['V_m']
                                        #'I_syn_AMPA': data['I_syn_AMPA'],
                                        #'I_syn_NMDA': data['I_syn_NMDA'],
                                        #'I_syn_GABA_A': data['I_syn_GABA_A'],
                                        #'I_syn_GABA_B': data['I_syn_GABA_B'],
                                        # 'g_AMPA': data['g_AMPA'],
                                        # 'g_NMDA': data['g_NMDA'],
                                        # 'g_GABA_A': data['g_GABA_A'],
                                        # 'g_GABA_B': data['g_GABA_B']
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