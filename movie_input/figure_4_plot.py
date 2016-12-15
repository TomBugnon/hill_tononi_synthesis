#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Simulation of figure 4 of Hill-Tononi paper
# Author: Pablo Martinez Cañada (pablomc@ugr.es)
#
# During the first time interval, peak conductances take their initial values.
# In the second time interval, the transition from the waking mode to the sleep
# mode is simulated by changing conductances to half of their maximum values.
# Conductances are set to their maximum values during the sleep mode.

import nest
import nest.topology as tp
import numpy as np
import os
import matplotlib.pyplot as plt

import pylab
import pickle
import os.path
import scipy.io
from nest import raster_plot

import plotting
reload(plotting)

def simulation(Params):


    #! =================
    #! Import network
    #! =================

    # NEST Kernel and Network settings
    nest.ResetKernel()
    nest.ResetNetwork()
    nest.SetKernelStatus({"local_num_threads": Params['threads'],'resolution': Params['resolution']})
    nest.SetStatus([0],{'print_time': True})

    # import network description
    # import network
    # reload(network)
    # models, layers, conns  = network.get_Network(Params)

    import network_full_keiko
    reload(network_full_keiko)
    models, layers, conns  = network_full_keiko.get_Network(Params)

    # Create models
    for m in models:
            nest.CopyModel(m[0], m[1], m[2])

    # Create layers, store layer info in Python variable
    for l in layers:
            exec '%s = tp.CreateLayer(l[1])' % l[0]

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


    # initiate network activity
    nest.Simulate(500.0)
    #nest.Simulate(100.0)


    #! =================
    #! Recording devices
    #! =================

    nest.CopyModel('multimeter', 'RecordingNode',
            params = {'interval'   : Params['resolution'],
            'record_from': ['V_m',
                            'I_syn_AMPA',
                            'I_syn_NMDA',
                            'I_syn_GABA_A',
                            'I_syn_GABA_B',
                            'g_AMPA',
                            'g_NMDA',
                            'g_GABAA',
                            'g_GABAB',
                            'I_h',
                            'I_NaP',
                            'I_KNa'],
            'record_to'  : ['memory'],
            'withgid'    : True,
            'withtime'   : False})

    recorders = []

    nest.CopyModel('multimeter', 'RecordingNodeIntrinsic',
            params = {'interval'   : Params['resolution'],
            'record_from': ['V_m',
                            'I_syn_AMPA',
                            'I_syn_NMDA',
                            'I_syn_GABA_A',
                            'I_syn_GABA_B',
                            'g_AMPA',
                            'g_NMDA',
                            'g_GABAA',
                            'g_GABAB',
                            'I_h',
                            'I_NaP',
                            'I_KNa'],
            'record_to'  : ['memory'],
            'withgid'    : True,
            'withtime'   : True})

    recorders2 = []

    for population, model in [(Retina_layer, 'Retina'),
                        (Tp_layer  , 'Tp_exc'),
                        (Rp_layer  , 'Rp'),
                        (Vp_vertical, 'L23_exc'),
                        (Vp_vertical, 'L23_inh'),
                        (Vp_vertical, 'L4_exc'),
                        (Vp_vertical, 'L56_exc')]:


            rec = nest.Create('RecordingNode')
            recorders.append([rec,population,model])
            if (model=='Retina'):
                    nest.SetStatus(rec,{'record_from': ['rate']})
            tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
            nest.Connect(rec, tgts)


    for population, model in [(Vp_vertical, 'L23_exc')]:
            rec = nest.Create('RecordingNodeIntrinsic')
            recorders2.append([rec,population,model])
            tgts = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]
            nest.Connect(rec, tgts)


    #! =================
    #! Spike detector
    #! =================
    detectors = []
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
                              (Vp_horizontal, 'L56_inh')]:

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

    nest.SetStatus([0],{'print_time': True})
    t_sim = 0

    for t in Params['intervals']:
        if (t_sim==0):
            gKL = 1.0
            NaP_g_peak = 0.5
            h_g_peak = 1.0
            T_g_peak = 1.0
            KNa_g_peak = 0.5
            I_syn_AMPA_gain = 1.0  # keiko added this parameter in ht_neuron.cpp
            print(t_sim)
        '''
        # Nov23 comment-out
        if (t_sim==1):
            gKL = 1.0 + 0.8/2.0
            NaP_g_peak = 0.5 + 0.75/2.0
            h_g_peak = 1.0 + 1.0/2.0
            T_g_peak = 1.0 + 0.25/2.0
            KNa_g_peak = 0.5 + 0.75/2.0
            I_syn_AMPA_gain = 1.0 + 0.5/2.0 # keiko added this parameter in ht_neuron.cpp
            #I_syn_AMPA_gain = 1.0 + 0.25/2.0 # keiko added this parameter in ht_neuron.cpp
            print(t_sim)
        '''
        if (t_sim==1):
            gKL = 1.0 + 0.8*(1.0/4.0)
            NaP_g_peak = 0.5 + 0.75*(1.0/4.0)
            h_g_peak = 1.0 + 1.0*(1.0/4.0)
            T_g_peak = 1.0 + 0.25*(1.0/4.0)
            KNa_g_peak = 0.5 + 0.75*(1.0/4.0)
            I_syn_AMPA_gain = 1.0 + 0.25*(1.0/4.0) # keiko added this parameter in ht_neuron.cpp
            #I_syn_AMPA_gain = 1.0 + 0.25/2.0 # keiko added this parameter in ht_neuron.cpp
            print(t_sim)

        if (t_sim==2):
            gKL = 1.0 + 0.8*(2.0/4.0)
            NaP_g_peak = 0.5 + 0.75*(2.0/4.0)
            h_g_peak = 1.0 + 1.0*(2.0/4.0)
            T_g_peak = 1.0 + 0.25*(2.0/4.0)
            KNa_g_peak = 0.5 + 0.75*(2.0/4.0)
            I_syn_AMPA_gain = 1.0 + 0.25*(2.0/4.0) # keiko added this parameter in ht_neuron.cpp
            #I_syn_AMPA_gain = 1.0 + 0.25/2.0 # keiko added this parameter in ht_neuron.cpp
            print(t_sim)

        if (t_sim==3):
            gKL = 1.0 + 0.8*(3.0/4.0)
            NaP_g_peak = 0.5 + 0.75*(3.0/4.0)
            h_g_peak = 1.0 + 1.0*(3.0/4.0)
            T_g_peak = 1.0 + 0.25*(3.0/4.0)
            KNa_g_peak = 0.5 + 0.75*(3.0/4.0)
            I_syn_AMPA_gain = 1.0 + 0.25*(3.0/4.0) # keiko added this parameter in ht_neuron.cpp
            #I_syn_AMPA_gain = 1.0 + 0.25/2.0 # keiko added this parameter in ht_neuron.cpp
            print(t_sim)

        if (t_sim==4):
            gKL = 1.0 + 0.8
            NaP_g_peak = 0.5 + 0.75
            h_g_peak = 1.0 + 1.0
            T_g_peak = 1.0 + 0.25
            KNa_g_peak = 0.5 + 0.75
            I_syn_AMPA_gain = 1.0 + 0.25  # keiko added this parameter in ht_neuron.cpp
            #I_syn_AMPA_gain = 1.0 + 0.25  # keiko added this parameter in ht_neuron.cpp
            print(t_sim)

        t_sim+=1


        # change conductances in all populations
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
                                nest.SetStatus([cell], {'g_KL': gKL})

                                if (nest.GetStatus([cell],'NaP_g_peak')[0])>0.0 :
                                    nest.SetStatus([cell], {'NaP_g_peak':NaP_g_peak})
                                if (nest.GetStatus([cell],'h_g_peak')[0])>0.0:
                                    nest.SetStatus([cell], {'h_g_peak':h_g_peak})
                                if (nest.GetStatus([cell],'T_g_peak')[0])>0.0:
                                    nest.SetStatus([cell], {'T_g_peak':T_g_peak})
                                if (nest.GetStatus([cell],'KNa_g_peak')[0])>0.0 :
                                    nest.SetStatus([cell], {'KNa_g_peak':KNa_g_peak})
                                #if((nest.GetStatus([cell],'I_syn_AMPA_gain')[0])>0.0):
                                #   nest.SetStatus([cell], {'I_syn_AMPA_gain': I_syn_AMPA_gain}) # keiko

                        if l[0].__contains__('Vp') or l[0].__contains__('Vs'):
                            #print('%s %s' %(l[0], sim_model))
                            for cell in pop:
                                nest.SetStatus([cell], {'I_syn_AMPA_gain': I_syn_AMPA_gain})  # keiko

        # simulate interval
        nest.Simulate(t)



    #! ====================
    #! Plot Results
    #! ====================

    print "plotting..."

    rows = 11
    cols = 2

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4)

    # Plot A: membrane potential rasters

    recorded_models = [(Retina_layer,'Retina'),
                        (Vp_vertical,'L23_exc'),
                        (Vp_vertical,'L4_exc'),
                        (Vp_vertical,'L56_exc'),
                        (Rp_layer,'Rp'),
                        (Tp_layer,'Tp_exc')]

    plotting.potential_raster(fig,recorders,recorded_models,100,3*Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,0)

    # Plot B: individual intracellular traces

    recorded_models =[(Vp_vertical,'L23_exc'),
                (Vp_vertical,'L23_inh'),
                (Rp_layer,'Rp'),
                (Tp_layer,'Tp_exc')]

    # original
    #plotting.intracellular_potentials(fig,recorders,recorded_models,100,rows,cols,6)
    # keiko
    total_time = 0.0
    for t in Params['intervals']:
        total_time += t
    plotting.intracellular_potentials(fig,recorders,recorded_models,100,rows,cols,6,total_time)

    # Plot C: topographical activity of the up- and down-states

    recorded_models = [(Vp_vertical,'L23_exc')]

    labels = ["Upstate"]
    start = 2500.0
    stop = 2510.0
    #start = 900.0
    #stop = 910.0

    plotting.topographic_representation(fig,recorders,recorded_models,labels,Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,start,stop,10,0)

    labels = ["Downstate"]
    start = 2600.0
    stop = 2610.0
    #start = 1550.0
    #stop = 1560.0

    plotting.topographic_representation(fig,recorders,recorded_models,labels,Params['Np'],np.sum(Params['intervals']),Params['resolution'],rows,cols,start,stop,10,1)

    # Plot D: Intrinsic currents of a selected cell

    recorded_models =[(Vp_vertical,'L23_exc')]
    plotting.intrinsic_currents(recorders2, recorded_models, 100)
    plotting.synaptic_currents(recorders2, recorded_models, 100)

    plt.show()

    # Plot E: movie

    #labels = ["Osc_Vp_L23_Vertical"]
    #recorded_models = [(Vp_vertical,'L23_exc')]
    #plotting.makeMovie(fig,recorders,recorded_models,labels,Params['Np'],np.sum(Params['intervals']),Params['resolution'])

    #! ====================
    #! Save Results
    #! ====================

    print('save recorders data')
    # Set folder
    rootdir = '/home/kfujii2/newNEST2/ht_model_pablo_based/data/'
    expdir = ''
    # expdir = 'random_full/'
    # expdir = 'structured_full/'
    data_folder = rootdir + expdir
    if not os.path.isdir(data_folder):
        os.makedirs(data_folder)

    # To save spike data, set pairs of population id and its name
    population_name = [ {'population': Retina_layer, 'name': 'Retina'},
                        {'population': Vp_vertical, 'name': 'Vp_v'},
                        {'population': Vp_horizontal, 'name': 'Vp_h'},
                        {'population': Rp_layer, 'name': 'Rp'},
                        {'population': Tp_layer, 'name': 'Tp'} ]

    # vertical
    population = population_name[1]['population']
    name = population_name[1]['name']
    model_list = ['L23_exc', 'L56_exc']
    for model in model_list:
        l = [nd for nd in np.arange(0,len(recorders)) if (recorders[nd][1] == population and recorders[nd][2] == model)][0]
        data = nest.GetStatus(recorders[l][0],keys='events')[0]
        scipy.io.savemat(data_folder + '/recorder_' + name + '_' + model + '.mat',
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

        raster = raster_plot.from_device(rec, hist=True)
        pylab.title( p_name + '_' + model )
        f = raster[0].figure
        f.set_size_inches(15, 9)
        f.savefig(data_folder + 'spikes_' + p_name + '_' + model + '.png', dpi=100)
        plt.close()
        '''
        # Set filename and save spike data
        filename = data_folder + 'spike_' + p_name + '_' + model + '.pickle'
        pickle.dump(spikes, open(filename, 'w'))

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
