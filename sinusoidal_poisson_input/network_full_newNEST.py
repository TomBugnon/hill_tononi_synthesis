#!/usr/bin/env python
# -*- coding: utf-8 -*-
import nest
import numpy as np

maskFactor = 1.0
#weight_gain = 1.0
weight_gain = 2.0/3.0

# Update dictionaries
def updateDicts(dict1, dict2):
    assert(isinstance(dict1, dict))
    assert(isinstance(dict2, dict))

    tmp = dict1.copy()
    tmp.update(dict2)
    return tmp

# Return lists of layers, models and connections
def get_Network(params):
    models = get_Models(params)
    layers = get_Layers(params)
    conns = get_Connections(params)
    return models,layers,conns

# keiko keep
'''
"AMPA_g_peak": 0.075,
"NMDA_g_peak": 0.01,
"GABA_A_g_peak": 0.15,
"GABA_B_g_peak": 0.01
'''
# original
'''
"AMPA_g_peak": 0.1,
"NMDA_g_peak": 0.075,
"GABA_A_g_peak": 0.33,
"GABA_B_g_peak": 0.0132
'''

# Configuration of the ht neuron model for each layer
def get_Models(params):

    nrnmod  = 'ht_neuron'

    #test Tom ### CRITICAL
    if params['synapse_old']:
        nest.SetDefaults('ht_synapse', {'tau_P' : 50., 'delta_P': 0.2})  #new defaults are tau_P = 500. and delta_P = 0.125
    if params['NMDA_old']:
        nest.SetDefaults('ht_neuron', params = { 'instant_unblock_NMDA' : True , 'S_act_NMDA' : 1./2.5, 'V_act_NMDA' : -58.}) #new defaults are instant_unblock_NMDA = False, S_act_NMDA = 0.081 and V_act_NMDA = -25.57


    cortical_excitatory = {
        "theta_eq": -51.0,
        "tau_theta": 2.0,
#        "spike_duration": 2.0,
        "tau_spike": 1.75,
        "tau_m": 16.0,
        "g_peak_NaP": 0.5,
        "g_peak_h": 0.0,
        "g_peak_T": 0.0,
        "g_peak_KNa": 0.5,
        "g_KL": 1.0,
        "E_rev_NaP": 55.0,  # nov23
        #
        "g_peak_AMPA": 0.1,
        "g_peak_NMDA": 0.075 * 0.2,
        #"NMDA_g_peak": 0.075, #original
        "g_peak_GABA_A": 0.33,
        "g_peak_GABA_B": 0.0132
        #"AMPA_g_peak": 0.075,
        #"NMDA_g_peak": 0.01,
        #"GABA_A_g_peak": 0.15,
        #"GABA_B_g_peak": 0.01
    }

    cortical_inhibitory = {
        "theta_eq": -53.0,
        "tau_theta": 1.0,
#        "spike_duration": 1.0,
        "tau_spike": 0.5,
        "tau_m": 8.0,
        "g_peak_NaP": 0.5,
        "g_peak_h": 0.0,
        "g_peak_T": 0.0,
        "g_peak_KNa": 0.5,
        "g_KL": 1.0,
        "E_rev_NaP": 55.0,  # nov23
        #
        "g_peak_AMPA": 0.1,
        "g_peak_NMDA": 0.075 * 0.2,
        #"NMDA_g_peak": 0.075, #original
        "g_peak_GABA_A": 0.33,
        "g_peak_GABA_B": 0.0132
        #"AMPA_g_peak": 0.075,
        #"NMDA_g_peak": 0.01,
        #"GABA_A_g_peak": 0.15,
        #"GABA_B_g_peak": 0.01
    }

    thalamic = {
        "theta_eq": -53.0,
        "tau_theta": 0.75,
#        "spike_duration": 1.0,
        "tau_spike": 0.75,
        "tau_m": 8.0,
        "E_rev_GABA_A": -70.0,
        "g_peak_NaP": 0.5,
        "g_peak_h": 1.0,
        "g_peak_T": 1.0,
        "g_peak_KNa": 0.0,
        "g_KL": 1.0,
        "E_rev_NaP": 55.0,  # nov23
        #
        "g_peak_AMPA": 0.1,
        "g_peak_NMDA": 0.075 * 0.2,
        #"NMDA_g_peak": 0.075, #original
        "g_peak_GABA_A": 0.33,
        "g_peak_GABA_B": 0.0132
        #"AMPA_g_peak": 0.075,
        #"NMDA_g_peak": 0.01,
        #"GABA_A_g_peak": 0.15,
        #"GABA_B_g_peak": 0.01
    }


    reticular = {
        "theta_eq": -53.0,
        "tau_theta": 0.75,
#        "spike_duration": 1.0,
        "tau_spike": 0.75,
        "tau_m": 8.0,
        "E_rev_GABA_A": -80.0,
        "g_peak_NaP": 0.5,
        "g_peak_h": 1.0,
        "g_peak_T": 1.0,
        "g_peak_KNa": 0.0,
        "g_KL": 1.0,
        "E_rev_NaP": 55.0,  # nov23
        #
        "g_peak_AMPA": 0.1,
        "g_peak_NMDA": 0.075 * 0.2,
        #"NMDA_g_peak": 0.075, #original
        "g_peak_GABA_A": 0.33,
        "g_peak_GABA_B": 0.0132
        #"AMPA_g_peak": 0.075,
        #"NMDA_g_peak": 0.01,
        #"GABA_A_g_peak": 0.15,
        #"GABA_B_g_peak": 0.01
    }

    models = [
        (nrnmod, 'Tp_exc', thalamic),
        (nrnmod, 'Tp_inh', thalamic),
        (nrnmod, 'Rp', reticular), #original
        #(nrnmod, 'Rp', thalamic),
        ('sinusoidal_poisson_generator', 'Retina', {'individual_spike_trains': True,
                                                    'amplitude': params['ret_amplitude'],
                                                    'rate': params['ret_rate'],
                                                    'frequency': params['temporal_frequency'],
                                                    'phase': 0.0})
    ]

    models += [(nrnmod, l+'_exc', cortical_excitatory) for l in ('L23','L4','L56')]
    models += [(nrnmod, l+'_inh' , cortical_inhibitory) for l in ('L23','L4','L56')]

    return models



# Layer properties
def get_Layers(params):

    if params.has_key('edge_wrap') and params['edge_wrap']:
        edge_wrap = params['edge_wrap']
    else:
        # edge_wrap = True #original
        edge_wrap = False #keiko


    # Primary pathway
    layerPropsP = {
    'rows'     : params['Np'],
    'columns'  : params['Np'],
    'extent'   : [params['visSize'], params['visSize']],
    'edge_wrap': edge_wrap
    }

    # Secondary pathway
    layerPropsS = {
    'rows'     : params['Ns'],
    'columns'  : params['Ns'],
    'extent'   : [params['visSize'], params['visSize']],
    'edge_wrap': edge_wrap
    }

    # Create layer dictionaries
    layers = [
    ('Retina_layer', updateDicts(layerPropsP, {'elements': 'Retina'})),
    ('Tp_layer'    , updateDicts(layerPropsP, {'elements': ['Tp_exc', 'Tp_inh']})),
    ('Rp_layer'    , updateDicts(layerPropsP, {'elements': 'Rp'})),
    ('Vp_horizontal', updateDicts(layerPropsP, {'elements': ['L23_exc', 2, 'L23_inh', 1,
                                                             'L4_exc' , 2, 'L4_inh' , 1,
                                                             'L56_exc', 2, 'L56_inh', 1]})),
    ('Vp_vertical', updateDicts(layerPropsP, {'elements': ['L23_exc', 2, 'L23_inh', 1,
                                                           'L4_exc' , 2, 'L4_inh' , 1,
                                                           'L56_exc', 2, 'L56_inh', 1]})),
    ('Ts_layer'    , updateDicts(layerPropsS, {'elements': ['Tp_exc', 'Tp_inh']})),
    ('Rs_layer'    , updateDicts(layerPropsS, {'elements': 'Rp'})),
    ('Vs_horizontal'  , updateDicts(layerPropsS, {'elements': ['L23_exc', 2, 'L23_inh', 1,
                                                               'L4_exc' , 2, 'L4_inh' , 1,
                                                               'L56_exc', 2, 'L56_inh', 1]})),
    ('Vs_cross'  , updateDicts(layerPropsS, {'elements': ['L23_exc', 2, 'L23_inh', 1,
                                                          'L4_exc' , 2, 'L4_inh' , 1,
                                                          'L56_exc', 2, 'L56_inh', 1]})),
    ('Vs_vertical'  , updateDicts(layerPropsS, {'elements': ['L23_exc', 2, 'L23_inh', 1,
                                                             'L4_exc' , 2, 'L4_inh' , 1,
                                                             'L56_exc', 2, 'L56_inh', 1]})) ]

    return layers


def scramble_connections(dict, params):

    v_size = params['visSize']
    if params['scrambled']:
        dict.update({"mask": {"rectangular": {"lower_left" : [-v_size/2., -v_size/2.],"upper_right": [-v_size/2., -v_size/2.]}}})
        dict.update({"mask": {"rectangular": {"lower_left" : [-v_size/2., -v_size/2.],"upper_right": [-v_size/2., -v_size/2.]}}})


# Connections between layers
def get_Connections(params):

    # Scaling parameters from grid elements to visual angle
    dpcP = maskFactor * params['visSize'] / (params['Np'] - 1)
    dpcS = params['visSize'] / (params['Ns'] - 1)

    # Mapping of receptor names to receptor indices from the ht_neuron
    receptors = nest.GetDefaults('ht_neuron')['receptor_types']

    # Synapses
    nest.CopyModel('ht_synapse', 'AMPA_syn', {'receptor_type': receptors['AMPA']})
    nest.CopyModel('ht_synapse', 'NMDA_syn', {'receptor_type': receptors['NMDA']})
    nest.CopyModel('ht_synapse', 'GABA_A_syn', {'receptor_type': receptors['GABA_A']})
    nest.CopyModel('ht_synapse', 'GABA_B_syn', {'receptor_type': receptors['GABA_B']})


    # build complete list of connections
    allconns = []

    # --------------------------------------------------------------------#
    # ---------- PRIMARY VISUAL AREA ------------------------------------ #
    # --------------------------------------------------------------------#

    # Cortico-cortical, same orientation
    ccConnections = []
    # Cortico-cortical, cross-orientation
    ccxConnections = []

    # Horizontal intralaminar connections
    Vp_horizontal_intralaminar_base = {
        "connection_type": "divergent",
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 12.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.05, "sigma": 7.5 * dpcP}},
        #"weights": 1.0, # original
        "weights": 1.0 * weight_gain,
        #"weights": 0.75 * weight_gain, # keiko
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    # Vp_horizontal_intralaminar_base = scramble_connections(Vp_horizontal_intralaminar_base,
    #                                                        params)
    for conn in [#{"sources": {"model": "L23_exc"}, "targets": {"model": "L23_exc"}, "synapse_model": "NMDA_syn"}, #original
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L23_exc"}, "synapse_model": "NMDA_syn", "weights": 1.0*weight_gain}, # keiko
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L23_exc"}},
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L23_inh"}},
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L4_exc"},
                  "mask"   : {"circular": {"radius": 7.0 * dpcP}}},
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L4_inh"},
                  "mask"   : {"circular": {"radius": 7.0 * dpcP}}},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L56_exc"}},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L56_inh"}}]:
        ndict = Vp_horizontal_intralaminar_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    # Vertical interlaminar connections
    Vp_vertical_interlaminar_base = {
        "connection_type": "divergent",
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 2.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 1.0, "sigma": 7.5 * dpcP}},
        #"weights": 2.0, # original
        "weights": 2.0 * weight_gain,
        #"weights": 1.5 * weight_gain, # keiko
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    for conn in [#{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "synapse_model": "NMDA_syn"},  #original
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "synapse_model": "NMDA_syn", "weights": 1.0*weight_gain}, # keiko
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "weights": 1.0},  # original
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_inh"}, "weights": 1.0},  # original
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "weights": 1.0 * weight_gain},  # keiko
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L56_inh"}, "weights": 1.0 * weight_gain},  # keiko
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L23_exc"}},
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L23_inh" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L23_exc"}},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L23_inh" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L4_exc" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L4_inh"  }}]:
        ndict = Vp_vertical_interlaminar_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    # Intracortical inhibitory connections
    Vp_intracortical_inhibitory_base = {
        "connection_type": "divergent",
        "synapse_model": "GABA_A_syn",
        "mask": {"circular": {"radius": 7.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.25, "sigma": 7.5 * dpcP}},
        # "weights": 1.0, # original
        "weights": 1.0 * weight_gain, #
        #"weights": 2.0 * weight_gain, # keiko
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    for conn in [{"sources": {"model": "L23_inh"},
                  "targets": {"model": "L23_exc"},
                  "synapse_model": "GABA_B_syn",
                  #"mask": {"circular": {"radius": 1.0 * dpcP}},  # original
                  "mask": {"circular": {"radius": 2.0 * dpcP}},  # keiko
                  "kernel": 0.3},
                 {"sources": {"model": "L23_inh" },
                  "targets": {"model": "L4_exc" },
                  "synapse_model": "GABA_B_syn",
                  #"mask": {"circular": {"radius": 1.0 * dpcP}}, # original
                  "mask": {"circular": {"radius": 2.0 * dpcP}},  # keiko
                  "kernel": 0.3},
                 {"sources": {"model": "L23_inh"},
                  "targets": {"model": "L56_exc" },
                  "synapse_model": "GABA_B_syn",
                  #"mask": {"circular": {"radius": 1.0 * dpcP}}, #original
                  "mask": {"circular": {"radius": 2.0 * dpcP}},  # keiko
                  "kernel": 0.3}]:
        ndict = Vp_intracortical_inhibitory_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    '''
    # nov23 comment-out
    for conn in [{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_exc"}},
                 #{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}},  #original
                 {"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}, "weights": 0.5 * weight_gain},  # keiko
                 {"sources": {"model": "L4_inh" }, "targets": {"model": "L4_exc" }},
                 #{"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}},  #original
                 {"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}, "weights": 0.5 * weight_gain },  # keiko
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_exc"}},
                 #{"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}}]:  #original
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh" }, "weights": 0.5 * weight_gain}]: # keiko
        ndict = Vp_intracortical_inhibitory_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)
        ccxConnections.append(ndict)
    '''

    # for same orientation, inh->inh = 0.5
    for conn in [{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_exc"}},
                 #{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}},  #original
                 {"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}, "weights": 0.5 * weight_gain},  # keiko
                 {"sources": {"model": "L4_inh" }, "targets": {"model": "L4_exc" }},
                 #{"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}},  #original
                 {"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}, "weights": 0.5 * weight_gain },  # keiko
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_exc"}},
                 #{"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}}]:  #original
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh" }, "weights": 0.5 * weight_gain}]: # keiko
        ndict = Vp_intracortical_inhibitory_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    # for different orientation, inh->inh = 0.75
    for conn in [{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_exc"}},
                 # {"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}},  #original
                 {"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}, "weights": 0.75 * weight_gain},
                 # keiko
                 {"sources": {"model": "L4_inh"}, "targets": {"model": "L4_exc"}},
                 # {"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}},  #original
                 {"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}, "weights": 0.75 * weight_gain},
                 # keiko
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_exc"}},
                 # {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}}]:  #original
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}, "weights": 0.75 * weight_gain}]:  # keiko
        ndict = Vp_intracortical_inhibitory_base.copy()
        ndict.update(conn)
        ccxConnections.append(ndict)

    # Leonardo:
    # For connections within populations with same orientation selectivity, connect twice so that it is
    # easier to scramble network and remove orientation selectivity if scramble = True.
    # Notice: probabilities have been divided by 2, so the final result should be the same
    # as before when scrambled == False, as long allow_multapses == False.

    #! Cortico-cortical, same orientation
    [allconns.append(['Vp_horizontal','Vp_horizontal',c]) for c in ccConnections]
    [allconns.append(['Vp_vertical','Vp_vertical',c]) for c in ccConnections]

    #! Cortico-cortical, cross-orientation
    [allconns.append(['Vp_horizontal','Vp_vertical',c]) for c in ccxConnections]
    [allconns.append(['Vp_vertical','Vp_horizontal',c]) for c in ccxConnections]


    # --------------------------------------------------------------------#
    # ---------- SECONDARY VISUAL AREA ---------------------------------- #
    # --------------------------------------------------------------------#

    # Cortico-cortical, same orientation
    ccConnections = []
    # Cortico-cortical, cross-orientation
    ccxConnections = []

    # Horizontal intralaminar connections
    Vs_horizontal_intralaminar_base = {
        "connection_type":"divergent",
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 12.0 * dpcS}},
        "kernel": {"gaussian": {"p_center": 0.05, "sigma": 7.5 * dpcP}},
        #"weights": 1.0, # original
        "weights": 1.0 * weight_gain,
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    for conn in [#{"sources": {"model": "L23_exc"}, "targets": {"model": "L23_exc"}, "synapse_model": "NMDA_syn"}, #original
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L23_exc"}, "synapse_model": "NMDA_syn", "weights": 1.0*weight_gain}, #keiko
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L23_exc"}},
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L23_inh" }},
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L4_exc" },
                  "mask"   : {"circular": {"radius": 7.0 * dpcS}}},
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L4_inh"  },
                  "mask"   : {"circular": {"radius": 7.0 * dpcS}}},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L56_exc" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L56_inh"  }}]:
        ndict = Vs_horizontal_intralaminar_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    # Vertical interlaminar connections
    Vs_vertical_interlaminar_base = {
        "connection_type":"divergent",
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 2.0 * dpcS}},
        "kernel": {"gaussian": {"p_center": 1.0, "sigma": 7.5 * dpcS}},
        #"weights": 2.0, # original
        #"weights": 2.0 * weight_gain,
        "weights": 1.0 * weight_gain, # keiko
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    for conn in [#{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "synapse_model": "NMDA_syn", "weights": 1.0}, #original
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "weights": 1.0}, #original
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_inh" }, "weights": 1.0}, #original
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "synapse_model": "NMDA_syn", "weights": 1.0 * weight_gain},
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "synapse_model": "NMDA_syn", "weights": 1.0*weight_gain},
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "weights": 1.0*weight_gain}, #original
                 {"sources": {"model": "L23_exc"}, "targets": {"model": "L56_inh"}, "weights": 1.0*weight_gain}, #original
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_exc"}, "weights": 0.75*weight_gain}, #keiko
                 #{"sources": {"model": "L23_exc"}, "targets": {"model": "L56_inh"}, "weights": 0.75*weight_gain}, #keiko
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L23_exc"}},
                 {"sources": {"model": "L4_exc" }, "targets": {"model": "L23_inh" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L23_exc"}},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L23_inh" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L4_exc" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "L4_inh"  }}]:
        ndict = Vs_vertical_interlaminar_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    # Intracortical inhibitory connections
    Vs_intracortical_inhibitory_base = {
        "connection_type":"divergent",
        "synapse_model": "GABA_A_syn",
        "mask": {"circular": {"radius": 7.0 * dpcS}},
        "kernel": {"gaussian": {"p_center": 0.25, "sigma": 7.5 * dpcS}},
        #"weights": 1.0,  # original
        "weights": 1.0 * weight_gain,
        #"weights": 1.5 * weight_gain, # keiko
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    for conn in [{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_exc"},
                  #"synapse_model": "GABA_B_syn", "mask": {"circular": {"radius": 1.0 * dpcS}}, "kernel": 0.3}, #original
                  "synapse_model": "GABA_B_syn", "mask": {"circular": {"radius": 2.0 * dpcS}}, "kernel": 0.3},  #keiko
                 {"sources": {"model": "L23_inh"}, "targets": {"model": "L4_exc"},
                  #"synapse_model": "GABA_B_syn", "mask": {"circular": {"radius": 1.0 * dpcS}}, "kernel": 0.3},  #original
                  "synapse_model": "GABA_B_syn", "mask": {"circular": {"radius": 2.0 * dpcS}}, "kernel": 0.3},  # keiko
                 {"sources": {"model": "L23_inh"}, "targets": {"model": "L56_exc"},
                  #"synapse_model": "GABA_B_syn", "mask": {"circular": {"radius": 1.0 * dpcS}}, "kernel": 0.3}]:
                  "synapse_model": "GABA_B_syn", "mask": {"circular": {"radius": 2.0 * dpcS}}, "kernel": 0.3}]:  # keiko

        ndict = Vs_intracortical_inhibitory_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)

    for conn in [{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_exc"}},
                 #{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}},  #original
                 #{"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}, "weights": 1.0 * weight_gain}, # keiko
                 {"sources": {"model": "L23_inh"}, "targets": {"model": "L23_inh"}, "weights": 0.5 * weight_gain}, #keiko nov22
                 {"sources": {"model": "L4_inh" }, "targets": {"model": "L4_exc" }},
                 #{"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}},  #original
                 #{"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}, "weights": 1.0 * weight_gain}, # keiko
                 {"sources": {"model": "L4_inh"}, "targets": {"model": "L4_inh"}, "weights": 0.5 * weight_gain}, #keiko nov22
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_exc"}},
                 #{"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}} ]: #original
                 #{"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}, "weights": 1.0 * weight_gain}]:  # keiko
                 {"sources": {"model": "L56_inh"}, "targets": {"model": "L56_inh"}, "weights": 0.5 * weight_gain}]:  # keiko nov22
        ndict = Vs_intracortical_inhibitory_base.copy()
        ndict.update(conn)
        ccConnections.append(ndict)
        ccxConnections.append(ndict)

    # Leonardo:
    # For connections within populations with same orientation selectivity, connect twice so that it is
    # easier to scramble network and remove orientation selectivity if scramble = True.
    # Notice: probabilities have been divided by 2, so the final result should be the same
    # as before when scrambled == False, as long allow_multapses == False.

    #! Cortico-cortical, same orientation
    [allconns.append(['Vs_horizontal','Vs_horizontal',c]) for c in ccConnections]
    [allconns.append(['Vs_vertical','Vs_vertical',c]) for c in ccConnections]
    [allconns.append(['Vs_cross','Vs_cross',c]) for c in ccConnections]

    #! Cortico-cortical, cross-orientation

    #--- original
    #[allconns.append(['Vs_horizontal','Vs_vertical',c]) for c in ccxConnections]
    #[allconns.append(['Vs_horizontal','Vs_cross',c]) for c in ccxConnections]
    #[allconns.append(['Vs_vertical','Vs_horizontal',c]) for c in ccxConnections]
    #[allconns.append(['Vs_vertical','Vs_cross',c]) for c in ccxConnections]
    #[allconns.append(['Vs_cross','Vs_vertical',c]) for c in ccxConnections]
    #[allconns.append(['Vs_cross','Vs_horizontal',c]) for c in ccxConnections]

    #--- keiko 2017/1/9
    [allconns.append(['Vs_horizontal', 'Vs_vertical', c]) for c in ccxConnections] # 1/2 H->V
    [allconns.append(['Vs_vertical', 'Vs_horizontal', c]) for c in ccxConnections] # 1/2 V->H

    [allconns.append(['Vs_horizontal', 'Vs_cross', c]) for c in ccxConnections]
    [allconns.append(['Vs_vertical', 'Vs_cross', c]) for c in ccxConnections]
    [allconns.append(['Vs_cross', 'Vs_vertical', c]) for c in ccxConnections]
    [allconns.append(['Vs_cross', 'Vs_horizontal', c]) for c in ccxConnections]


    # --------------------------------------------------------------------#
    # ---------- INTERAREAL CONNECTIONS --------------------------------- #
    # --------------------------------------------------------------------#

    # Forward interareal connections (mask is in Vp, thus scaling with dpcP)
    Vp_Vs_forward_interareal_base = {
        "connection_type":"convergent",
        "sources": {"model": "L23_exc"},
        "kernel": 0.8,
        #"weights": 1.0, # original
        "weights": 1.0 * weight_gain,
        "delays": {"uniform": {"min": 2.75, "max": 3.25}},
        "synapse_model": "AMPA_syn",
        "allow_autapses": False,
        "allow_multapses": False
    }

    fwdInterConns = [('Vp_horizontal', 'Vs_horizontal', {"targets": {"model": "L4_exc"},
                                                         "mask": {"rectangular": {"lower_left" : [-9.5*dpcP, -1.5*dpcP],"upper_right": [ 9.5*dpcP,  1.5*dpcP]}}}),
                     ('Vp_horizontal', 'Vs_horizontal', {"targets": {"model": "L4_inh"},
                                                         "mask": {"rectangular": {"lower_left" : [-9.5*dpcP, -1.5*dpcP],"upper_right": [ 9.5*dpcP,  1.5*dpcP]}}}),
                     ('Vp_vertical', 'Vs_vertical', {"targets": {"model": "L4_exc"},
                                                     "mask": {"rectangular": {"lower_left" : [-1.5*dpcP, -9.5*dpcP],"upper_right": [ 1.5*dpcP,  9.5*dpcP]}}}),
                     ('Vp_vertical', 'Vs_vertical', {"targets": {"model": "L4_inh"},
                                                     "mask": {"rectangular": {"lower_left" : [-1.5*dpcP, -9.5*dpcP],"upper_right": [ 1.5*dpcP,  9.5*dpcP]}}}),
                     ('Vp_horizontal', 'Vs_cross', {"targets": {"model": "L4_exc"},
                                                    "mask": {"rectangular": {"lower_left" : [-4.5*dpcP, -1.5*dpcP], "upper_right": [ 4.5*dpcP,  1.5*dpcP]}}}),
                     ('Vp_horizontal', 'Vs_cross', {"targets": {"model": "L4_inh"},
                                                    "mask": {"rectangular": {"lower_left" : [-4.5*dpcP, -1.5*dpcP], "upper_right": [ 4.5*dpcP,  1.5*dpcP]}}}),
                     ('Vp_vertical', 'Vs_cross', {"targets": {"model": "L4_exc"},
                                                  "mask": {"rectangular": {"lower_left" : [-1.5*dpcP, -4.5*dpcP], "upper_right": [ 1.5*dpcP,  4.5*dpcP]}}}),
                     ('Vp_vertical', 'Vs_cross', {"targets": {"model": "L4_inh"},
                                                  "mask": {"rectangular": {"lower_left" : [-1.5*dpcP, -4.5*dpcP], "upper_right": [ 1.5*dpcP,  4.5*dpcP]}}}) ]

    allconns += [ (c[0], c[1], updateDicts(Vp_Vs_forward_interareal_base, c[2])) for c in fwdInterConns ]

    # Backward interareal connections (mask in Vp, thus scaling with dpcP)
    Vs_Vp_backward_interareal_base = {
        "connection_type":"divergent",
        "sources": {"model": "L56_exc"},
        "mask": {"circular": {"radius": 12.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.1, "sigma": 7.5 * dpcP}},
        #"weights": 1.0, # original
        "weights": 1.0 * weight_gain,
        "delays": {"uniform": {"min": 5.5, "max": 6.5}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    bckInterConns = [#('Vs_vertical', 'Vp_vertical', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn'}), #original
                     ('Vs_vertical', 'Vp_vertical', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), #keiko
                     #('Vs_vertical', 'Vp_vertical', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), #keiko
                     ('Vs_vertical', 'Vp_vertical', {'targets': {'model': 'L23_inh'}, 'synapse_model': 'AMPA_syn'}),
                     #('Vs_horizontal', 'Vp_horizontal', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn'}), #original
                     ('Vs_horizontal', 'Vp_horizontal', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), #keiko
                     #('Vs_horizontal', 'Vp_horizontal', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), #keiko
                     ('Vs_horizontal', 'Vp_horizontal', {'targets': {'model': 'L23_inh'}, 'synapse_model': 'AMPA_syn'}),
                     #('Vs_cross', 'Vp_vertical', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn'}), #original
                     ('Vs_cross', 'Vp_vertical', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), # keiko
                     #('Vs_cross', 'Vp_vertical', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), # keiko
                     ('Vs_cross', 'Vp_vertical', {'targets': {'model': 'L23_inh'}, 'synapse_model': 'AMPA_syn'}),
                     #('Vs_cross', 'Vp_horizontal', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn'}), # original
                     ('Vs_cross', 'Vp_horizontal', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), # keiko
                     #('Vs_cross', 'Vp_horizontal', {'targets': {'model': 'L23_exc'}, 'synapse_model': 'NMDA_syn', "weights": 1.0 * weight_gain}), # keiko
                     ('Vs_cross', 'Vp_horizontal', {'targets': {'model': 'L23_inh'}, 'synapse_model': 'AMPA_syn'})]

    allconns += [ (c[0], c[1], updateDicts(Vs_Vp_backward_interareal_base, c[2])) for c in bckInterConns ]

    # --------------------------------------------------------------------#
    # ---------- THALAMUS ----------------------------------------------- #
    # --------------------------------------------------------------------#

    # Thalamic connections
    thal_base = {
        "connection_type": "divergent",
        "delays": {"uniform": {"min": 1.75, "max": 2.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    for src, tgt, conn in [('Tp_layer', 'Rp_layer', {"sources": {"model": "Tp_exc"},
                                                     "targets": {"model": "Rp"},
                                                     "synapse_model": "AMPA_syn",
                                                     "mask": {"circular": {"radius": 2.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 1.0, "sigma": 7.5 * dpcP}},
                                                     #"weights": 2.0}),  # original
                                                     "weights": 2.0*weight_gain}),
                                                     #"weights": 1.0}),  # keiko
                                                     #"weights": 2.5*weight_gain}),  # keiko
                           ('Tp_layer', 'Tp_layer', {"sources": {"model": "Tp_inh"},
                                                     "targets": {"model": "Tp_exc"},
                                                     "synapse_model": "GABA_A_syn",
                                                     #"weights": 1.0, # original
                                                     "weights": 1.0*weight_gain,
                                                     "mask": {"circular": {"radius": 2.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 0.25, "sigma": 7.5 * dpcP}}}),
                           ('Tp_layer', 'Tp_layer', {"sources": {"model": "Tp_inh"},
                                                     "targets": {"model": "Tp_inh"},
                                                     "synapse_model": "GABA_A_syn",
                                                     #"weights": 1.0, #original
                                                     "weights": 1.0 *weight_gain,
                                                     "mask": {"circular": {"radius": 2.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 0.25, "sigma": 7.5 * dpcP}}}),
                           ('Rp_layer', 'Tp_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_exc"},
                                                     "synapse_model": "GABA_A_syn",
                                                     #"weights": 1.0,  # original
                                                     "weights": 1.0*weight_gain,
                                                     #"weights": 1.5 * weight_gain, #keiko
                                                     "mask": {"circular": {"radius": 12.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 0.15, "sigma": 7.5 * dpcP}}}),
                           ('Rp_layer', 'Tp_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_exc"},
                                                     "synapse_model": "GABA_B_syn",
                                                     #"weights": 1.0, #original
                                                     "weights": 1.0*weight_gain,
                                                     "mask": {"circular": {"radius": 12.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 0.05, "sigma": 7.5 * dpcP}}}),  # GABA_B
                           ('Rp_layer', 'Tp_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_inh"},
                                                     "synapse_model": "GABA_A_syn",
                                                     #"weights": 1.0, #original
                                                     "weights": 1.0 * weight_gain,  # keiko
                                                     #"weights": 1.5 * weight_gain,  # keiko
                                                     "mask": {"circular": {"radius": 12.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 0.15, "sigma": 7.5 * dpcP}}}),
                           ('Rp_layer', 'Rp_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Rp"},
                                                     "synapse_model": "GABA_B_syn",
                                                     #"weights": 1.0, # original
                                                     "weights": 1.0*weight_gain,
                                                     #"weights": 0.75, # keiko
                                                     #"weights": 0.5*weight_gain,  # keiko
                                                     "mask": {"circular": {"radius": 12.0 * dpcP}},
                                                     "kernel": {"gaussian": {"p_center": 0.5, "sigma": 7.5 * dpcP}}}) ]:
        thal = thal_base.copy()
        thal.update(conn)
        allconns.append([src, tgt, thal])
        # NOTE: Pablo's script has the follwoing connectin (Rp -> Tp_inh, GABA_B ), but this network doesn't.
        #('Rp_layer', 'Tp_layer', {"sources": {"model": "Rp"},"targets": {"model": "Tp_inh"},
        #"mask": {"rectangular": {"lower_left" : [-12.0* dpcP, -12.0* dpcP],"upper_right": [ 12.0* dpcP,  12.0* dpcP]}},
        #"synapse_model": "GABA_B_syn",
        #"weights": 1.0,
        #"kernel":  0.05}),

    for src, tgt, conn in [('Ts_layer', 'Rs_layer', {"sources": {"model": "Tp_exc"},
                                                     "targets": {"model": "Rp"},
                                                     "synapse_model": "AMPA_syn",
                                                     "mask": {"circular": {"radius": 2.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 1.0, "sigma": 7.5 * dpcS}},
                                                     #"weights": 2.0}), #original
                                                     "weights": 2.0*weight_gain}),
                           ('Ts_layer', 'Ts_layer', {"sources": {"model": "Tp_inh"},
                                                     "targets": {"model": "Tp_exc"},
                                                     "synapse_model": "GABA_A_syn",
                                                     "mask": {"circular": {"radius": 2.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.25, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0}), #origial
                                                     "weights": 1.0*weight_gain}),
                           ('Ts_layer', 'Ts_layer', {"sources": {"model": "Tp_inh"},
                                                     "targets": {"model": "Tp_inh"},
                                                     "synapse_model": "GABA_A_syn",
                                                     "mask": {"circular": {"radius": 2.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.25, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0}), # original
                                                     "weights": 1.0*weight_gain}),
                           ('Rs_layer', 'Ts_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_exc"},
                                                     "synapse_model": "GABA_A_syn",
                                                     "mask": {"circular": {"radius": 12.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.15, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0}), # original
                                                     "weights": 1.0*weight_gain}),
                           ('Rs_layer', 'Ts_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_exc"},
                                                     "synapse_model": "GABA_B_syn",
                                                     "mask": {"circular": {"radius": 12.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.05, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0}), #original
                                                     "weights": 1.0*weight_gain}),
                           ('Rs_layer', 'Ts_layer',{ "sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_inh"},
                                                     "synapse_model": "GABA_A_syn",
                                                     "mask": {"circular": {"radius": 12.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.15, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0}), # original
                                                     "weights": 1.0*weight_gain}),
                           ('Rs_layer', 'Ts_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Tp_inh"},
                                                     "synapse_model": "GABA_B_syn",
                                                     "mask": {"circular": {"radius": 12.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.05, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0}), #original
                                                     "weights": 1.0*weight_gain}),
                           ('Rs_layer', 'Rs_layer', {"sources": {"model": "Rp"},
                                                     "targets": {"model": "Rp"},
                                                     "synapse_model": "GABA_B_syn",
                                                     "mask": {"circular": {"radius": 12.0 * dpcS}},
                                                     "kernel": {"gaussian": {"p_center": 0.5, "sigma": 7.5 * dpcS}},
                                                     #"weights": 1.0})]: #original
                                                     "weights": 1.0*weight_gain})]:
        thal = thal_base.copy()
        thal.update(conn)
        allconns.append([src,tgt,thal])


    # Thalamocortical connections
    Vp_Thalamocortical_base = {
        "connection_type": "convergent",
        "sources": {"model": "Tp_exc"},
        "synapse_model": "AMPA_syn",
        #"weights": 5.0, #original
        "weights": 5.0 * weight_gain,
        #"weights": 6.5 * weight_gain,
        "delays": {"uniform": {"min": 2.75, "max": 3.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    # Horizontally tuned
    Vp_Thalamocortical_base.update({"mask": {"rectangular": {"lower_left" : [-4.0*dpcP, -1.0*dpcP],
                                                             "upper_right": [ 4.0*dpcP,  1.0*dpcP]}}})
    for conn in [{"targets": {"model": "L4_exc" }, "kernel": 0.5},
                 {"targets": {"model": "L56_exc"}, "kernel": 0.3}]:
        Vp_Thalamocortical_base.update(conn)
        allconns.append(['Tp_layer','Vp_horizontal', Vp_Thalamocortical_base.copy()])

    # Vertically tuned
    Vp_Thalamocortical_base.update({"mask": {"rectangular": {"lower_left" : [-1.0*dpcP, -4.0*dpcP],
                                                             "upper_right": [ 1.0*dpcP,  4.0*dpcP]}}})
    for conn in [{"targets": {"model": "L4_exc" }, "kernel": 0.5},
                 {"targets": {"model": "L56_exc"}, "kernel": 0.3}]:
        Vp_Thalamocortical_base.update(conn)
        allconns.append(['Tp_layer','Vp_vertical', Vp_Thalamocortical_base.copy()])

    #! Diffuse connections
    Vp_Thalamocortical_diff_base = {
        "connection_type": "convergent",
        "sources": {"model": "Tp_exc"},
        "synapse_model": "AMPA_syn",
        #"weights": 5.0,  # original
        "weights": 5.0 * weight_gain,
        #"weights": 6.5 * weight_gain,
        "mask": {"circular": {"radius": 5.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.1, "sigma": 7.5 * dpcP}},
        "delays": {"uniform": {"min": 2.75, "max": 3.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    for conn in [{"targets": {"model": "L4_inh" }},
                 {"targets": {"model": "L56_inh"}}]:
        Vp_Thalamocortical_diff_base.update(conn)
        allconns.append(['Tp_layer','Vp_horizontal', Vp_Thalamocortical_diff_base.copy()])
        allconns.append(['Tp_layer','Vp_vertical', Vp_Thalamocortical_diff_base.copy()])

    # thalamocortical connections: inferred from text p 1677 and Lumer
    #! Diffuse connections
    TsVsBase = {
        "connection_type": "divergent",
        "sources": {"model": "Tp_exc"},
        "synapse_model": "AMPA_syn",
        #"weights": 5.0, # original
        "weights": 5.0 * weight_gain,
        #"weights": 2.5 * weight_gain,
        "mask": {"rectangular": {"lower_left" : [-2.0* dpcS, -2.0* dpcS],"upper_right": [ 2.0* dpcS,  2.0* dpcS]}},
        "delays": {"uniform": {"min": 2.75, "max": 3.25}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    for c in [({"targets": {"model": "L4_exc"}}, 0.5),
              ({"targets": {"model": "L4_inh" }}, 0.5),
              ({"targets": {"model": "L56_exc"}}, 0.3),
              ({"targets": {"model": "L56_inh" }}, 0.3)]:
        cdict = c[0]
        c[0]['kernel'] =  c[1]
        c[0].update(TsVsBase)
        allconns += [ ('Ts_layer', t, c[0] ) for t in ['Vs_horizontal', 'Vs_cross', 'Vs_vertical'] ]


    # Corticothalamic connections
    ctConnections = []

    Vp_Corticothalamic_base = {
        "connection_type": "divergent",
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 5.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.5, "sigma": 7.5 * dpcP}},
        #"weights": 1.0, #original
        "weights": 1.0 * weight_gain,
        "delays": {"uniform": {"min": 7.5, "max": 8.5}},
        "allow_autapses": False,
        "allow_multapses": False
    }
    for conn in [{"sources": {"model": "L56_exc"}, "targets": {"model": "Tp_exc" }},
                 {"sources": {"model": "L56_exc"}, "targets": {"model": "Tp_inh" }}]:
        ndict = Vp_Corticothalamic_base.copy()
        ndict.update(conn)
        ctConnections.append(ndict)

    corRet = Vp_Corticothalamic_base.copy()
    #corRet.update({"sources": {"model": "L56_exc"}, "targets": {"model": "Rp"}, "weights": 2.5}) # original
    corRet.update({"sources": {"model": "L56_exc"}, "targets": {"model": "Rp"}, "weights": 2.5*weight_gain})
    #corRet.update({"sources": {"model": "L56_exc"}, "targets": {"model": "Rp"}, "weights": 5.0})  # keiko
    #corRet.update({"sources": {"model": "L56_exc"}, "targets": {"model": "Rp"}, "weights": 5.0*weight_gain})  # keiko

    [allconns.append(['Vp_horizontal','Tp_layer',c]) for c in ctConnections]
    [allconns.append(['Vp_vertical','Tp_layer',c]) for c in ctConnections]

    [allconns.append(['Vp_horizontal','Rp_layer',c]) for c in [corRet]]
    [allconns.append(['Vp_vertical','Rp_layer',c]) for c in [corRet]]


    ctConnections = []

    Vs_Corticothalamic_base = {
        "connection_type":"divergent",
        "sources": {"model": "L56_exc"},
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 5.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.5, "sigma": 7.5 * dpcP}},
        #"weights": 1.0, # original
        "weights": 1.0*weight_gain,
        "delays": {"uniform": {"min": 7.5, "max": 8.5}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    for conn in [{"targets": {"model": "Tp_exc" }},
                 {"targets": {"model": "Tp_inh" }}]:
        ndict = Vs_Corticothalamic_base.copy()
        ndict.update(conn)
        ctConnections.append(ndict)

    corRet = Vs_Corticothalamic_base.copy()
    #corRet.update({"targets": {"model": "Rp"}, "weights": 2.5}) # original
    corRet.update({"targets": {"model": "Rp"}, "weights": 2.5*weight_gain})

    [allconns.append(['Vs_horizontal','Ts_layer',c]) for c in ctConnections]
    [allconns.append(['Vs_cross','Ts_layer',c]) for c in ctConnections]
    [allconns.append(['Vs_vertical','Ts_layer',c]) for c in ctConnections]

    [allconns.append(['Vs_horizontal','Rs_layer',c]) for c in [corRet]]
    [allconns.append(['Vs_cross','Rs_layer',c]) for c in [corRet]]
    [allconns.append(['Vs_vertical','Rs_layer',c]) for c in [corRet]]



    Vp_Ts_Corticothalamic_base = {
        "connection_type":"divergent",
        "sources": {"model": "L56_exc"},
        "synapse_model": "AMPA_syn",
        "mask": {"circular": {"radius": 2.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 1.0, "sigma": 7.5 * dpcP}},
        #"weights": 5.0,  # original
        "weights": 5.0 * weight_gain,  #
        "delays": {"uniform": {"min": 7.5, "max": 8.5}},
        "allow_autapses": False,
        "allow_multapses": False
    }

    VpTsConns = [('Vp_vertical', {'targets': {'model': 'Tp_exc'}}),
                 ('Vp_vertical', {'targets': {'model': 'Tp_inh'}}),
                 ('Vp_horizontal', {'targets': {'model': 'Tp_exc'}}),
                 ('Vp_horizontal', {'targets': {'model': 'Tp_inh'}})]

    allconns += [ (c[0], 'Ts_layer', updateDicts(Vp_Ts_Corticothalamic_base, c[1]) ) for c in VpTsConns ]


    # --------------------------------------------------------------------#
    # ---------- THALAMIC INPUT ----------------------------------------- #
    # --------------------------------------------------------------------#
    nest.CopyModel('static_synapse', 'reticular_projection',
                   params={'receptor_type':receptors['AMPA']})
    retThal = {
        "connection_type": "divergent",
        "synapse_model": "reticular_projection",
        "mask": {"circular": {"radius": 1.0 * dpcP}},
        "kernel": {"gaussian": {"p_center": 0.75, "sigma": 2.5 * dpcP}},
        #"weights": 10.0,  # original
        "weights": 10.0*weight_gain,  # original
        #"weights": 7.5 * weight_gain,
        "delays": 1.0,
        "allow_autapses": False,
        "allow_multapses": False
    }
    for conn in [{"targets": {"model": "Tp_exc"}},
                 {"targets": {"model": "Tp_inh"}}]:
        retThal.update(conn)
        allconns.append(['Retina_layer','Tp_layer', retThal.copy()])

    '''
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

    ridx_vertical = np.random.randint(num_neuron, size=(1,num_ib))[0].tolist()
    ridx_horizontal = np.random.randint(num_neuron, size=(1,num_ib))

    for i in range(1,num_ib,1):
        nest.SetStatus([L56_vertical_idx[ridx_vertical[i]]], {'h_g_peak': 1.0})
        nest.SetStatus([L56_horizontal_idx[ridx_horizontal[i]]], {'h_g_peak': 1.0})
    '''

    return allconns

