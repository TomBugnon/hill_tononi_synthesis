import numpy as np
import pickle

# only for keiko's environment
import sys
sys.path.append('/usr/lib/python2.7/dist-packages')

import matplotlib.pyplot as plt


filename = 'spikes_Vp_h L4pyr.pickle'
#filename = 'spikes_Vp_v L4pyr.pickle'

#--- random
dir_to_load = 'sim_1_phi_dg_0.00_retAC_30.00_simtime_40.00_lambda_dg_-1.00_sim_interval_2.00_visSize_8.00_retDC_30.00_f_dg_20.00_N_40.00_detectors/'
title_str = 'Histgram of ISI (Random ' + filename + ')'

#--- structured
#dir_to_load = 'sim_1_phi_dg_0.00_retAC_30.00_simtime_40.00_lambda_dg_2.00_sim_interval_2.00_visSize_8.00_retDC_30.00_f_dg_20.00_N_40.00_detectors/'
#title_str = 'Histgram of ISI (Structured ' + filename + ')'

root_dir = '/home/kfujii2/newNEST2/iaf_model/data/'
data = pickle.load(open(root_dir+dir_to_load+filename,'r'))

senders = data['senders']
times = data['times']

min_neuron_idx = np.min(senders)
max_neuron_idx = np.max(senders)

ftime = []
isi = []
isi_all = [] # to draw a histogram of ISI
for i in range(min_neuron_idx,(max_neuron_idx+1),1):

    # Find positions for i-th neuron
    pos = np.where(senders==i)

    # Time stamp when i-th neuron fired
    tmp_ftime = times[pos]
    tmp_isi = tmp_ftime[1:len(tmp_ftime)] - tmp_ftime[0:(len(tmp_ftime)-1)]

    # Connect data
    ftime.append(tmp_ftime.tolist())
    isi.append(tmp_isi.tolist())
    isi_all = isi_all + tmp_isi.tolist()


#plt.hist(isi_all, bins=range(0,30,1))
#plt.title(title_str)
#plt.show()


# Segregate states based on ISI
# --- 1: Fired but not bursting
# --- 2: Bursting
# First, create ones matrix ( size=len(times)+len(isi[i]+1) ) <- Set 'fired but not bursting' value as a default
# Second, if ISI is smaller than threshold, change values to 'bursting_val'

threshold = 3 # if ISI is smaller than 3 msec, the neuron is bursting
bursting_val = 2

num_neurons = max_neuron_idx - min_neuron_idx + 1
firing_segregation = []

for i in range(0,num_neurons,1):

    tmp_state = np.ones( (len(isi[i])+1) )

    for t in range(0,len(isi[i]),1):

        if isi[i][t]<threshold:
            tmp_state[t] = bursting_val
            tmp_state[t+1] = bursting_val

    firing_segregation.append(tmp_state.tolist())


# convert time to step
# Because 'times' is calculated in 0.1msec resolution (check nest.GetKernelStatus()['resolution']),
#  to get the state matrix, do int(ftime*10) 10=1/resolution

import nest
resolution = nest.GetKernelStatus()['resolution']
sim_time = 40

# 1 msec resolution version
num_step = sim_time
states = np.zeros((num_neurons,num_step))
states_all = [] # to draw a histogram of state labeling
for i in range(0,num_neurons,1):

    tmp_step = (np.asarray(ftime[i])+0.5).astype(int)  # rounding, ex) 0.5-1.4msec -> 1msec
    states[i][tmp_step] = firing_segregation[i]
    states_all = states_all + states[i].astype(int).tolist()

'''
# 0.1 mec resolution version
num_step = int( sim_time * (1/resolution) )
states = np.zeros((num_neurons,num_step))
states_all = [] # to draw a histogram of state labeling
for i in range(0,num_neurons,1):

    tmp_step = ((1/resolution)*np.asarray(ftime[i])).astype(int)
    states[i][tmp_step] = firing_segregation[i]
    states_all = states_all + states[i].astype(int).tolist()
'''


# count the number of states (i.e. D1)
unique_states = np.vstack({tuple(row) for row in states})
d1 = len(unique_states)

#plt.hist(states_all)
#plt.title('states: 1 = Regular spiking, 2 = Bursting')
#plt.show()

print('end')




