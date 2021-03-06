import nest
import nest.topology as topo
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import nsdm_run_params as rp

# ====================
#  Load an image
# ====================
# Load an image, and convert from RGB to gray
test_img_tmp = Image.open(rp.input_figure).convert('L')





# Convert PIL format to numpy format
# Convert int to float
# Reshape to a row vector
test_img = np.array(test_img_tmp).astype(float)
N = len(test_img)
input_img = test_img.reshape(1, N*N)[0]


# ====================
#  Set nest kernel
# ====================
nest.ResetKernel()
nest.SetKernelStatus({"local_num_threads": 16})
nest.SetStatus([0],{'print_time': True})


# ======================
#  Create retina layer
# ======================
nest.CopyModel('ht_neuron', 'RetinaNode',
               params={"Theta_eq": -51.0,
                       "Tau_theta": 2.0,
                       "spike_duration": 2.0,
                       "Tau_spike": 1.75,
                       "Tau_m": 16.0} )

layerProps = {'rows': N,
              'columns': N,
              'elements': 'RetinaNode'}
retina = topo.CreateLayer(layerProps)


# ==============================
#  Connect recorder & detector
# ==============================
recorder = nest.Create('multimeter',
                      params={'interval': 0.1,
                              'record_from': ['V_m', 'spike_input_AMPA']})

detector = nest.Create('spike_detector',
                       params={"withgid": True, "withtime": True})

tgts = nest.GetLeaves(retina)[0]
nest.Connect(recorder, tgts)
#nest.Connect(tgts, detector)


# ===========================
#  Connect spike generators
# ===========================
#dcg =  nest.Create('dc_generator')

# Create spike_generator for all retina neurons
sg = nest.Create('spike_generator', N*N)
nest.Connect(sg, detector)

# Connect using AMPA synapses
receptors = nest.GetDefaults('ht_neuron')['receptor_types']
w = 20.0
for i in range(0, N*N, 1):
    nest.Connect([sg[i]], [tgts[i]], syn_spec={'receptor_type': receptors['AMPA'], 'weight': w})


# ==============================
#  Simulation
# ==============================

# if the intensity of pixel is smaller than 'intensity_threshold',
# set spike_generator ON
# In this simulation, retina receives spikes only when t = 100 ms
# TODO check spike timing (=nest.GetStatus(detector)[0]['events']['times']). This should be 101 (but 105.5 and 112.6)

t_end = 500 # int
sim_interval = 1 # int
intensity_threshold = 200.0

# Set duration (float value)
# Without duration, the network didn't have any spike
interval = np.int(t_end/sim_interval)
l = np.random.poisson(100, N*N)
for i in range(0, N*N, 1):
   if input_img[i] < intensity_threshold:
       t = np.random.randint(1, t_end, l[i])
       t.sort()
       nest.SetStatus([sg[i]], {'spike_times': t.astype(float).tolist()})

nest.Simulate(t_end)

spikes = nest.GetStatus(detector, "events")[0]
screen_time = np.zeros((1600,500))
screen_time[spikes['senders']-min(sg), spikes['times'].astype(np.int)] = 1
screen = np.mean(screen_time, 1)
s = screen.reshape((40,40))
plt.imshow(s)
plt.colorbar()

print('end')
