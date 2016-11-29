

if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

import numpy as np
import glob
import pickle

# git clone https://github.com/noreun/pypci
from pypci import pci

import time
from joblib import Parallel, delayed

#downsample_neuron = 0
#downsample_neuron = 4
downsample_neuron = 1000

#downsample_time = 0
#downsample_time = 4
downsample_time = 100


#root_dir = '/home/leonardo/projects/nsdm/data'
root_dir = '/home/leonardo/projects/nsdm/data/from_allen'

files_to_load = '/mov_2_scrbl_normal_*.pickle'
#files_to_load = '/z_*_Vp*.pickle'
#files_to_load = '/any_*_Vp*.pickle'

all_files = glob.glob(root_dir + files_to_load)

# load epochs and calculate
def processInput(idx, next_file):

    t = time.time()

    print ('Calculating PCI for' + next_file)
    with open(next_file, 'r') as ft:
        data = pickle.load(ft)

    data_to_pci = data

    #data_to_pci = data[::10,::10]
    #data_to_pci = np.array([1 * (np.any(data[x:x + 9, :], 0)) for x in range(0, data.shape[0], 10)])

    # if downsample_neuron > 0:
    #     data = np.array([1 * (np.mean(data[x:x + downsample_neuron-1, :], 0) > .5) for x in range(0, data.shape[0], downsample_neuron)])

    # if downsample_time > 0:
    #    data = np.array([1 * (np.mean(data[:, x:x + downsample_time-1], 0) > .5) for x in range(0, data.shape[1], downsample_time)])

    # quick speed test
    #data_dsampled = data[10000:16000:10,501:]; data_to_pci = data_dsampled; print(data_to_pci.shape)
    #t = time.time(); data_c = pci.pci(data_to_pci); elapsed = time.time() - t; print("Elapsed time : %.2f s" % elapsed)

    #print('only calculating LIF Neurons (10000 to 20000) from 500 ms to 3000 ms')
    #data_dsampled = data[10000:20001,500:3000]

    print('only calculating bio-physically realistic neurons (1 to 10000) from 500 ms to 5000 ms')
    data_dsampled = data[:10000,500:3000]

    data_to_pci = data_dsampled

    data_c = pci.calculate(data_to_pci)

    elapsed = time.time() - t

    log = open(root_dir + '/pci_log.txt', 'a+')

    m = "Complexity of data for %s : %d" % (next_file, data_c)
    print(m)
    log.write(m+'\n')

    m = "Elapsed time : %.2f s" % elapsed
    print(m)
    log.write(m+'\n')

    log.close()
    return data_c

log = open(root_dir + '/pci_log.txt', 'a+')

log.write('\n\n')
log.write("Start: " + time.strftime("%c"))
log.write('\n\n')
log.close()

num_cores = 16
results = Parallel(n_jobs=num_cores, verbose=100)(delayed(processInput)(fileid, next_file) for fileid, next_file in enumerate(all_files))

#results = []
#results.append(processInput(0, all_files[0]))

output_file = 'pci_' + files_to_load.replace('*', '').replace('/', '')

with open(root_dir + '/' + output_file, 'w') as f:
    pickle.dump(results, f)

print('Done PCI!')

log = open(root_dir + '/pci_log.txt', 'a+')

log.write('\n\n')
log.write("End: " + time.strftime("%c"))
log.write('\n\n')
log.close()

