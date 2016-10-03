

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

root_dir = '/home/leonardo/projects/nsdm/hill_tononi_synthesis/data'

files_to_load = '/z_*_Vp*.pickle'

all_files = glob.glob(root_dir + files_to_load)

# load epochs and calculate
def processInput(idx, next_file):

    t = time.time()

    print ('Calculating PCI for' + next_file)
    with open(next_file, 'r') as ft:
        data = pickle.load(ft)

    data_to_pci = data[::10,::10]
    #data_to_pci = data
    data_c = pci.pci(data_to_pci)

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


num_cores = 16
results = Parallel(n_jobs=num_cores, verbose=100)(delayed(processInput)(fileid, next_file) for fileid, next_file in enumerate(all_files))

output_file = 'results_pci_' + files_to_load.replace('*', '').replace('/', '')

with open(root_dir + '/' + output_file, 'w') as f:
    pickle.dump(results, f)

print('Done PCI!')

log = open(root_dir + '/pci_log.txt', 'a+')

log.write('\n\n')
log.write("End: " + time.strftime("%c"))
log.write('\n\n')
log.close()

