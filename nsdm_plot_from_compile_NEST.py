

if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

import os
import numpy as np
import glob
import pickle
from pypci import pci
import time

import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update_line(num, data, line):
    line.set_data(data[..., :num])
    return line,

#

root_dir = '/home/leonardo/projects/nsdm/hill_tononi_synthesis/data'

#files_to_load = '/spikes_Vp*L4*.pickle'
files_to_load = '/*.pickle'

all_files = glob.glob(root_dir + files_to_load)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


for idx, next_file in enumerate(all_files):

    figure_name = path.split(next_file)[-1].replace('.pickle', '')
    figure_folder = root_dir + '/' + figure_name + '_plot'
    if not os.path.isdir(figure_folder):
        os.makedirs(figure_folder)

    print ('Loading ' + next_file)
    with open(next_file, 'r') as f:
        results = pickle.load(f)

    #print data.items()

    for data_idx, this_data in enumerate(results['z']):
        fig = plt.figure()
        #plt.set_cmap('gray')

        thresh = .1
        this_data[np.isnan(this_data)] = 0
        this_data[this_data < thresh] = 0

        plot_data = this_data
        masked_array = np.ma.array (plot_data, mask=np.isnan(this_data))
        cmap = matplotlib.cm.gray
        cmap.set_bad('black',1.)
        plt.imshow(masked_array, interpolation='nearest', aspect='auto', cmap=cmap, vmax = thresh)
        #plt.imshow(this_data, aspect='auto')
        plt.colorbar()
        fig.savefig(figure_folder + '/' + figure_name + '_z_%d.png' % data_idx )
        plt.clf()

        fig = plt.figure()
        plt.plot(plot_data.flatten(), 'bs')
        fig.savefig(figure_folder + '/' + figure_name + '_z_dotplot_%d.png' % data_idx )
        plt.clf()

        # save for PCI
        with open((root_dir + '/z_' + figure_name + '_%d.pickle') % data_idx, 'w') as fz:
            z_data = plot_data
            pickle.dump(z_data, fz)


    for data_idx, this_data in enumerate(results['mean']):
        fig = plt.figure()
        plt.set_cmap('gray')
        plt.imshow(this_data, aspect='auto')
        plt.colorbar()
        fig.savefig(figure_folder + '/' + figure_name + '_mean_%d.png' % data_idx )
        plt.clf()

    #ims = []
    #for t in range(0,this_neurons.shape[1]-1):
    #    this_screen = this_neurons[:,t].reshape(40,40)
    #    ims.append((plt.imshow(this_screen, aspect='auto'),))

    #im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
    #                               blit=True)
    #im_ani.save('test_st.mp4', writer=writer)

    # for t in range(0,this_neurons.shape[1]-1):
    #     this_screen = this_neurons[:,t].reshape(40,40)
    #     plt.imshow(1 * (this_screen > 0), aspect='auto')
    #     fig.savefig(figure_folder + '/' + figure_name + '%d.png' % t )


