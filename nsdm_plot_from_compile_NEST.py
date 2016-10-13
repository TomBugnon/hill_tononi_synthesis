

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

from scipy import stats

def update_line(num, data, line):
    line.set_data(data[..., :num])
    return line,

#

root_dir = '/home/leonardo/projects/nsdm/hill_tononi_synthesis/data'

#files_to_load = '/spikes_Vp*L4*.pickle'
#files_to_load = '/results_D16*.pickle'
files_to_load = '/results_D0*VpL4.pickle'
P = 2

#files_to_load = '/results_D0*Retina.pickle'
#P = 1

all_files = glob.glob(root_dir + files_to_load)

# Set up formatting for the movie files
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


for idx, next_file in enumerate(all_files):

    figure_name = path.split(next_file)[-1].replace('.pickle', '')
    figure_folder = root_dir + '/' + figure_name + '_plot'
    if not os.path.isdir(figure_folder):
        os.makedirs(figure_folder)

    movie_folder = figure_folder + '/animation'
    if not os.path.isdir(movie_folder):
        os.makedirs(movie_folder)

    print ('Loading ' + next_file)
    with open(next_file, 'r') as f:
        results = pickle.load(f)

    #print data.items()

    for data_idx, this_data in enumerate(results['z']):
        fig = plt.figure()
        #plt.set_cmap('gray')

        this_data = np.absolute(this_data)
        #thresh = 0.00001
        #thresh = 1
        thresh = stats.t.isf(0.05, results['N'])
        this_data[np.isnan(this_data)] = 0

        fig = plt.figure()
        plt.plot(this_data.flatten(), 'bs')
#        axes = plt.gca()
#        axes.set_ylim([0,10])
        fig.savefig(figure_folder + '/' + figure_name + '_z_dotplot_%d.png' % data_idx )
        plt.clf()

        #downsample = 20
        #this_data = np.array([1 * (np.any(this_data[x:x + downsample-1, :], 0)) for x in range(0, this_data.shape[0], downsample)])

        plot_data = this_data
        cmap = matplotlib.cm.hot
        plt.imshow(plot_data, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 3*thresh, vmin=0)
        #plt.imshow(this_data, aspect='auto')
        plt.colorbar()
        fig.savefig(figure_folder + '/' + figure_name + '_z_not_thresh_%d.png' % data_idx )
        plt.clf()

        cbaron = False
        D = np.sqrt(this_data.shape[0]/P)
        if D%1 == 0.0:
            D = np.int64(D)
            #if np.mod(this_data.shape[0],2) == 0:
            #D = np.int64(np.round(np.sqrt(this_data.shape[0]/P)))
            N = this_data.shape[0]
            T = this_data.shape[1]
            _all_data = np.zeros((T,P,D,D))
            print('saving figures for animation')
            for t in range(0,T):
                for p in range(1,P+1,1):
                    p1 = (p-1)*N/P
                    p2 = p1+(N/P)
                    this_screen = this_data[p1:p2,t].reshape(D,D)
                    _all_data[t,p-1,:,:] = this_screen
                    #plt.imshow(this_screen, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 3*thresh, vmin=0)
                    #if not cbaron:
                    #    plt.colorbar()
                    #    cbaron = True
                    #fig.savefig(movie_folder + '/' + figure_name + '_p_%d_vh_%d_t_%d.png' % (p, data_idx, t) )

            for p in range(1,P+1,1):
                mean_screen = np.mean(_all_data[:,p-1,:,:],0)
                cmap = matplotlib.cm.jet
                plt.clf()
                plt.imshow(mean_screen, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 0.4, vmin=0.2)
                plt.colorbar()
                #plt.imshow(mean_screen, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 3*thresh, vmin=0)
                fig.savefig(figure_folder + '/' + figure_name + '_screen_z_p_%d_vh_%d.png' % (p, data_idx) )

        this_data[this_data < thresh] = 0

        plt.clf()
        plot_data = this_data
        masked_array = np.ma.array (plot_data, mask=np.isnan(this_data))
        cmap = matplotlib.cm.gray
        cmap.set_bad('black',1.)
        plt.imshow(masked_array, interpolation='nearest', aspect='auto', cmap=cmap)
        #plt.imshow(this_data, aspect='auto')
        plt.colorbar()
        fig.savefig(figure_folder + '/' + figure_name + '_z_%d.png' % data_idx )
        plt.clf()

       # save for PCI
        with open((root_dir + '/z_' + figure_name + '_%d.pickle') % data_idx, 'w') as fz:
            z_data = plot_data
            pickle.dump(z_data, fz)



    for data_idx, this_data in enumerate(results['any']):
        fig = plt.figure()
        #plt.set_cmap('gray')

        #downsample = 20
        #this_data = np.array([1 * (np.any(this_data[x:x + downsample-1, :], 0)) for x in range(0, this_data.shape[0], downsample)])

        plot_data = this_data
        cmap = matplotlib.cm.gray
        plt.imshow(plot_data, interpolation='nearest', aspect='auto', cmap=cmap, vmax=1, vmin=0)
        #plt.imshow(this_data, aspect='auto')
        plt.colorbar()
        fig.savefig(figure_folder + '/' + figure_name + '_any_%d.png' % data_idx )
        plt.clf()

        #cbaron = False
        #if np.mod(this_data.shape[0],2) == 0:
        #    D = np.int64(np.round(np.sqrt(this_data.shape[0]/P)))
        #    N = this_data.shape[0]
        #    print('saving figures for animation')
        #    for t in range(0,this_data.shape[1]-1):
        #        for p in range(1,P+1,1):
        #            p1 = (p-1)*N/P
        #            p2 = p1+(N/P)
        #            this_screen = this_data[p1:p2,t].reshape(D,D)
        #            plt.imshow(this_screen, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 1, vmin=0)
        #            if not cbaron:
        #                plt.colorbar()
        #                cbaron = True
        #            fig.savefig(movie_folder + '/' + figure_name + '_p_%d_vh_%d_t_%d.png' % (p, data_idx, t) )


        # save for PCI
        with open((root_dir + '/any_' + figure_name + '_%d.pickle') % data_idx, 'w') as fz:
            z_data = plot_data
            pickle.dump(z_data, fz)


    for data_idx, this_data in enumerate(results['mean']):
        fig = plt.figure()
        plt.set_cmap('gray')
        plt.imshow(this_data, aspect='auto')
        plt.colorbar()
        fig.savefig(figure_folder + '/' + figure_name + '_mean_%d.png' % data_idx )
        plt.clf()

        # If spatial coarse graining in the compile script left a square, plot/save the screen
        D = np.sqrt(this_data.shape[0]/P)
        if D%1 == 0.0:
            D = np.int64(D)
            N = this_data.shape[0]
            print('saving figures for animation')
            for t in range(0,this_data.shape[1]-1):
                for p in range(1,P+1,1):
                    p1 = (p-1)*N/P
                    p2 = p1+(N/P)
                    this_screen = this_data[p1:p2,t].reshape(D,D)
                    _all_data[t,p-1,:,:] = this_screen
                    #plt.imshow(this_screen > 0, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 1, vmin=0)
                    #fig.savefig(movie_folder + '/' + figure_name + '_p_%d_vh_%d_t_%d.png' % (p, data_idx, t) )

            for p in range(1,P+1,1):
                mean_screen = np.mean(_all_data[:,p-1,:,:],0)
                plt.clf()
                plt.imshow(mean_screen, interpolation='nearest', aspect='auto', cmap=cmap)
                plt.colorbar()
                #plt.imshow(mean_screen, interpolation='nearest', aspect='auto', cmap=cmap, vmax = 1, vmin=0)
                fig.savefig(figure_folder + '/' + figure_name + '_screen_mean_p_%d_vh_%d.png' % (p, data_idx) )

    #ims = []
    #for t in range(0,this_neurons.shape[1]-1):
    #    this_screen = this_neurons[:,t].reshape(40,40)
    #    ims.append((plt.imshow(this_screen, aspect='auto'),))

    #im_ani = animation.ArtistAnimation(fig, ims, interval=50, repeat_delay=3000,
    #                               blit=True)
    #im_ani.save('test_st.mp4', writer=writer)



