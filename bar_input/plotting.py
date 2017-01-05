#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Plot functions
# Author: Pablo Martinez Cañada (pablomc@ugr.es)

import nest
import nest.topology as tp
import numpy as np
import matplotlib.pyplot as plt
import os
import pylab as pl
from sys import stdout
from time import sleep

## A: membrane potential rasters

def potential_raster(fig,recorders,recorded_models,starting_neuron,number_cells,simtime,resolution,rows,cols,starting_pos):

    pos = starting_pos
    rec_from = ["rate","V_m"]

    for population, model in recorded_models:

            l = [nd for nd in np.arange(0,len(recorders)) if (recorders[nd][1] == population and recorders[nd][2] == model)][0]
            data = nest.GetStatus(recorders[l][0],keys='events')
            pop = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]

            # Retina layer always the first (l=0)
            if l==0:
                    rf = rec_from[0]
            else:
                    rf = rec_from[1]

            raster = np.zeros((number_cells, round((simtime - resolution)/resolution)))
            senders = data[0]['senders']

            for neuron in range(0, len(raster)):
                    selected_senders = np.where(senders==pop[neuron + starting_neuron])
                    raster[neuron,:] = (data[0][rf])[selected_senders[0]]

            Vax = plt.subplot2grid((rows,cols), (pos,0), colspan=cols)


            if l>0:
                cax = Vax.matshow(raster,aspect='auto',vmin=-70.0,vmax=-45.0)
                Vax.axes.get_xaxis().set_ticks([])
            else:
                # cax = Vax.matshow(raster,aspect='auto') # original
                cax = Vax.matshow(raster, aspect='auto',vmin=0.0,vmax=200.0) # keiko
                fig.colorbar(cax,ticks=[0.0, 100.0, 200.0],orientation='horizontal') # keiko

            Vax.xaxis.tick_top()
            plt.setp(Vax, yticks=[0, number_cells], yticklabels=['0', str(number_cells-1)])
            Vax.set_ylabel(model)
            pos+=1


    fig.colorbar(cax,ticks=[-70.0, -57.5, -45.0],orientation='horizontal')


## B: intracellular potentials

#def intracellular_potentials(fig, recorders, recorded_models, starting_neuron, rows, cols, starting_pos): #original
def intracellular_potentials(fig, recorders, recorded_models, starting_neuron, rows, cols, starting_pos,total_time): #keiko

    pos = starting_pos
    counter = 0

    for population, model in recorded_models:

            l = [nd for nd in np.arange(0,len(recorders)) if (recorders[nd][1] == population and recorders[nd][2] == model)][0]
            data = nest.GetStatus(recorders[l][0],keys='events')
            pop = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]

            senders = data[0]['senders']
            selected_senders = np.where(senders==pop[starting_neuron])

            Vax = plt.subplot2grid((rows,cols), (pos,0), colspan=cols)
            Vax.plot( (data[0]['V_m'])[selected_senders[0]] )
            if(counter!=(len(recorded_models)-1)):
                Vax.axes.get_xaxis().set_ticks([])
            plt.setp(Vax, yticks=[-80, 0], yticklabels=['-80', '0'])
            Vax.set_ylabel(model)
            Vax.set_xlabel('time (ms)')
            pos+=1
            counter+=1

            # keiko
            #min_x = np.min(data[0]['times'])
            #max_x = np.max(data[0]['times'])
            plt.xlim(0, total_time)


## C: time-averaged topographic representation of the membrane potential

def topographic_representation(fig,recorders,recorded_models,labels,number_cells,simtime,resolution,rows,cols,start,stop,starting_pos,col_paint):

    col_ind = col_paint
    pos = starting_pos
    counter = 0

    for population, model in recorded_models:

            l = [nd for nd in np.arange(0,len(recorders)) if (recorders[nd][1] == population and recorders[nd][2] == model)][0]
            data = nest.GetStatus(recorders[l][0],keys='events')
            pop = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]

            raster = np.zeros((number_cells, number_cells))
            mult_pos = 0
            senders = data[0]['senders']
            for x in range(0, number_cells):
                    for y in range(0, number_cells):

                            selected_senders = np.where(senders==pop[mult_pos])
                            ind_rec = (data[0]['V_m'])[selected_senders[0]]
                            raster[y,x] = np.sum( ind_rec[int(start/resolution):int(stop/resolution)] ) / ((stop-start)/resolution)
                            mult_pos+=1


            Iax = plt.subplot2grid((rows,cols), (pos,col_ind), colspan=1)
            cax2 = Iax.matshow(raster,aspect='auto',vmin=-70.0,vmax=-45.0)
            Iax.xaxis.tick_bottom()
            Iax.axes.get_xaxis().set_ticks([])
            Iax.axes.get_yaxis().set_ticks([])
            Iax.set_ylabel(model)
            Iax.set_xlabel(labels[counter])
            counter+=1
            col_ind+=1

    fig.colorbar(cax2,ticks=[-70.0, -57.5, -45.0])


## D: Intrinsic currents

def intrinsic_currents(recorders,recorded_models,starting_neuron):

    fig = plt.figure()
    counter = 0

    for population, model in recorded_models:

            l = [nd for nd in np.arange(0,len(recorders)) if (recorders[nd][1] == population and recorders[nd][2] == model)][0]
            data = nest.GetStatus(recorders[l][0],keys='events')
            pop = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]

            senders = data[0]['senders']
            selected_senders = np.where(senders==pop[starting_neuron])

            Vax = plt.subplot2grid((1,1), (0,0), colspan=1)
            Vax.plot( (data[0]['I_h'])[selected_senders[0]] ,label='I_h')
            Vax.plot( (data[0]['I_KNa'])[selected_senders[0]] ,label='I_KNa')
            Vax.plot( (data[0]['I_NaP'])[selected_senders[0]],label='I_NaP' )

            Vax.legend(fancybox=True)
            Vax.set_ylabel(model)
            Vax.set_xlabel('time (ms)')


# Keiko
def synaptic_currents(recorders, recorded_models, starting_neuron):

    fig = plt.figure()
    counter = 0

    for population, model in recorded_models:
        l = [nd for nd in np.arange(0, len(recorders)) if
             (recorders[nd][1] == population and recorders[nd][2] == model)][0]
        data = nest.GetStatus(recorders[l][0], keys='events')
        pop = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0] == model]

        senders = data[0]['senders']
        selected_senders = np.where(senders == pop[starting_neuron])

        Vax = plt.subplot2grid((1, 1), (0, 0), colspan=1)
        Vax.plot((data[0]['I_syn_AMPA'])[selected_senders[0]], label='AMPA')
        Vax.plot((data[0]['I_syn_NMDA'])[selected_senders[0]], label='NMDA')
        Vax.plot((data[0]['I_syn_GABA_A'])[selected_senders[0]], label='GABA_A')
        Vax.plot((data[0]['I_syn_GABA_B'])[selected_senders[0]], label='GABA_B')

        Vax.legend(fancybox=True)
        Vax.set_ylabel(model)
        Vax.set_xlabel('time (ms)')


## E: movie

def makeMovie(fig,recorders,recorded_models,labels,number_cells,simtime,resolution):

    counter = 0

    print "creating movies..."

    for population, model in recorded_models:

        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.set_xlabel('neuron ID')
        ax.set_ylabel('neuron ID')

        im = ax.matshow(np.zeros((number_cells, number_cells)),vmin=-70.0, vmax=-45.0)

        if counter==0:
            plt.subplot(111)
            cbar = fig.colorbar(im)

        os.system("rm -r movies/"+labels[counter]+"/*")
        os.system("mkdir movies/"+labels[counter])

        l = [nd for nd in np.arange(0,len(recorders)) if (recorders[nd][1] == population and recorders[nd][2] == model)][0]
        data = nest.GetStatus(recorders[l][0],keys='events')
        pop = [nd for nd in nest.GetLeaves(population)[0] if nest.GetStatus([nd], 'model')[0]==model]

        senders = data[0]['senders']

        for t in np.arange(0.0,simtime-resolution,resolution):
            mult_pos = 0
            raster = np.zeros((number_cells, number_cells))
            for x in range(0, number_cells):
                for y in range(0, number_cells):
                    selected_senders = np.where(senders==pop[mult_pos])
                    ind_rec = (data[0]['V_m'])[selected_senders[0]]
                    raster[y,x] = ind_rec[int(t)]
                    mult_pos+=1

            im=ax.matshow(raster,vmin=-70.0, vmax=-45.0)
            plt.savefig('movies/'+labels[counter]+'/'+str(t)+'.png')

        counter+=1

def showMovie(label,simtime,resolution):

    img = None

    for t in np.arange(0.0,simtime-resolution,resolution):
        stdout.write("\r Time: %d ms" % t)
        stdout.flush()
#        sleep(0.001)
        im=pl.imread('movies/'+label+'/'+str(t)+'.png')
        if img is None:
            img = pl.imshow(im)
        else:
            img.set_data(im)
        pl.pause(.001)
        plt.axis("off")
        pl.draw()



