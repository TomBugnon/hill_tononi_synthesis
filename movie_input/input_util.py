import os
import numpy as np
from PIL import Image


def create_movie_input(params):

    #dirname = '/home/kfujii2/newNEST2/ht_model_movie/' + params['movie'] + '/dst/'
    #dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/movie/' + params['movie'] + '/'
    dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/movie/' + params['movie'] + '_rate' + params['ret_rate'] + '/dst/'
    #dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/movie/' + params['movie'] + '_rate' + params['ret_rate'] + '/scrbl_xy/'
    #dirname = '/media/kfujii2/TOSHIBA EXT/experimental_data/movie/' + params['movie'] + '_rate' + params['ret_rate'] + '/scrbl_t/'
    num_files = len(os.listdir(dirname))

    input_img = np.zeros([num_files, params['Np']*params['Np']])
    spikes = np.zeros([num_files, params['Np']*params['Np']])


    # ==========================
    #  Load image data
    # ==========================
    for t in range(0, num_files, 1):

        # Load an image
        #image_filename = dirname + params['movie'] + str(t+1).zfill(4) + '.png'
        #image_filename = dirname + params['movie'] + '_' + str(t).zfill(2) + '.png'
        image_filename = dirname + params['movie'] + '_' + str(t+1).zfill(4) + '.png'

        # Convert from RGB to gray
        test_img_tmp = Image.open(image_filename).convert('L')

        # Convert PIL format to numpy format
        # Convert int to float
        test_img = np.array(test_img_tmp).astype(float)
        if len(test_img) != params['Np']:
            print ('Check the size of images, that should be equal to params[Np]')

        # Reshape to a row vector
        # Put data into 'input_img' which contains all time-serise image data
        input_img[t] = test_img.reshape(1, params['Np']*params['Np'])[0]

    # ==========================
    #  Calc intensity threshold
    # ==========================
    # calculate intensity_threshold based on params['ret_rate'](<= mean firing rate of the retina layer)
    '''
    tmp_intensity = np.sort(input_img, axis=1)
    tmp_threshold = np.zeros([num_files])
    for t in range(0, num_files, 1):
        tmp_threshold = tmp_intensity[t][params['Np']*params['Np']-params['ret_rate']]

    intensity_threshold = np.mean(tmp_threshold)
    '''
    intensity_threshold = 255*(4.0/5.0)

    # ==========================
    #  Set spike_times
    # ==========================
    for t in range(0, num_files, 1):

        # Set an image data to Retina nodes
        #intensity_threshold = np.sort(input_img)[Params['Np']*Params['Np']-Params['ret_rate']]
        #intensity_threshold = 255/3
        for i in range(0, params['Np']*params['Np'], 1):
            if input_img[t][i] > intensity_threshold:
                spikes[t][i] = 1


    return spikes