import os
import numpy as np
import random
from matplotlib import pyplot as plt
import pickle


def save_new_image_dic(savepath, image_dic):
    '''If an image is already saved under that name, save under savepath+'_'+str(k) '''
    savedir, savename = (os.path.dirname(savepath), os.path.basename(savepath))
    # print savedir, savename
    try:
        os.makedirs(savedir)
    except:
        pass
    nfiles = sum([savename in file for file in os.listdir(savedir)])
    if nfiles == 0:
        savepathfull = savepath + '.pickle'
    else:
        savepathfull = savepath + '_'+str(nfiles+1) + '.pickle'
    with open(savepathfull, 'wb') as handle:
        pickle.dump(image_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)

def create_default_image(new_image_elements, Np):
    #new_image_elements = {'horizontal':{'size':3, 'number':1}, 'vertical':{'size':3, 'number':1}, 'cross':{'size':3, 'number':1}}

    output = {}
    image = np.zeros((Np,Np))
    for elemtype in new_image_elements.keys():
        for i in range(new_image_elements[elemtype]['number']):
            x,y = (random.randint(0,Np-1), random.randint(0,Np-1))
            size = new_image_elements[elemtype]['size']
            if elemtype == 'cross':
                size = sorted(size)
                draw_bar(image, (x, y), size[::1], Np)
                draw_bar(image, (x, y), size[::-1], Np)
            else:
                draw_bar(image, (x,y), size, Np)
            output[elemtype]=(x,y) #Eventually returns the position of one elem of each type

    # plt.imshow(image, interpolation='nearest')
    # plt.show()
    output['luminances']=image
    return output


def draw_bar(image, pos, size, Np):
    # Bar is centered on x, y
    x,y = pos
    s1, s2 = size
    for i in range(s1):
        for j in range(s2):
            image[(x - s2/2 +j)%Np, (y - s1/2+i)%Np]=1
    return

# #
# new_image_elements = {'horizontal': {'size': (6,2), 'number': 1}, 'vertical': {'size': (3,50), 'number': 1},
#                           'cross': {'size': (8,2), 'number': 1}}
#
# i = create_default_image(new_image_elements, 40)
#
# print(i)
#
# from matplotlib import pyplot as plt
# plt.imshow(i['luminances'], interpolation='nearest')
# print(i['vertical'])
# plt.show()

