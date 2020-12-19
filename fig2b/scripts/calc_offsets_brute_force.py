from __future__ import division
import sys
import os
import numpy as np
from PIL import Image
from tifffile import imsave
from scipy import ndimage
from skimage.feature import register_translation
import cv2
import pickle

# parameters for image correction
N_CHANNELS = 3
ALIGN_USING_CHANNEL = 1
NORM_TO = 1000
MIN_DELTA_TO_CORRECT = 1 
PERCENTILE_TO_NORM = 2

if os.path.exists(sys.argv[2]):
    sys.exit()

# read in the timelapse
timelapse = Image.open(sys.argv[1])
IMG_SIZE = np.array(timelapse).shape[0]

def get_offset(img1, img2, r=100):
    min_x_offset = 0
    min_y_offset = 0
    minMSD = float('inf')
    for ii in np.linspace(-r, r, num=r*2+1):
        for jj in np.linspace(-r, r, num=r*2+1):
            x_offset = int(ii)
            y_offset = int(jj)
            diff = None
            if x_offset < 0 and y_offset < 0:
                diff = img1[0:x_offset,0:y_offset] - \
                         img2[abs(x_offset):,abs(y_offset):]
            elif x_offset == 0 and y_offset < 0:
                diff = img1[:,0:y_offset] - \
                         img2[:,abs(y_offset):]
            elif x_offset > 0 and y_offset < 0:
                diff = img1[x_offset:,0:y_offset] - \
                         img2[0:-x_offset,abs(y_offset):]

            elif x_offset < 0 and y_offset == 0:
                diff = img1[0:x_offset,:] - \
                         img2[abs(x_offset):,:]
            elif x_offset == 0 and y_offset == 0:
                diff = img1[:,:] - \
                         img2[:,:]
            elif x_offset > 0 and y_offset == 0:
                diff = img1[x_offset:,:] - \
                         img2[0:-x_offset,:]

            elif x_offset < 0 and y_offset > 0:
                diff = img1[0:x_offset,y_offset:] - \
                        img2[abs(x_offset):,0:-y_offset]
            elif x_offset == 0 and y_offset > 0:
                diff = img1[:,y_offset:] - \
                        img2[:,0:-y_offset]
            elif x_offset > 0 and y_offset > 0:
                diff = img1[x_offset:,y_offset:] - \
                        img2[0:-x_offset,0:-y_offset]
            MSD = np.mean(np.power(diff,2))
            if MSD < minMSD:
                minMSD = MSD
                min_x_offset = x_offset
                min_y_offset = y_offset

    return -min_x_offset, -min_y_offset

################################################################################

with open(sys.argv[2], 'w') as fout:
    running_x_offset = 0
    running_y_offset = 0
    print('Calculating offsets...', end='')
    previous_frame = None
    for ii in range(timelapse.n_frames):
        print('\rCalculating offsets... %.2f%%      ' % (ii / 
                timelapse.n_frames * 100), end='')
        sys.stdout.flush()
        channel = ii % N_CHANNELS
        if not channel == ALIGN_USING_CHANNEL:
            continue
        t = ii // N_CHANNELS
        timelapse.seek(ii)
        frame = np.array(timelapse)
        if previous_frame is not None:
            x_delta, y_delta = get_offset(previous_frame, frame)
            running_x_offset += x_delta
            running_y_offset += y_delta
            fout.write('%d,%d,%d\n' % (t+1, 
                                       -running_x_offset, 
                                       -running_y_offset))
        previous_frame = frame
    print('\rCalculating offsets... Done.             ')
################################################################################

