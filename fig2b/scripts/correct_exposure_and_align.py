from __future__ import division
import sys
import os
import numpy as np
from PIL import Image
from tifffile import imsave
from scipy import ndimage
from skimage.feature import register_translation
import cv2

# parameters for image correction
N_CHANNELS = 3
ALIGN_USING_CHANNEL = 2
NORM_TO = 1000
MIN_DELTA_TO_CORRECT = 1 
PERCENTILE_TO_NORM = 2

# read in the timelapse
timelapse = Image.open(sys.argv[1])
IMG_SIZE = np.array(timelapse).shape[0]

offsets_file = None
if len(sys.argv) > 4:
    offsets_file = sys.argv[4]

################################################################################
# read in all of the empty field images to 
# calculate the illumination profile
empty_field_images_dir = sys.argv[3]
empty_field_images = []
print('Calculating illumination profile... ', end='')
for d in os.listdir(empty_field_images_dir):
    img_path = os.path.join(empty_field_images_dir, d, 'MMStack_Pos0.ome.tif')
    empty_field_images.append(Image.open(img_path))

# load the images into a numpy array
empty_field_images_array = np.zeros((len(empty_field_images), 
                                    empty_field_images[0].n_frames, 
                                    IMG_SIZE, IMG_SIZE), 'int16')
for ii, img in enumerate(empty_field_images):
    for ch in range(img.n_frames):
        img.seek(ch)
        empty_field_images_array[ii, ch, :, :] = img

# average over all of the empty field to build the illumination profile
mean_empty_field = np.mean(empty_field_images_array, axis=0)
print('Done.')
################################################################################

################################################################################
if offsets_file is not None:
    max_x_offset = 0
    min_x_offset = 0
    max_y_offset = 0
    min_y_offset = 0
    offsets = []
    print('Reading offsets...', end='')
    with open(offsets_file) as fin:
        for line in fin:
            line= line.strip()
            if line[0] == '#':
                continue
            e = line.split(',')
            to_frame = int(e[0])
            x_offset = int(e[1])
            y_offset = int(e[2])
            if x_offset > max_x_offset:
                max_x_offset = x_offset
            if x_offset < min_x_offset:
                min_x_offset = x_offset
            if y_offset > max_y_offset:
                max_y_offset = y_offset
            if y_offset < min_y_offset:
                min_y_offset = y_offset
    with open(offsets_file) as fin:
        for line in fin:
            line= line.strip()
            if line[0] == '#':
                continue
            e = line.split(',')
            to_frame = int(e[0])
            x_offset = int(e[1])
            y_offset = int(e[2])
            if min_x_offset < 0:
                x_offset += abs(min_x_offset)
            if min_y_offset < 0:
                y_offset += abs(min_y_offset)
            offsets.append((to_frame, x_offset, y_offset))
    print('Done.')
################################################################################

################################################################################
image_out = None
if offsets_file is not None:
    image_out = np.zeros((timelapse.n_frames, 
                          IMG_SIZE + (max_x_offset - min_x_offset), 
                          IMG_SIZE + (max_y_offset - min_y_offset)), 'uint16')
else:
    image_out = np.zeros((timelapse.n_frames, 
                          IMG_SIZE, 
                          IMG_SIZE), 'uint16')

offset_idx = 0
x_offset = None
y_offset = None
print('Correcting frames... 0%', end='')
for ii in range(timelapse.n_frames):
    print('\rCorrecting frames... %.2f%%      ' % (ii / 
            timelapse.n_frames * 100), end='')
    sys.stdout.flush()
    channel = ii % N_CHANNELS
    t = ii // N_CHANNELS
    if offsets_file is not None:
        while offsets[offset_idx][0] <= t:
            offset_idx += 1
        _, x_offset, y_offset = offsets[offset_idx]
    else:
        x_offset = 0
        y_offset = 0
    timelapse.seek(ii)
    frame = np.array(timelapse)
    frame = np.multiply(frame, np.mean(np.mean(mean_empty_field[channel, :, :]))
            / mean_empty_field[channel, :, :])
    frame = frame * (NORM_TO / np.percentile(frame, 1))
    frame = frame - np.percentile(frame, 1)
    frame[frame < 0] = 0
    image_out[ii,x_offset:x_offset+IMG_SIZE,y_offset:y_offset+IMG_SIZE] = frame
print('\rCorrecting frames... saving results.        ', end='')
sys.stdout.flush()
imsave(sys.argv[2], image_out)
print('\rCorrecting frames... Done.                  ')
################################################################################
