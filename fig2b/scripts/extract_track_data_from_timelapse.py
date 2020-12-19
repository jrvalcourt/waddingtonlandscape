from __future__ import division
import sys
from PIL import Image
import os
import numpy as np
import math
import pickle
#import matplotlib
#matplotlib.use('Qt5Agg')
#import matplotlib.pyplot as plt

filein = sys.argv[1]
tracks_dir = sys.argv[2]
img = Image.open(filein)

PIXEL_WIDTH  = 1
PIXEL_HEIGHT = 1
MINUTES_PER_TIMESTEP = 15
IMAGE_WIDTH = 708
IMAGE_HEIGHT = 711

tracks = {}
for traj in os.listdir(tracks_dir):
    if not traj.endswith('.csv'):
        continue
    key = '.'.join(traj.split('.')[:-1])
    with open(os.path.join(tracks_dir, traj)) as fin:
        fin.readline()
        tracks[key] = []
        for line in fin:
            e = line.split(',')
            tracks[key].append((int(e[8]), float(e[5]) / PIXEL_WIDTH, float(e[6]) / PIXEL_HEIGHT))

for track in tracks:
    tracks[track].sort(key=lambda x: x[0])

def get_channel_vals(t, x, y, oct4_channel, sox2_channel):
    x = int(round(x))
    y = int(round(y))
    if y < 2:
        y = 2
    if x < 2:
        x = 2
    if y > IMAGE_HEIGHT - 2:
        y = IMAGE_HEIGHT - 2
    if x > IMAGE_WIDTH - 2:
        x = IMAGE_WIDTH - 2
    return np.mean(oct4_channel[y-2:y+3,x-2:x+3].flatten()), np.mean(sox2_channel[y-2:y+3,x-2:x+3].flatten())

data = {}
for track in tracks:
    data[track] = []

print("Getting data from stack...")
print("    Starting...", end='')
for ii in np.linspace(1, img.n_frames / 3, num=int(img.n_frames / 3)):
    idx = int(ii)
    print("\r%.3f%%         " % (idx*3 / float(img.n_frames) * 100), end='')
    sys.stdout.flush()
    img.seek((idx - 1) * 3 + 1)
    oct4_channel = np.array(img)
    img.seek((idx - 1) * 3 + 2)
    sox2_channel = np.array(img)
    for track in tracks:
        for frame_num, x, y in tracks[track]:
            if idx == frame_num:
                oct4, sox2 = get_channel_vals(frame_num, x, y, oct4_channel, sox2_channel)
                #if track == '0009':
                #    if frame_num == 400:
                #        print((frame_num, x, y, oct4, sox2))
                #        #plt.subplot(121)
                #        #plt.imshow(oct4_channel)
                #        #plt.subplot(122)
                #        #plt.imshow(sox2_channel)
                #        #plt.show()
                data[track].append((frame_num, (frame_num - 1) * MINUTES_PER_TIMESTEP, x, y, oct4, sox2))
print('\rDone.         ')

pickle.dump(data, open("output/cell_data.pkl", 'wb'))
