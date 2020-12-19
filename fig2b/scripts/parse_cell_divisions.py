import sys
import os
import pickle

cell_div_dir = sys.argv[1]
data = {}

for f in os.listdir(cell_div_dir):
    if not f.endswith('.txt'):
        continue
    key = '.'.join(f.split('.')[:-1])
    filename = os.path.join(cell_div_dir, f)
    frame_nums = []
    with open(filename) as fin:
        for line in fin:
            frame_nums.append(int(line))
    data[key] = sorted(frame_nums)

pickle.dump(data, open("output/cell_div_frame_nums.pkl", "wb"))
